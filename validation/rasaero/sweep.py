"""Run the analyses of the RASAero II comparison over the full grid, through the same drivers the nTop blocks run.

    sh build/py.sh validation/rasaero/sweep.py flights|models|wave|parasite|vspaero NAME...

The blocks write a settings file and run ``python -m vsp_ntop.<driver>`` on it; here the settings files are the ones
the blocks wrote for the check runs of ``build_models.py``, with the flight condition, case and the setting under
study changed, so each run is exactly what the block would run. Grid: Mach 1.5, 2, 2.5, 3, 4, 5 (sea level, as the
RASAero II runs), alpha -15 to 15 degrees, beta -4 to 4 degrees. VSPAERO crosses every Mach number with every
Reynolds number of its flight condition, so it runs one flight condition per Mach number.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rocket as rk  # noqa: E402

MACHS = (1.5, 2.0, 2.5, 3.0, 4.0, 5.0)
ALPHAS = (-15.0, -12.0, -8.0, -4.0, -2.0, 0.0, 2.0, 4.0, 8.0, 12.0, 15.0)
BETAS = (-4.0, -2.0, 0.0, 2.0, 4.0)
MODELS = os.path.join(rk.RUNS, "models").replace("\\", "/")
WORK = os.path.join(rk.RUNS, "sweep").replace("\\", "/")
TEMPLATE_FLIGHT = MODELS + "/thin_exposed_check_flight.ntvsp.settings.txt"
TEMPLATE_VSPAERO = MODELS + "/thin_exposed_check_vspaero.ntvsp.settings.txt"
TEMPLATE_MODEL = MODELS + "/thin_exposed_model.ntvsp.settings.txt"
MODEL_SETS = {  # variant: (body set, fin set, source model settings)
    "thin_exposed": None,
    "thin_exposed_body": ("thick", "thin"),
    "thick_exposed": ("none", "thick"),
    "thick_body": ("thick", "thick"),
}


def settings(template, changes, path):
    text = open(template, encoding="utf-8").read()
    for key, value in changes.items():
        pat = re.compile(rf"^{re.escape(key)} = .*$", re.M)
        line = f"{key} = {value}"
        text = pat.sub(lambda _: line, text) if pat.search(text) else text + line + "\n"
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    return path


def driver(module, path, log=None):
    p = subprocess.run([rk.PYTHON, "-m", f"vsp_ntop.{module}", path], capture_output=True, text=True, env=rk.env(),
                       cwd=rk.ROOT)
    if p.returncode:
        raise RuntimeError(f"{module} {path}: {p.stderr.strip()[-500:]}")
    return p


def flight_base(machs, name):
    return f"{WORK}/{name}"


def flights():
    """The full grid, and one flight condition per Mach number (for VSPAERO)."""
    os.makedirs(WORK, exist_ok=True)
    jobs = [("flight_all", MACHS)] + [(f"flight_M{m:g}", (m,)) for m in MACHS]
    for name, machs in jobs:
        path = settings(TEMPLATE_FLIGHT, {"work_dir": WORK, "case": name, "mach": ";".join(f"{m:g}" for m in machs),
                                          "alpha": ";".join(f"{a:g}" for a in ALPHAS),
                                          "beta": ";".join(f"{b:g}" for b in BETAS)}, f"{WORK}/{name}.ntvsp.settings.txt")
        driver("run_flight", path)
        print("flight", name, flush=True)


def models(names):
    for name in names:
        body_set, fin_set = MODEL_SETS[name]
        text = open(TEMPLATE_MODEL, encoding="utf-8").read()
        text, n1 = re.subn(r"(kind: fuselage;[^\n]*?; set: )\w+", rf"\g<1>{body_set}", text)
        text, n2 = re.subn(r"(kind: wing;[^\n]*?; set: )\w+", rf"\g<1>{fin_set}", text)
        if (n1, n2) != (1, 1):
            raise RuntimeError(f"{name}: set lines not found ({n1}, {n2})")
        text = re.sub(r"^case = .*$", f"case = {name}_model", text, flags=re.M)
        path = f"{MODELS}/{name}_model.ntvsp.settings.txt"
        open(path, "w", encoding="utf-8").write(text)
        driver("run_model", path)
        print("model", name, flush=True)


def wave():
    path = settings(MODELS + "/thin_exposed_check_wave.ntvsp.settings.txt",
                    {"work_dir": WORK, "case": "wave", "flight": f"{WORK}/flight_all"}, f"{WORK}/wave.ntvsp.settings.txt")
    driver("run_wave", path)
    print("wave done", flush=True)


def parasite():
    path = settings(MODELS + "/thin_exposed_check_parasite.ntvsp.settings.txt",
                    {"work_dir": WORK, "case": "parasite", "flight": f"{WORK}/flight_all"},
                    f"{WORK}/parasite.ntvsp.settings.txt")
    driver("run_parasite", path)
    print("parasite done", flush=True)


RELAXED_MACHS = (1.5, 3.0, 5.0)  # relaxed wakes: compared with the fixed wake at these Mach numbers
WAKES = {  # name: (fixed_wake, wake_iterations)
    "fixed": (1, 5),
    "relaxed3": (0, 3),
    "relaxed5": (0, 5),
    "relaxed8": (0, 8),
}


def vspaero_case(model, wake, mach):
    return f"vsp_{model}_{wake}_M{mach:g}"


def vspaero(model, wake, ncpu=8, machs=None, parallel=2):
    if machs is None:  # SWEEP_MACHS=1.5,3,5 limits a study to some Mach numbers
        env = os.environ.get("SWEEP_MACHS", "")
        machs = tuple(float(x) for x in env.split(",")) if env else (MACHS if wake == "fixed" else RELAXED_MACHS)
    fixed, iters = WAKES[wake]

    def one(m):
        case = vspaero_case(model, wake, m)
        if os.path.exists(f"{WORK}/{case}.ntvsp.polar.csv"):
            return case
        path = settings(TEMPLATE_VSPAERO, {"work_dir": WORK, "case": case, "model": f"{MODELS}/{model}_model",
                                           "flight": f"{WORK}/flight_M{m:g}", "ncpu": ncpu, "fixed_wake": fixed,
                                           "wake_iterations": iters, "convergence_check": 0},
                        f"{WORK}/{case}.ntvsp.settings.txt")
        try:
            driver("run_vspaero", path)
            print("vspaero", case, flush=True)
        except RuntimeError as e:
            print("FAILED", case, str(e)[-300:], flush=True)
        return case

    with ThreadPoolExecutor(parallel) as pool:
        return list(pool.map(one, machs))


def main():
    what, names = sys.argv[1], sys.argv[2:]
    if what == "flights":
        flights()
    elif what == "models":
        models(names or [n for n, v in MODEL_SETS.items() if v])
    elif what == "wave":
        wave()
    elif what == "parasite":
        parasite()
    elif what == "probe":
        for n in names:
            probe(n)
    elif what == "vspaero":
        model, wakes = names[0], names[1:] or ["fixed"]
        for w in wakes:
            vspaero(model, w)



def probe(model, mach=2.0, alphas=(0.0, 4.0, 8.0), ncpu=16):
    """A few conditions only (the thick-body panel method stalls supersonic: every case runs to GMRES's cap)."""
    name = f"flight_probe_M{mach:g}"
    path = settings(TEMPLATE_FLIGHT, {"work_dir": WORK, "case": name, "mach": f"{mach:g}",
                                      "alpha": ";".join(f"{a:g}" for a in alphas), "beta": "0"},
                    f"{WORK}/{name}.ntvsp.settings.txt")
    driver("run_flight", path)
    case = f"probe_{model}_M{mach:g}"
    path = settings(TEMPLATE_VSPAERO, {"work_dir": WORK, "case": case, "model": f"{MODELS}/{model}_model",
                                       "flight": f"{WORK}/{name}", "ncpu": ncpu, "fixed_wake": 1,
                                       "convergence_check": 0}, f"{WORK}/{case}.ntvsp.settings.txt")
    driver("run_vspaero", path)
    print("probe", case, flush=True)


if __name__ == "__main__":
    main()
