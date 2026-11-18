"""Score combinations of nTop_vsp analyses against RASAero II on the hobby rocket.

    sh build/py.sh validation/rasaero/compare.py [--out results.json]

Every source enters as the Aero Database does it (vsp_ntop.run_database.Contribution: VSPAERO polars, Wave and
Parasite Drag tables, Modified Newtonian coefficient tables, with Use = all / drag / nodrag), on the sweep's flight
condition (Mach 1.5 to 5 at sea level, alpha -15 to 15, beta -4 to 4; Sref = pi d^2 / 4, every moment over Sref d,
about x = 0.7 m on the axis). RASAero II's values at each condition come from rasaero_ref (axisymmetric: normal
force along the crossflow at the total angle of attack, centre of pressure on the axis).

Newtonian runs with Cp_max = 1 (per unit Cp_max); a Cp_max law multiplies them: Modified (Rayleigh pitot), Classical
(2), or "fit": the constant per component and Mach number that best matches RASAero II (least squares on the
normalised errors below), which is what a User Defined Cp Max Field holding those constants would give.

Errors at each Mach number, over its 55 (alpha, beta) conditions: for CN, CY, Cm, Cn the RMS error over the RMS
of RASAero's values; for CA the RMS error over RASAero's mean CA. The score is the mean of the five (Cl is
reported, not scored: RASAero II gives none).
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rocket as rk  # noqa: E402
import sweep  # noqa: E402
from rasaero_ref import RASAero  # noqa: E402

sys.path.insert(0, rk.ROOT)
from vsp_ntop import run_database as db  # noqa: E402
from vsp_ntop.flight import Flight  # noqa: E402
from vsp_ntop.run_newtonian_mesh import cp_max_stagnation  # noqa: E402

KEYS = ("CA", "CY", "CN", "Cl", "Cm", "Cn")
SCORED = ("CA", "CY", "CN", "Cm", "Cn")
X_REF, D = 0.7, 0.054


def quiet(*_):
    pass


class Data:
    """Every source's six coefficients at every condition of the sweep, as arrays (n_conditions, 6)."""

    def __init__(self):
        self.flight = Flight(f"{sweep.WORK}/flight_all")
        self.conds = self.flight.conditions
        self.ref = self.flight.ref
        self.mach = np.array([c["Mach"] for c in self.conds])
        r = RASAero()
        self.rasaero = np.array([[r.coefficients(c["Mach"], c["AoA"], c["Beta"], X_REF, D)[k] for k in KEYS]
                                 for c in self.conds])
        self.base_drag = np.array([r.breakdown(c["Mach"]).get("base", 0.0) for c in self.conds])
        self.cache = {}

    def _eval(self, contribs):
        out = np.zeros((len(self.conds), 6))
        for i, c in enumerate(self.conds):
            for con in contribs:
                if not con.applies(c["Mach"]):
                    continue
                got = con.on(self.ref, c)
                if got is not None:
                    out[i] = np.add(out[i], list(got[0]) + list(got[1]))
        return out

    def source(self, kind, results, use="all"):
        """A source's coefficients; ``results`` a base or a list of per-Mach bases (VSPAERO)."""
        key = (kind, tuple(results) if isinstance(results, list) else results, use)
        if key not in self.cache:
            bases = results if isinstance(results, list) else [results]
            missing = [b for b in bases if not os.path.exists(b + {"vspaero": ".ntvsp.polar.csv", "wave": ".ntvsp.wave.csv",
                                                                       "parasite": ".ntvsp.parasite.csv",
                                                                       "coefficients": ".ntvsp.coefficients.csv"}[kind])]
            if missing:
                self.cache[key] = None
            else:
                cons = [db.Contribution({"kind": kind, "results": b, "use": use, "name": os.path.basename(b)}, quiet)
                        for b in bases]
                self.cache[key] = self._eval(cons)
        return self.cache[key]


def cp_law(name, mach):
    if name == "modified":
        return np.array([cp_max_stagnation(m) for m in mach])
    if name == "classical":
        return np.full(len(mach), 2.0)
    raise ValueError(name)


def target(data: Data, base=True):
    """RASAero II's coefficients; without its base drag (taken off CA) when ``base`` is False."""
    if base:
        return data.rasaero
    t = data.rasaero.copy()
    t[:, 0] -= data.base_drag
    return t


MACHS = None  # Mach numbers scored (None: all)


def errors(data: Data, total, base=True):
    """Per Mach: normalised errors per coefficient and the score."""
    out = {}
    tgt = target(data, base)
    for m in sorted(set(data.mach)):
        if MACHS and not any(abs(m - x) < 1e-9 for x in MACHS):
            continue
        sel = data.mach == m
        ref, got = tgt[sel], total[sel]
        e = {}
        for j, k in enumerate(KEYS):
            rms = math.sqrt(np.mean((got[:, j] - ref[:, j]) ** 2))
            if k == "CA":
                e[k] = rms / abs(np.mean(ref[:, j]))
            elif k == "Cl":
                e[k] = rms  # absolute
            else:
                e[k] = rms / math.sqrt(np.mean(ref[:, j] ** 2))
        e["score"] = float(np.mean([e[k] for k in SCORED]))
        out[float(m)] = e
    return out


def fit_scales(data: Data, fixed, units, base=True):
    """Per Mach, the scales k_i >= 0 of ``units`` (arrays) that minimise the normalised squared errors of
    fixed + sum k_i units against RASAero II. Returns (total, {mach: [k_i]})."""
    total = fixed.copy()
    ks = {}
    tgt = target(data, base)
    for m in sorted(set(data.mach)):
        sel = data.mach == m
        ref = tgt[sel]
        w = []
        for j, k in enumerate(KEYS):
            if k not in SCORED:
                w.append(0.0)
            elif k == "CA":
                w.append(1.0 / abs(np.mean(ref[:, j])))
            else:
                w.append(1.0 / math.sqrt(np.mean(ref[:, j] ** 2)))
        w = np.array(w)
        a = np.column_stack([(u[sel] * w).ravel() for u in units])
        b = ((ref - fixed[sel]) * w).ravel()
        if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
            total[sel] = np.nan  # a source failed at this Mach number
            ks[float(m)] = [float("nan")] * len(units)
            continue
        k_, *_ = np.linalg.lstsq(a, b, rcond=None)
        k_ = np.maximum(k_, 0.0)
        for kk, u in zip(k_, units):
            total[sel] += kk * u[sel]
        ks[float(m)] = [float(x) for x in k_]
    return total, ks


def vsp_bases(model, wake):
    """Per-Mach VSPAERO results. ``thin_extend_whole`` is the Exposed Strip run's whole-wing polar (the driver keeps
    it as ``.ntvsp.polar_whole.csv``), copied under its own case name."""
    if model == "thin_extend_whole":
        import shutil
        out = []
        for m in sweep.MACHS:
            src = f"{sweep.WORK}/{sweep.vspaero_case('thin_extend_exposed', wake, m)}"
            dst = f"{sweep.WORK}/{sweep.vspaero_case(model, wake, m)}"
            if os.path.exists(src + ".ntvsp.polar_whole.csv") and not os.path.exists(dst + ".ntvsp.polar.csv"):
                shutil.copy(src + ".ntvsp.polar_whole.csv", dst + ".ntvsp.polar.csv")
                shutil.copy(src + ".ntvsp.settings.txt", dst + ".ntvsp.settings.txt")
            out.append(dst)
        return out
    return [f"{sweep.WORK}/{sweep.vspaero_case(model, wake, m)}" for m in sweep.MACHS]


def candidates(data: Data, base=True):
    """(name, parts, total array) for every combination whose data are all present.

    A combination is: a lift source (VSPAERO with Use = All or No drag, or none), a Newtonian body term (none, Use All
    or No drag), a pressure-drag term (OpenVSP Wave Drag, Newtonian fin drag (Use Drag), Newtonian whole-rocket drag
    (Use Drag), or none), and a friction term (Parasite Drag with either compressible Cf equation, or none). Every
    Newtonian term takes one Cp_max law (Modified or Classical) or is fitted (a constant per term and Mach number).
    The all-Newtonian rocket (whole, or body + fins, Use All) is a lift source of its own.
    """
    W = sweep.WORK
    zero = np.zeros_like(data.rasaero)
    newt = {(c, u): data.source("coefficients", f"{W}/newt_{c}_unit", u) for c in ("whole", "body", "fins")
            for u in ("all", "nodrag", "drag")}
    wave = data.source("wave", f"{W}/wave", "drag")
    frictions = {"Parasite (Schlichting Compressible)": data.source("parasite", f"{W}/parasite", "drag"),
                 "Parasite (Roughness Schlichting Avg Compressible)": data.source("parasite", f"{W}/parasite_cf19", "drag"),
                 "no friction": zero}
    lifts = {"no VSPAERO": zero}
    for model in ("thin_exposed", "thin_axis", "thin_exposed_body", "thick_exposed", "thick_body", "thin_extend_exposed", "thin_extend_carry",
                  "thin_extend_whole"):
        for wake in sweep.WAKES:
            for use in ("all", "nodrag"):
                got = data.source("vspaero", vsp_bases(model, wake), use)
                if got is not None:
                    lifts[f"VSPAERO {model} {wake} (Use {'All' if use == 'all' else 'No drag'})"] = got
    bodies = {"": None, "Newtonian body (All)": ("body", "all"), "Newtonian body (No drag)": ("body", "nodrag")}
    body_lift = {"Body Lift (slender body)": data.source("coefficients", f"{W}/bodylift_pot"),
                 "Body Lift (slender body + crossflow)": data.source("coefficients", f"{W}/bodylift")}
    drags = {"wave": ("fixed", wave), "Newtonian fin drag": ("newt", ("fins", "drag")),
             "Newtonian whole drag": ("newt", ("whole", "drag")), "no pressure drag": ("fixed", zero)}
    thick_fin_drag = data.source("vspaero", vsp_bases("thick_exposed", "fixed"), "drag")
    if thick_fin_drag is not None:
        drags["VSPAERO thick fin drag"] = ("fixed", thick_fin_drag)
        drags["VSPAERO thick fin drag + Newtonian body drag"] = ("newt+", (thick_fin_drag, ("body", "drag")))
    out = []
    laws = ("modified", "classical", "fitted")

    def combine(fixed, units, law):
        if any(u is None for u in units) or fixed is None:
            return None, None
        if not units:
            return fixed, None
        if law == "fitted":
            return fit_scales(data, fixed, units, base)
        k = cp_law(law, data.mach)[:, None]
        return fixed + k * sum(units), None

    for fname, fr in frictions.items():
        if fr is None:
            continue
        # all-Newtonian rockets
        for law in laws:
            for label, units in (("Newtonian whole rocket (All)", [newt[("whole", "all")]]),
                                 ("Newtonian body + fins (All)", [newt[("body", "all")], newt[("fins", "all")]])):
                tot, ks = combine(fr, units, law)
                if tot is not None:
                    out.append((f"{label} [{law} Cp_max] + {fname}", {"cp_max_fit": ks} if ks else {}, tot))
        for lname, lift in lifts.items():
            # Body Lift (no Cp_max) with any drag source
            for blname, bl in body_lift.items():
                if bl is None or lname == "no VSPAERO":
                    continue
                for dname, (dkind, dval) in drags.items():
                    if dkind == "fixed":
                        if dval is not None:
                            out.append((f"{lname} + {blname} + {dname} + {fname}", {}, lift + bl + dval + fr))
                    elif dkind == "newt":
                        for law in laws:
                            tot, ks = combine(lift + bl + fr, [newt[dval]], law)
                            if tot is not None:
                                out.append((f"{lname} + {blname} + {dname} + {fname} [{law} Cp_max]",
                                            {"cp_max_fit": ks} if ks else {}, tot))
            for bname, bkey in bodies.items():
                if lname == "no VSPAERO" and not bkey:
                    continue
                for dname, (dkind, dval) in drags.items():
                    if dkind == "newt+":  # a fixed drag plus a Newtonian one
                        if bkey == ("body", "all"):
                            continue  # the body's drag is in its All term already
                        units = ([newt[bkey]] if bkey else []) + [newt[dval[1]]]
                        fixed = lift + fr + dval[0]
                    else:
                        units = ([newt[bkey]] if bkey else []) + ([newt[dval]] if dkind == "newt" else [])
                        fixed = lift + fr + (dval if dkind == "fixed" else 0.0) if (dkind != "fixed" or dval is not None) \
                            else None
                    for law in (laws if units else ("-",)):
                        tot, ks = combine(fixed, units, law)
                        if tot is None:
                            continue
                        parts = [lname] + ([bname] if bname else []) + [dname, fname]
                        name = " + ".join(parts) + (f" [{law} Cp_max]" if units else "")
                        out.append((name, {"cp_max_fit": ks} if ks else {}, tot))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default=os.path.join(HERE, "results", "scores.json"))
    p.add_argument("--top", type=int, default=25)
    p.add_argument("--no-base", action="store_true", help="score CA against RASAero's CA less its base drag")
    p.add_argument("--machs", default="", help="score only these Mach numbers, e.g. 2,2.5,3,4,5")
    args = p.parse_args()
    global MACHS
    if args.machs:
        MACHS = [float(x) for x in args.machs.split(",")]
        args.out = args.out.replace(".json", "_M" + "_".join(f"{m:g}" for m in MACHS) + ".json")
    base = not args.no_base
    if args.no_base:
        args.out = args.out.replace("scores", "scores_no_base", 1)
    data = Data()
    cands = candidates(data, base)
    rows = []
    for name, parts, total in cands:
        e = errors(data, total, base)
        mean = {k: float(np.mean([v[k] for v in e.values()])) for k in list(KEYS) + ["score"]}
        if not np.isfinite(mean["score"]):
            continue
        rows.append({"name": name, "parts": parts, "per_mach": e, "mean": mean})
    rows.sort(key=lambda r: r["mean"]["score"])
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    json.dump(rows, open(args.out, "w"), indent=1)
    print(f"{len(rows)} combinations; RASAero base drag (not modelled by any source) is "
          + ", ".join(f"M{m:g} {data.base_drag[data.mach == m][0]:.3f}" for m in sorted(set(data.mach))))
    print(f"{'score':>6} {'CA':>6} {'CN':>6} {'CY':>6} {'Cm':>6} {'Cn':>6} {'Cl':>7}  combination")
    for r in rows[:args.top]:
        m = r["mean"]
        print(f"{m['score']:6.3f} {m['CA']:6.3f} {m['CN']:6.3f} {m['CY']:6.3f} {m['Cm']:6.3f} {m['Cn']:6.3f} {m['Cl']:7.4f}  "
              f"{r['name']}")
    print("\nbest per Mach, Cp_max from a law (no fitting):")
    report_per_mach(rows)
    print("\nbest per Mach, fitted Cp_max allowed:")
    report_per_mach(rows, fitted=True)



def report_per_mach(rows, top=3, fitted=False):
    """The best combinations at each Mach number (without fitted Cp_max unless ``fitted``)."""
    pool = [r for r in rows if fitted or "fitted" not in r["name"]]
    for mach in sorted(pool[0]["per_mach"], key=float):
        best = sorted(pool, key=lambda r: r["per_mach"][mach]["score"])[:top]
        print(f"  M{float(mach):g}:")
        for b in best:
            e = b["per_mach"][mach]
            print(f"     {e['score']:.3f} (CA {e['CA']:.3f} CN {e['CN']:.3f} Cm {e['Cm']:.3f})  {b['name']}")


if __name__ == "__main__":
    main()
