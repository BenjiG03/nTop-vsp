"""Over what Mach range is OpenVSP's Wave Drag (supersonic area rule) a good estimate?

    python validation/wave_drag_mach.py          -> validation/results/wave_drag_mach.json, wave_drag_mach.png

OpenVSP's Wave Drag cuts the configuration with Mach planes at several roll angles, turns the cut areas into
equivalent bodies of revolution, and averages their von Karman slender-body wave drag: the supersonic area rule.
That is linear (small-disturbance) theory at zero lift, so its error has two parts, measured separately here:

1. OpenVSP against linear theory, which it implements:
   * a Sears-Haack body (fineness 10): linear theory gives D/q = 9 pi A_max^2 / (2 L^2), the same at every Mach;
   * rectangular wings of 4% biconvex sections, aspect ratio 12 and 3: Ackeret's theory gives
     CD = (16/3) (t/c)^2 / sqrt(M^2 - 1), which a fine roll-angle integration of the area rule reproduces;
     OpenVSP's result depends strongly on its number of rotation sections, scanned at Mach 2.
2. Linear theory against exact inviscid solutions, which shows where the small-disturbance assumption fails:
   * two-dimensional biconvex airfoils: shock-expansion theory (oblique shock at the leading edge, Prandtl-Meyer
     expansion along the surface), exact for an attached shock, against Ackeret;
   * cones: the Taylor-Maccoll solution against the slender-body cone pressure
     Cp = delta^2 (2 ln(2 / (delta sqrt(M^2 - 1))) - 1).

Run with a Python that has matplotlib (for the plot); OpenVSP 3.56 is found through the VSP_PATH environment
variable.
"""
from __future__ import annotations

import json
import math
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VSPSCRIPT = os.path.join(os.environ.get("VSP_PATH", os.path.join(ROOT, "OpenVSP-3.56.0-win64")), "vspscript.exe")
WORK = os.path.join(ROOT, "runs", "wave_study")
GAMMA = 1.4
MACHS = [1.05, 1.1, 1.2, 1.4, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]


# --- exact two-dimensional flow ----------------------------------------------------------------------------------------

def prandtl_meyer(m):
    g = GAMMA
    a = math.sqrt((g + 1) / (g - 1))
    return a * math.atan(math.sqrt((m * m - 1) / a / a)) - math.atan(math.sqrt(m * m - 1))


def inverse_pm(nu):
    lo, hi = 1.0, 100.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if prandtl_meyer(mid) < nu:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def oblique_shock(m, theta):
    """Weak shock angle for a deflection ``theta`` (rad), or None if detached."""
    g = GAMMA

    def f(beta):
        return math.tan(theta) - 2 / math.tan(beta) * (m * m * math.sin(beta) ** 2 - 1) / (
            m * m * (g + math.cos(2 * beta)) + 2)

    mu = math.asin(1 / m)
    # theta(beta) rises from 0 at the Mach angle to its maximum; find the maximum, then the weak root below it
    best, b_max = -1.0, mu
    for k in range(1, 2000):
        b = mu + (math.pi / 2 - mu) * k / 2000
        th = math.atan(2 / math.tan(b) * (m * m * math.sin(b) ** 2 - 1) / (m * m * (g + math.cos(2 * b)) + 2))
        if th > best:
            best, b_max = th, b
    if theta > best:
        return None
    lo, hi = mu, b_max
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def shock_relations(m, beta):
    g = GAMMA
    mn = m * math.sin(beta)
    p_ratio = 1 + 2 * g / (g + 1) * (mn * mn - 1)
    mn2 = math.sqrt((1 + (g - 1) / 2 * mn * mn) / (g * mn * mn - (g - 1) / 2))
    return p_ratio, mn2


def biconvex_exact(m, tc, n=400):
    """Wave drag coefficient of a 2D circular-arc biconvex airfoil at zero incidence, shock-expansion theory.
    None when the leading-edge shock detaches."""
    g = GAMMA
    # one surface: y = 2 tc x (1 - x) (chord 1), slope 2 tc (1 - 2x)
    theta0 = math.atan(2 * tc)
    beta = oblique_shock(m, theta0)
    if beta is None:
        return None
    p_ratio, mn2 = shock_relations(m, beta)
    m2 = mn2 / math.sin(beta - theta0)
    p0_2 = p_ratio * (1 + (g - 1) / 2 * m2 * m2) ** (g / (g - 1))  # total pressure behind the shock / p_inf
    nu2 = prandtl_meyer(m2)
    q = 0.5 * g * m * m
    cd = 0.0
    for k in range(n):
        x = (k + 0.5) / n
        theta = math.atan(2 * tc * (1 - 2 * x))
        mx = inverse_pm(nu2 + theta0 - theta)
        p = p0_2 / (1 + (g - 1) / 2 * mx * mx) ** (g / (g - 1))
        cd += (p - 1) / q * math.tan(theta) / n  # drag of one surface per unit chord
    return 2 * cd


def biconvex_ackeret(m, tc):
    return 16.0 / 3.0 * tc * tc / math.sqrt(m * m - 1)


# --- exact cone flow (Taylor-Maccoll) ----------------------------------------------------------------------------------

def taylor_maccoll(m, beta, steps=4000):
    """Integrate from the shock at ``beta`` towards the body; returns (cone half-angle, Cp on it)."""
    g = GAMMA
    p_ratio, mn2 = shock_relations(m, beta)
    delta = math.atan(2 / math.tan(beta) * (m * m * math.sin(beta) ** 2 - 1) / (m * m * (g + math.cos(2 * beta)) + 2))
    m2 = mn2 / math.sin(beta - delta)
    v2 = 1 / math.sqrt(2 / ((g - 1) * m2 * m2) + 1)  # V / V_max
    vr, vt = v2 * math.cos(beta - delta), -v2 * math.sin(beta - delta)

    def rhs(th, vr, vt):
        a = (g - 1) / 2 * (1 - vr * vr - vt * vt)
        dvt = (vt * vt * vr - a * (2 * vr + vt / math.tan(th))) / (a - vt * vt)
        return vt, dvt

    th, h = beta, -beta / steps
    while vt < 0 and th > 1e-6:
        k1 = rhs(th, vr, vt)
        k2 = rhs(th + h / 2, vr + h / 2 * k1[0], vt + h / 2 * k1[1])
        k3 = rhs(th + h / 2, vr + h / 2 * k2[0], vt + h / 2 * k2[1])
        k4 = rhs(th + h, vr + h * k3[0], vt + h * k3[1])
        nvr = vr + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        nvt = vt + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        if nvt >= 0:
            frac = -vt / (nvt - vt)
            th += h * frac
            vr = vr + frac * (nvr - vr)
            vt = 0.0
            break
        th, vr, vt = th + h, nvr, nvt
    # isentropic from behind the shock to the cone surface
    v2sq = v2 * v2
    p_shock = p_ratio
    p_cone = p_shock * ((1 - vr * vr) / (1 - v2sq)) ** (g / (g - 1))
    return th, (p_cone - 1) / (0.5 * g * m * m)


def cone_exact(m, delta):
    """Cone surface Cp for half-angle ``delta`` (rad), or None if the shock detaches."""
    mu = math.asin(1 / m)
    lo, hi = mu + 1e-6, math.pi / 2 - 1e-4
    d_hi = taylor_maccoll(m, hi)[0]
    best = max((taylor_maccoll(m, mu + (hi - mu) * k / 60)[0], mu + (hi - mu) * k / 60) for k in range(1, 60))
    if delta > best[0]:
        return None
    hi = best[1]
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if taylor_maccoll(m, mid, 1500)[0] < delta:
            lo = mid
        else:
            hi = mid
    return taylor_maccoll(m, 0.5 * (lo + hi))[1]


def cone_linear(m, delta):
    b = math.sqrt(m * m - 1)
    return delta * delta * (2 * math.log(2 / (delta * b)) - 1)


# --- OpenVSP ----------------------------------------------------------------------------------------------------------

def vsp(script_text, name):
    os.makedirs(WORK, exist_ok=True)
    path = os.path.join(WORK, name + ".vspscript")
    open(path, "w", newline="\n").write(script_text)
    p = subprocess.run([VSPSCRIPT, "-script", path], capture_output=True, text=True, cwd=WORK)
    out = p.stdout + p.stderr
    if "NTVSP_OK" not in out:
        raise RuntimeError(out[-2000:])
    return out


def wave_lines(sref, machs, slices, sections):
    lines = ['    string wc = FindContainer( "WaveDragSettings", 0 );',
             '    SetParmVal( wc, "RefFlag", "WaveDrag", 0 );',
             f'    SetParmVal( wc, "Sref", "WaveDrag", {sref!r} );',
             '    SetIntAnalysisInput( "WaveDrag", "Set", array<int> = {0} );',
             f'    SetIntAnalysisInput( "WaveDrag", "NumSlices", array<int> = {{{slices}}} );',
             f'    SetIntAnalysisInput( "WaveDrag", "NumRotSects", array<int> = {{{sections}}} );',
             '    SetIntAnalysisInput( "WaveDrag", "SymmFlag", array<int> = {1} );']
    for m in machs:
        lines += [f'    SetDoubleAnalysisInput( "WaveDrag", "Mach", array<double> = {{{m!r}}} );',
                  '    { string rid = ExecAnalysis( "WaveDrag" );',
                  f'      Print( "NTVSP_CD {m!r} " + GetDoubleResults( rid, "CDWave" )[0] ); }}']
    return lines


def sears_haack(machs, fineness=10.0, length=10.0, slices=40, sections=12, n=41):
    r = length / fineness / 2
    lines = ["void main()", "{", '    string g = AddGeom( "FUSELAGE", "" );', f'    SetParmVal( g, "Length", "Design", {length!r} );',
             "    string xs = GetXSecSurf( g, 0 );",
             f"    while ( GetNumXSec( xs ) < {n} ) {{ InsertXSec( g, GetNumXSec( xs ) - 2, XS_CIRCLE ); }}",
             f"    for ( int i = 1; i < {n - 1}; i++ ) {{ ChangeXSecShape( xs, i, XS_CIRCLE ); }}",
             '    SetParmVal( g, "Tess_W", "Shape", 33 );']
    xs = [0.5 * (1 - math.cos(math.pi * i / (n - 1))) for i in range(n)]
    for order in (range(1, n - 1), range(n - 2, 0, -1)):
        for i in order:
            lines.append(f'    SetParmVal( GetXSecParm( GetXSec( xs, {i} ), "XLocPercent" ), {xs[i]!r} ); Update();')
    for i in range(1, n - 1):
        ri = r * (4 * xs[i] * (1 - xs[i])) ** 0.75
        lines.append(f'    SetParmVal( GetXSecParm( GetXSec( xs, {i} ), "Circle_Diameter" ), {2 * ri!r} );')
        lines.append(f'    SetParmVal( GetXSecParm( GetXSec( xs, {i} ), "SectTess_U" ), 4 );')
    lines += ["    Update();"]
    a_max = math.pi * r * r
    lines += wave_lines(a_max, machs, slices, sections) + ['    Print( "NTVSP_OK" );', "}"]
    out = vsp("\n".join(lines) + "\n", f"sears_haack_{slices}_{sections}")
    cd = {float(l.split()[1]): float(l.split()[2]) for l in out.splitlines() if l.strip().startswith("NTVSP_CD")}
    theory = 9 * math.pi * a_max / (2 * length ** 2)  # D/q over A_max
    return cd, theory


def rect_wing(machs, tc=0.04, aspect=12.0, chord=1.0, slices=60, sections=12):
    span = aspect * chord
    lines = ["void main()", "{", '    string w = AddGeom( "WING", "" );',
             "    SetDriverGroup( w, 1, SPAN_WSECT_DRIVER, ROOTC_WSECT_DRIVER, TIPC_WSECT_DRIVER );",
             f'    SetParmVal( w, "Span", "XSec_1", {span / 2!r} );', f'    SetParmVal( w, "Root_Chord", "XSec_1", {chord!r} );',
             f'    SetParmVal( w, "Tip_Chord", "XSec_1", {chord!r} );', '    SetParmVal( w, "Sweep", "XSec_1", 0.0 );',
             '    SetParmVal( w, "Tess_W", "Shape", 41 );', '    SetParmVal( w, "SectTess_U", "XSec_1", 30 );',
             "    string xs = GetXSecSurf( w, 0 );",
             "    ChangeXSecShape( xs, 0, XS_BICONVEX ); ChangeXSecShape( xs, 1, XS_BICONVEX );",
             f'    SetParmVal( GetXSecParm( GetXSec( xs, 0 ), "ThickChord" ), {tc!r} );',
             f'    SetParmVal( GetXSecParm( GetXSec( xs, 1 ), "ThickChord" ), {tc!r} );', "    Update();"]
    lines += wave_lines(span * chord, machs, slices, sections) + ['    Print( "NTVSP_OK" );', "}"]
    out = vsp("\n".join(lines) + "\n", f"rect_wing_{slices}_{sections}")
    return {float(l.split()[1]): float(l.split()[2]) for l in out.splitlines() if l.strip().startswith("NTVSP_CD")}


def rotation_scan(mach=2.0, sections=(10, 20, 30, 60, 120)):
    out = {}
    for ar in (12.0, 3.0):
        out[str(ar)] = [(n, rect_wing([mach], aspect=ar, slices=60, sections=n)[mach]) for n in sections]
    return {"mach": mach, "ackeret": biconvex_ackeret(mach, 0.04), "by_aspect_ratio": out}


def plot(r, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    m = r["machs"]
    th = r["sears_haack"]["theory"]
    ax[0].plot(m, [v / th for v in r["sears_haack"]["openvsp_40x12"]], "o-", label="40 slices")
    ax[0].plot(m, [v / th for v in r["sears_haack"]["openvsp_20x10"]], "s--", label="20 slices")
    ax[0].axhline(1.0, color="k", lw=0.8)
    ax[0].set(xlabel="Mach", ylabel="OpenVSP / linear theory", title="Sears-Haack body, fineness 10")
    ax[0].legend()
    rs = r["rotation_scan"]
    for ar, rows in rs["by_aspect_ratio"].items():
        ax[1].plot([n for n, _ in rows], [v / rs["ackeret"] for _, v in rows], "o-", label=f"aspect ratio {float(ar):g}")
    ax[1].axhline(1.0, color="k", lw=0.8)
    ax[1].set(xscale="log", xlabel="rotation sections", ylabel="OpenVSP / Ackeret",
              title=f"4% biconvex rectangular wing, Mach {rs['mach']:g}")
    ax[1].legend()
    for tc, rows in r["biconvex"].items():
        ax[2].plot([mm * float(tc) for mm, a, e in rows if e], [a / e for mm, a, e in rows if e], "o-", label=f"biconvex t/c {tc}")
    for deg, rows in r["cone"].items():
        d = math.radians(float(deg))
        ax[2].plot([mm * d for mm, a, e in rows if e], [a / e for mm, a, e in rows if e], "s--", label=f"cone {deg} deg")
    ax[2].axhline(1.0, color="k", lw=0.8)
    ax[2].axhline(0.9, color="grey", lw=0.6, ls=":")
    ax[2].set(xlabel="Mach x thickness (t/c, or cone half-angle in rad)", ylabel="linear theory / exact",
              title="Linear theory against exact inviscid flow", ylim=(0, 1.1))
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=110)


def main():
    results = {"machs": MACHS}
    # 1. OpenVSP vs linear theory
    sh, sh_theory = sears_haack(MACHS)
    sh_coarse, _ = sears_haack(MACHS, slices=20, sections=10)
    results["sears_haack"] = {"theory": sh_theory, "openvsp_40x12": [sh.get(m) for m in MACHS],
                              "openvsp_20x10": [sh_coarse.get(m) for m in MACHS]}
    wing = rect_wing(MACHS)
    results["rect_wing"] = {"tc": 0.04, "aspect": 12.0, "openvsp": [wing.get(m) for m in MACHS],
                            "ackeret": [biconvex_ackeret(m, 0.04) for m in MACHS],
                            "exact_2d": [biconvex_exact(m, 0.04) for m in MACHS]}
    results["rotation_scan"] = rotation_scan()
    # 2. linear vs exact
    results["biconvex"] = {}
    for tc in (0.02, 0.04, 0.06, 0.10):
        results["biconvex"][str(tc)] = [(m, biconvex_ackeret(m, tc), biconvex_exact(m, tc)) for m in MACHS]
    results["cone"] = {}
    for deg in (5.0, 7.5, 10.0, 15.0):
        d = math.radians(deg)
        results["cone"][str(deg)] = [(m, cone_linear(m, d), cone_exact(m, d)) for m in MACHS]
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(results, open(os.path.join(HERE, "results", "wave_drag_mach.json"), "w"), indent=1)
    report(results)
    plot(results, os.path.join(HERE, "results", "wave_drag_mach.png"))


def report(r):
    print("Sears-Haack, fineness 10: linear theory CD (on A_max) =", round(r["sears_haack"]["theory"], 5))
    for m, a, b in zip(r["machs"], r["sears_haack"]["openvsp_40x12"], r["sears_haack"]["openvsp_20x10"]):
        print(f"  M {m:4.2f}: OpenVSP 40x12 {a:.5f} ({100 * (a / r['sears_haack']['theory'] - 1):+.1f}%)  20x10 {b:.5f}")
    print("Rectangular biconvex wing, t/c 4%, AR 12:")
    for m, v, ack, ex in zip(r["machs"], r["rect_wing"]["openvsp"], r["rect_wing"]["ackeret"], r["rect_wing"]["exact_2d"]):
        print(f"  M {m:4.2f}: OpenVSP {v:.5f}  Ackeret 2D {ack:.5f}  exact 2D {ex if ex is None else round(ex, 5)}")
    rs = r["rotation_scan"]
    for ar, rows in rs["by_aspect_ratio"].items():
        print(f"rectangular wing AR {ar}, Mach {rs['mach']}: OpenVSP / Ackeret by rotation sections:",
              " ".join(f"{n}:{v / rs['ackeret']:.2f}" for n, v in rows))
    for tc, rows in r["biconvex"].items():
        print(f"biconvex t/c {tc}: Ackeret / exact - 1 =",
              " ".join(f"M{m:g}:{'detached' if e is None else f'{100 * (a / e - 1):+.0f}%'}" for m, a, e in rows))
    for deg, rows in r["cone"].items():
        print(f"cone {deg} deg: linear / exact - 1 =",
              " ".join(f"M{m:g}:{'detached' if e is None else f'{100 * (a / e - 1):+.0f}%'}" for m, a, e in rows))


if __name__ == "__main__":
    main()
