"""Dev benchmark: gray box hessian options on the softening-acid-RO flowsheet
from reaktoro_enabled_watertap (stability sweep style, warm started recovery sweep).

Not part of the package. Each case runs in its own process so hessian defaults can be
overridden before reaktoro_pse is imported, and failures/hangs are isolated.

Usage:
    python dev/bench_softening_ro.py --out results.csv [--cases a,b] [--waters BGW,SW_HPRO]
        [--rec-workers 4] [--retw PATH_TO_reaktoro_enabled_watertap/src]
    python dev/bench_softening_ro.py --run-case JSON   (internal)
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

DEFAULT_RETW = r"D:\github\reaktoro_enabled_watertap\src"

# name -> (rkt_hessian_type, bfgs_initialization_type, extra HESSIAN_DEFAULTS overrides)
CASES = {
    "lmt_sc1": ("limited-memory", "scalar1", {}),
    "lbfgs_gn": ("LBFGS", "GaussNewton", {}),
    "lbfgs_sc1": ("LBFGS", "scalar1", {}),
    "bfgs_ipopt_gn": ("BFGS_ipopt", "GaussNewton", {}),
    "bfgs_ipopt_sc1": ("BFGS_ipopt", "scalar1", {}),
    "sr1_zero": ("SR1", "zero", {}),
    "sr1_gn": ("SR1", "GaussNewton", {}),
    "sr1_sc1": ("SR1", "scalar1", {}),
    "psb_zero": ("PSB", "zero", {}),
    "psb_gn": ("PSB", "GaussNewton", {}),
    "psb_sc1": ("PSB", "scalar1", {}),
    "sr1lag_zero": ("SR1_lagrangian", "zero", {}),
    "sr1lag_sc1": ("SR1_lagrangian", "scalar1", {}),
    # scaled variants
    "lbfgs_sc1_scaled": ("LBFGS", "scalar1", {"hessian_input_scaling": True}),
    "sr1_zero_scaled": ("SR1", "zero", {"hessian_input_scaling": True}),
    "sr1_sc1_scaled": ("SR1", "scalar1", {"hessian_input_scaling": True}),
    "psb_sc1_scaled": ("PSB", "scalar1", {"hessian_input_scaling": True}),
    "sr1lag_sc1_scaled": ("SR1_lagrangian", "scalar1", {"hessian_input_scaling": True}),
}

WATERS = {
    "BGW": ("USDA_brackish.yaml", False),
    "SW_RO": ("Seawater.yaml", False),
    "SW_HPRO": ("Seawater.yaml", True),
    "BGW_500": ("sample_500_hardness.yaml", False),
    "BGW_1500": ("sample_1500_hardness.yaml", False),
}

RECOVERIES = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]


def run_case(spec):
    """runs inside subprocess, prints one JSON line per solve prefixed with RESULT"""
    import reaktoro_pse
    from reaktoro_pse.core.util_classes import hessian_functions

    if spec["overrides"]:
        # must happen before modules that read the defaults at import are imported
        hessian_functions.HESSIAN_DEFAULTS.update(spec["overrides"])
    print(f"reaktoro_pse from {reaktoro_pse.__file__}", file=sys.stderr)

    import reaktoro_enabled_watertap.flowsheets.softening_acid_ro.softening_acid_ro as sar
    from pyomo.environ import value

    water_file, hpro = WATERS[spec["water"]]

    def emit(stage, recovery, ok, t, m=None, err=None):
        rec = {
            "tag": spec["tag"],
            "case": spec["case"],
            "water": spec["water"],
            "stage": stage,
            "recovery": recovery,
            "ok": ok,
            "time": round(t, 2),
            "iters": None,
            "lcow": None,
            "error": err,
        }
        if m is not None and ok:
            rec["iters"] = int(value(m.fs.ipopt_iterations["Number of iterations"]))
            try:
                rec["lcow"] = value(m.fs.costing.LCOW)
            except Exception:
                pass
        print("RESULT " + json.dumps(rec), flush=True)

    hess_type, init_type, _ = CASES[spec["case"]]
    m = None
    try:
        t0 = time.time()
        m = sar.build_model(
            water_file,
            multi_process_reaktoro=True,
            hpro=hpro,
            rkt_hessian_type=hess_type,
            bfgs_initialization_type=init_type,
            softening_reagents=["CaO", "Na2CO3"],
            acidification_reagents=["HCl", "H2SO4"],
            system_costing="Amusat_et_al_2024",
        )
        sar.initialize(m, linear_solver=spec["linear_solver"], tee=False)
        emit("init", None, True, time.time() - t0, m)
    except Exception as e:
        emit("init", None, False, time.time() - t0, err=repr(e)[:300])
        if m is not None and m.find_component("reaktoro_manager") is not None:
            m.reaktoro_manager.terminate_workers()
        return
    for r in spec["recoveries"]:
        if not feasible_recovery(spec["water"], r):
            continue
        t0 = time.time()
        try:
            m.fs.water_recovery.fix(r)
            sar.solve_model(m, linear_solver=spec["linear_solver"], tee=False)
            emit("solve", r, True, time.time() - t0, m)
        except Exception as e:
            emit("solve", r, False, time.time() - t0, err=repr(e)[:300])
    m.reaktoro_manager.terminate_workers()


def feasible_recovery(water, r):
    """same limits as sar.test_func: seawater max 64% recovery with RO, 86% with HPRO"""
    if water == "SW_RO":
        return r <= 0.64
    if water == "SW_HPRO":
        return r <= 0.86
    return True


def launch(spec, timeout):
    t0 = time.time()
    cmd = [sys.executable, os.path.abspath(__file__), "--run-case", json.dumps(spec)]
    # PYTHONPATH (not sys.path) so reaktoro worker processes see the same code
    env = dict(os.environ)
    paths = [p for p in [spec["rktpse_src"], spec["retw"]] if p]
    env["PYTHONPATH"] = os.pathsep.join(paths + [env.get("PYTHONPATH", "")])
    # Ipopt loads HSL (libhsl.dll) by name at runtime; without conda activation the env's
    # Library/bin is not on PATH and MA27 fails with "Library loading failure"
    lib_bin = os.path.join(os.path.dirname(sys.executable), "Library", "bin")
    if os.path.isdir(lib_bin):
        env["PATH"] = os.pathsep.join([lib_bin, env.get("PATH", "")])
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=os.getcwd(),
            env=env,
        )
        out = proc.stdout
        timed_out = False
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        timed_out = True
    rows = [json.loads(l[7:]) for l in out.splitlines() if l.startswith("RESULT ")]
    if timed_out:
        rows.append(
            {
                "tag": spec["tag"],
                "case": spec["case"],
                "water": spec["water"],
                "stage": "timeout",
                "recovery": None,
                "ok": False,
                "time": round(time.time() - t0, 2),
                "iters": None,
                "lcow": None,
                "error": f"timeout {timeout}s",
            }
        )
    elif not rows:
        rows.append(
            {
                "tag": spec["tag"],
                "case": spec["case"],
                "water": spec["water"],
                "stage": "crash",
                "recovery": None,
                "ok": False,
                "time": round(time.time() - t0, 2),
                "iters": None,
                "lcow": None,
                "error": (proc.stderr or "")[-300:],
            }
        )
    status = sum(r["ok"] for r in rows if r["stage"] == "solve")
    print(
        f"{spec['tag']:>8} {spec['case']:>18} {spec['water']:>9}: {status} solves ok, "
        f"{time.time() - t0:.0f}s",
        flush=True,
    )
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-case")
    p.add_argument("--out", default="bench_softening_ro.csv")
    p.add_argument("--cases", default=",".join(CASES))
    p.add_argument("--waters", default="BGW,SW_HPRO,BGW_1500")
    p.add_argument("--recoveries", default=",".join(map(str, RECOVERIES)))
    p.add_argument("--linear-solver", default="mumps")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--timeout", type=int, default=3600)
    p.add_argument("--retw", default=DEFAULT_RETW)
    p.add_argument(
        "--rktpse-src",
        default="",
        help="reaktoro_pse src dir to use instead of installed one (e.g. a baseline worktree)",
    )
    p.add_argument("--tag", default="new", help="label for this code version")
    args = p.parse_args()
    if args.run_case:
        run_case(json.loads(args.run_case))
        return
    specs = [
        {
            "case": c,
            "water": w,
            "overrides": CASES[c][2],
            "recoveries": [float(r) for r in args.recoveries.split(",")],
            "linear_solver": args.linear_solver,
            "retw": args.retw,
            "rktpse_src": args.rktpse_src,
            "tag": args.tag,
        }
        for w in args.waters.split(",")
        for c in args.cases.split(",")
    ]
    rows = []
    with ThreadPoolExecutor(args.workers) as ex:
        for result in ex.map(lambda s: launch(s, args.timeout), specs):
            rows.extend(result)
    fields = ["tag", "case", "water", "stage", "recovery", "ok", "time", "iters", "lcow", "error"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
