"""Summarize bench_softening_ro.py csv outputs into a markdown table.

Usage: python dev/summarize_bench.py results1.csv [results2.csv ...] [--ref new:lbfgs_gn]
"""

import argparse
import csv
from collections import defaultdict
from statistics import median


def load(paths):
    rows = []
    for p in paths:
        with open(p) as f:
            for r in csv.DictReader(f):
                r["ok"] = r["ok"] == "True"
                r["iters"] = int(r["iters"]) if r["iters"] else None
                r["time"] = float(r["time"]) if r["time"] else 0.0
                rows.append(r)
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("paths", nargs="+")
    p.add_argument("--ref", default="baseline:lbfgs_gn", help="tag:case used as reference")
    args = p.parse_args()
    rows = load(args.paths)

    keys = sorted({(r["tag"], r["case"]) for r in rows})
    waters = sorted({r["water"] for r in rows})
    solves = defaultdict(dict)  # key -> {(water, recovery): row}
    inits = defaultdict(dict)
    for r in rows:
        k = (r["tag"], r["case"])
        if r["stage"] == "solve":
            solves[k][(r["water"], r["recovery"])] = r
        elif r["stage"] in ("init", "timeout", "crash"):
            inits[k][r["water"]] = r
    # points attempted: union over keys of solve points (+ points never reached
    # because init failed count as failures)
    all_points = sorted({pt for k in keys for pt in solves[k]})
    def solved(k, pt):
        return pt in solves[k] and solves[k][pt]["ok"]

    # iteration comparison uses points solved by every case that solved everything
    # (cases that miss one of those points show "-")
    robust = [k for k in keys if all(solved(k, pt) for pt in all_points)] or keys
    common = [pt for pt in all_points if all(solved(k, pt) for k in robust)]
    ref = tuple(args.ref.split(":"))

    def stats(k):
        ok = sum(1 for pt in all_points if pt in solves[k] and solves[k][pt]["ok"])
        init_ok = sum(1 for w in waters if inits[k].get(w, {}).get("ok"))
        has_common = all(solved(k, pt) for pt in common)
        it_common = [solves[k][pt]["iters"] for pt in common] if has_common else None
        init_iters = sum(r["iters"] or 0 for r in inits[k].values() if r["ok"])
        all_iters = [r["iters"] for r in solves[k].values() if r["ok"]]
        t = sum(r["time"] for r in solves[k].values()) + sum(
            r["time"] for r in inits[k].values()
        )
        failed_waters = defaultdict(int)
        for w in waters:
            if not inits[k].get(w, {}).get("ok"):
                failed_waters[w] += sum(1 for pt in all_points if pt[0] == w)
                continue
            for pt in all_points:
                if pt[0] == w and not (pt in solves[k] and solves[k][pt]["ok"]):
                    failed_waters[w] += 1
        return dict(
            ok=ok,
            init_ok=init_ok,
            sum_common=sum(it_common) if has_common else None,
            init_iters=init_iters,
            med_all=median(all_iters) if all_iters else float("nan"),
            max_all=max(all_iters) if all_iters else float("nan"),
            time=t,
            failed=dict(failed_waters),
        )

    s = {k: stats(k) for k in keys}
    ref_sum = s[ref]["sum_common"] if ref in s else None
    print(
        f"{len(all_points)} solve points over {len(waters)} waters, "
        f"{len(common)} solved by every case (iteration comparison uses those)\n"
    )
    print(
        "| tag | case | init ok | init iters | solves ok | iters on common pts | vs ref | median iters | max iters | total time (s) | failures by water |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    order = sorted(
        keys, key=lambda k: (-s[k]["ok"], s[k]["sum_common"] if s[k]["sum_common"] is not None else 1e18)
    )
    for k in order:
        st = s[k]
        rel = f"{st['sum_common'] / ref_sum:.2f}" if ref_sum and st["sum_common"] else "-"
        sc = st["sum_common"] if st["sum_common"] is not None else "-"
        fails = ", ".join(f"{w}:{n}" for w, n in sorted(st["failed"].items())) or "-"
        print(
            f"| {k[0]} | {k[1]} | {st['init_ok']}/{len(waters)} | {st['init_iters']} | "
            f"{st['ok']}/{len(all_points)} | {sc} | {rel} | {st['med_all']:.0f} | {st['max_all']:.0f} | "
            f"{st['time']:.0f} | {fails} |"
        )


if __name__ == "__main__":
    main()
