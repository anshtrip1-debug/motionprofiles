#!/usr/bin/env python3
"""
plot_profile - draw the CSV that profile.cpp writes.

    ./profile --csv --distance 400 > move.csv
    python3 plot_profile.py move.csv -o move.png

Four stacked panels sharing a time axis: position, velocity, acceleration,
jerk. The trapezoid and the jerk-limited profile are drawn together so the
difference is visible rather than described. Limit lines are dashed, which
makes it obvious when a profile is riding a limit and when it never gets
near one.
"""

import argparse
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    cols = {k: [float(r[k]) for r in rows] for k in rows[0]}
    return cols


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("csvfile")
    ap.add_argument("-o", "--out", default="move.png")
    ap.add_argument("--vmax", type=float, default=None)
    ap.add_argument("--amax", type=float, default=None)
    ap.add_argument("--jmax", type=float, default=None)
    args = ap.parse_args(argv)

    d = read_csv(args.csvfile)
    t = d["t"]

    trap = "#B03A2E"   # brick
    ds = "#1E3AA0"     # anodized blue
    lim = "#9AA39A"

    fig, ax = plt.subplots(4, 1, figsize=(8.2, 9.4), sharex=True)
    fig.subplots_adjust(hspace=0.16, left=0.12, right=0.97, top=0.94,
                        bottom=0.07)

    panels = [
        ("position", "trap_q", "ds_q", None, None),
        ("velocity", "trap_v", "ds_v", args.vmax, "vmax"),
        ("acceleration", "trap_a", "ds_a", args.amax, "amax"),
        ("jerk", None, "ds_j", args.jmax, "jmax"),
    ]

    for a, (label, ka, kb, limit, limname) in zip(ax, panels):
        if ka:
            a.plot(t, d[ka], color=trap, lw=1.6, label="trapezoidal")
        a.plot(t, d[kb], color=ds, lw=1.6, label="jerk limited")
        if limit:
            for sign in (1, -1):
                a.axhline(sign * limit, color=lim, lw=0.9, ls="--")
            a.text(t[-1], limit, " " + limname, va="center", fontsize=8,
                   color=lim)
        a.set_ylabel(label)
        a.grid(True, color="#E2E5E1", lw=0.7)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)

    ax[3].set_xlabel("time (s)")
    ax[3].text(0.015, 0.5,
               "the trapezoid has no jerk trace: its acceleration steps,\n"
               "so jerk is unbounded at four instants",
               transform=ax[3].transAxes, fontsize=8, color="#5A625A",
               va="center",
               bbox=dict(fc="white", ec="none", alpha=0.85, pad=3))
    ax[0].legend(loc="lower right", frameon=False, fontsize=9)
    ax[0].set_title("Rest to rest move, same limits, two profiles",
                    fontsize=11, loc="left")

    fig.savefig(args.out, dpi=150)
    print("wrote %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
