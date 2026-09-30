"""Capstone figure: mutation-mover retrieval AUC across the model gates.
Values are the recorded outputs of gate2 / loop_gate / powered_loop_gate /
gate3_3d (see RESULTS.md). CIs are the ORIGINAL per-pair bootstrap; the
scripts now use a residue-cluster bootstrap -- update these after a re-run."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# label, AUC, ci_low, ci_high (None if not yet bootstrapped), n_movers, colour
LOCAL, THREED = "#4C78A8", "#F58518"
tests = [
    ("All T4L pairs\n(local, Gate 2)",  0.524, None, None, 72, LOCAL),
    ("Loops, T4L\n(local, n=10)",       0.700, 0.49, 0.89, 10, LOCAL),
    ("Loops, powered\n(local, n=33)",   0.591, 0.49, 0.69, 33, LOCAL),
    ("All T4L, post-QC\n(local+3D, Gate 5)", 0.590, 0.52, 0.66, 70, THREED),
]
fig, ax = plt.subplots(figsize=(7.6, 4.8))
for i, (lab, auc, lo, hi, n, c) in enumerate(tests):
    if lo is not None:
        ax.errorbar(i, auc, yerr=[[auc - lo], [hi - auc]], fmt="o",
                    color=c, capsize=6, ms=9, lw=2)
    else:
        ax.plot(i, auc, "o", color=c, ms=9, mfc="white", mew=2)
        ax.annotate("no CI", (i, auc), textcoords="offset points",
                    xytext=(0, 12), ha="center", fontsize=8, color="#666666")
    ax.annotate(f"{auc:.2f}", (i, auc), textcoords="offset points",
                xytext=(12, -2), fontsize=10)
ax.axhline(0.5, color="#E45756", ls="--", lw=1.2, label="chance (0.50)")
ax.plot([], [], "o", color=LOCAL, label="local sequence only")
ax.plot([], [], "o", color=THREED, label="+ crude 3D contact composition")
ax.set_xticks(range(len(tests)))
ax.set_xticklabels([t[0] for t in tests], fontsize=9)
ax.set_ylim(0.40, 0.95)
ax.set_ylabel("mutation-mover retrieval AUC (90% CI)")
ax.set_title("Local sequence never robustly clears chance\n"
             "(every bootstrapped local-only CI touches 0.5; the 3D lift itself is untested)")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
fig.savefig("summary_auc.png", dpi=130)
print("wrote summary_auc.png")
