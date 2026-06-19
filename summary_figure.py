"""Capstone figure: mutation-mover retrieval AUC across the model gates.
Values are the established outputs of gate2 / loop_gate / powered_loop_gate."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# label, AUC, ci_low, ci_high (None if not bootstrapped), n_movers
tests = [
    ("Core\n(global)",        0.524, None, None, 72),
    ("Loops\n(T4L, n=10)",    0.700, 0.49, 0.89, 10),
    ("Loops\n(powered, n=33)", 0.591, 0.49, 0.69, 33),
]
fig, ax = plt.subplots(figsize=(6.4, 4.6))
for i, (lab, auc, lo, hi, n) in enumerate(tests):
    if lo is not None:
        ax.errorbar(i, auc, yerr=[[auc - lo], [hi - auc]], fmt="o",
                    color="#4C78A8", capsize=6, ms=9, lw=2)
    else:
        ax.plot(i, auc, "o", color="#888888", ms=9)
    ax.annotate(f"{auc:.2f}", (i, auc), textcoords="offset points",
                xytext=(12, -2), fontsize=10)
ax.axhline(0.5, color="#E45756", ls="--", lw=1.2, label="chance (0.50)")
ax.set_xticks(range(len(tests)))
ax.set_xticklabels([t[0] for t in tests])
ax.set_ylim(0.40, 0.95)
ax.set_ylabel("mutation-mover retrieval AUC")
ax.set_title("Local sequence never robustly clears chance\n"
             "(every CI touches 0.5; the only point above it was n=10)")
ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig("summary_auc.png", dpi=130)
print("wrote summary_auc.png")
