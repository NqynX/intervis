# =============================================================================
# Figure 5 - recall vs attC divergence on synthetic genomes (known truth)
# Self-contained: matplotlib + numpy only (preinstalled in Colab). Paste & run.
#
# Two panels, same axes. Left: attC repeats diverged by random substitution
# (secondary structure broken). Right: diverged by the same amount with the
# imperfect inverted repeat preserved (compensatory substitutions). Each line is
# one detector; y = fraction of the known inserted attC sites recovered.
# IntegronFinder holds up under structure-preserving divergence but collapses
# under naive divergence; HattCI declines the same way in both panels.
# =============================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.linewidth": 0.9, "savefig.dpi": 300,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

X = [0, 5, 10, 15, 20, 25, 30]

naive = {
 "IFloc":  [1.000, 0.688, 0.338, 0.062, 0.087, 0.025, 0.025],
 "IFdef":  [0.400, 0.300, 0.113, 0.000, 0.025, 0.000, 0.000],
 "HattCI": [0.738, 0.500, 0.338, 0.200, 0.163, 0.075, 0.037],
 "BacAnt": [0.400, 0.300, 0.113, 0.000, 0.025, 0.000, 0.000],
}
structure = {
 "IFloc":  [1.000, 0.787, 0.613, 0.475, 0.400, 0.362, 0.200],
 "IFdef":  [0.400, 0.287, 0.200, 0.150, 0.087, 0.100, 0.000],
 "HattCI": [0.738, 0.487, 0.312, 0.250, 0.138, 0.163, 0.037],
 "BacAnt": [0.400, 0.287, 0.200, 0.150, 0.087, 0.100, 0.000],
}

# label, colour, marker  (BacAnt drawn as open markers: its recall == IF default)
STYLE = [
 ("IFloc",  "IntegronFinder, local search", "#17465c", "o", True),
 ("IFdef",  "IntegronFinder, default",      "#8fb7c7", "s", True),
 ("HattCI", "HattCI",                        "#d98c3f", "^", True),
 ("BacAnt", "BacAnt (recall = IF default)",  "#9aa0a6", "D", False),
]
INK, MUTE = "#26303a", "#6b7680"

fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.0), sharey=True)
panels = [("Random substitutions", naive),
          ("Compensatory substitutions", structure)]

for ax, (title, data) in zip(axes, panels):
    for key, lab, c, mk, filled in STYLE:
        y = data[key]
        if filled:
            ax.plot(X, y, color=c, lw=1.9, marker=mk, ms=6.5, mfc=c, mec="white",
                    mew=0.6, zorder=4, label=lab)
        else:
            ax.plot(X, y, color=c, lw=0, marker=mk, ms=7.0, mfc="none", mec=c,
                    mew=1.2, zorder=5, label=lab)
    ax.set_xlim(-1.5, 31.5); ax.set_ylim(-0.03, 1.05)
    ax.set_xticks(X); ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_xlabel("attC sequence divergence (%)", fontsize=11.5, color=INK, labelpad=6)
    ax.set_title(title, fontsize=11, color=INK, pad=8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color("#c3ccd2"); ax.spines["bottom"].set_color("#c3ccd2")
    ax.tick_params(colors=MUTE, labelsize=9.5, length=4, width=0.9, color="#c3ccd2")

axes[0].set_ylabel("Recall  (known attC sites recovered)", fontsize=11.5, color=INK, labelpad=8)

leg = axes[1].legend(loc="upper right", bbox_to_anchor=(0.99, 0.985), fontsize=9,
                     frameon=True, borderpad=0.8, edgecolor="#c3ccd2",
                     facecolor="white", framealpha=1.0, handlelength=1.6,
                     handletextpad=0.6, labelspacing=0.55)
leg.get_frame().set_linewidth(0.8)
for t in leg.get_texts(): t.set_color(INK)

fig.tight_layout()
fig.savefig("fig5_synthetic_recall.png", bbox_inches="tight")
fig.savefig("fig5_synthetic_recall.pdf", bbox_inches="tight")
print("wrote fig5_synthetic_recall.png / .pdf")
try:
    from IPython.display import Image, display
    display(Image(filename="fig5_synthetic_recall.png"))
except Exception:
    pass
