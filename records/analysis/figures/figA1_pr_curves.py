# =============================================================================
# Figure 4 - precision-recall curves for HattCI, threshold-independent
# Self-contained: matplotlib + numpy only (preinstalled in Colab). Paste & run.
#
# Each curve = HattCI calls on one genome, ranked by Viterbi score and swept
# from the highest-scoring call downward. A point at recall r, precision p means:
# taking every call at or above some Viterbi cut-off recovers fraction r of the
# IntegronFinder --local-max reference attC sites, and fraction p of those calls
# are real (matched within 100 bp). Average Precision (AP) = area under the
# curve; it summarises accuracy across every possible score cut-off, so it does
# not depend on the Viterbi >= 2 threshold used for the F1 table.
#
# Each species has its own marker shape as well as its own colour, so the curves
# separate in colour, in greyscale, and for colour-blind readers.
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

# ---- PR points (from score_ap.py; already downsampled) ----------------------
CURVES = {
 "vc569B": [(0.0085,1.0000),(0.0593,0.8750),(0.1102,0.8667),(0.1610,0.8261),
            (0.2119,0.8621),(0.2627,0.8378),(0.3136,0.8222),(0.3644,0.8269),
            (0.4153,0.8305),(0.4661,0.7857),(0.5169,0.7262),(0.5678,0.7204),
            (0.6186,0.7019),(0.6695,0.6752),(0.7203,0.6296),(0.7712,0.5833),
            (0.8220,0.5389),(0.8729,0.5228),(0.9237,0.4802),(0.9492,0.3446)],
 "vp178": [(0.0070,1.0000),(0.0490,0.8750),(0.0909,0.9286),(0.1329,0.9500),
           (0.1748,0.8929),(0.2168,0.9118),(0.2587,0.9024),(0.3007,0.8958),
           (0.3427,0.8167),(0.3846,0.7051),(0.4266,0.6354),(0.4685,0.5583),
           (0.5105,0.5214),(0.5524,0.4675),(0.5944,0.4208),(0.6364,0.3699),
           (0.6783,0.3489),(0.7203,0.2543),(0.7622,0.2266),(0.7762,0.2143)],
 "vv275": [(0.0046,1.0000),(0.0413,0.6429),(0.0780,0.7727),(0.1147,0.8333),
           (0.1514,0.8684),(0.1881,0.7736),(0.2248,0.7424),(0.2615,0.7215),
           (0.2982,0.6771),(0.3349,0.6697),(0.3716,0.6328),(0.4083,0.5894),
           (0.4450,0.5419),(0.4817,0.5024),(0.5183,0.4689),(0.5550,0.4276),
           (0.5917,0.3874),(0.6284,0.3723),(0.6651,0.3396)],
 "xcc":   [(0.0455,1.0),(0.0909,1.0),(0.1364,1.0),(0.1818,1.0),(0.2273,1.0),
           (0.2727,1.0),(0.3182,1.0),(0.3636,1.0),(0.4091,1.0),(0.4545,1.0),
           (0.5000,1.0),(0.5455,1.0),(0.5909,1.0),(0.6364,1.0),(0.6818,1.0),
           (0.7273,1.0),(0.7727,1.0),(0.8182,1.0),(0.8636,1.0),(0.9091,1.0),
           (0.9545,1.0)],
 "afDSM": [(0.0079,0.0222),(0.0159,0.0274),(0.0238,0.0227),(0.0317,0.0301),
           (0.0397,0.0298),(0.0476,0.0286),(0.0556,0.0312),(0.0635,0.0349),
           (0.0714,0.0383),(0.0794,0.0377),(0.0873,0.0414),(0.0952,0.0438),
           (0.1032,0.0469),(0.1111,0.0504),(0.1190,0.0466)],
}

# species label, strain, colour, marker, Average Precision (from score_ap.py)
META = {
 "xcc":    (r"$\mathit{X.\ campestris}$",      "ATCC 33913", "#2e8b57", "s", 0.955),
 "vc569B": (r"$\mathit{V.\ cholerae}$",        "569B",       "#1f6f8b", "o", 0.692),
 "vp178":  (r"$\mathit{V.\ parahaemolyticus}$","ATCC 17802", "#e08a1e", "^", 0.516),
 "vv275":  (r"$\mathit{V.\ vulnificus}$",      "ATCC 27562", "#8e6fb3", "D", 0.414),
 "afDSM":  (r"$\mathit{A.\ fischeri}$",        "DSM 507",    "#c0392b", "X", 0.004),
}
ORDER = ["xcc", "vc569B", "vp178", "vv275", "afDSM"]   # legend top->bottom by AP
INK, MUTE = "#26303a", "#6b7680"

fig, ax = plt.subplots(figsize=(6.6, 6.0))

for g in ORDER:
    sp, st, c, mk, ap = META[g]
    xs = [r for r, _ in CURVES[g]]
    ys = [p for _, p in CURVES[g]]
    lab = "%s %s  (AP = %s)" % (sp, st, ("0.004" if ap < 0.01 else "%.2f" % ap))
    ms = 6.5 if mk in ("D",) else (7.5 if mk in ("s", "^", "X") else 6.0)
    ax.plot(xs, ys, color=c, lw=1.8, zorder=4,
            marker=mk, ms=ms, mfc=c, mec="white", mew=0.6, label=lab)

# axes
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.03, 1.05)
ax.set_xticks(np.arange(0, 1.01, 0.2)); ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.set_xlabel("Recall  (reference attC sites recovered)", fontsize=12, color=INK, labelpad=6)
ax.set_ylabel("Precision  (calls matching the reference)", fontsize=12, color=INK, labelpad=6)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.spines["left"].set_color("#c3ccd2"); ax.spines["bottom"].set_color("#c3ccd2")
ax.tick_params(colors=MUTE, labelsize=9.5, length=4, width=0.9, color="#c3ccd2")

leg = ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.0), fontsize=8.6,
                frameon=True, borderpad=0.8, edgecolor="#c3ccd2", facecolor="white",
                framealpha=1.0, handlelength=1.8, handletextpad=0.7, labelspacing=0.7)
leg.get_frame().set_linewidth(0.8)
for t in leg.get_texts(): t.set_color(INK)

fig.tight_layout()
fig.savefig("fig4_pr_curves.png", bbox_inches="tight")
fig.savefig("fig4_pr_curves.pdf", bbox_inches="tight")
print("wrote fig4_pr_curves.png / .pdf")
try:
    from IPython.display import Image, display
    display(Image(filename="fig4_pr_curves.png"))
except Exception:
    pass
