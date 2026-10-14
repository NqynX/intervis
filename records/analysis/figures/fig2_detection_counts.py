# =============================================================================
# Figure 1 - attC detection per tool (publication-quality, self-contained)
# matplotlib + numpy only (preinstalled in Colab). Paste & run.
# =============================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
import numpy as np

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.linewidth": 0.9, "savefig.dpi": 300,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

# ---- data (attC sites detected) ---------------------------------------------
species = [r"$\mathit{V.\ cholerae}$", r"$\mathit{A.\ fischeri}$",
           r"$\mathit{V.\ vulnificus}$", r"$\mathit{V.\ parahaemolyticus}$",
           r"$\mathit{X.\ campestris}$"]
strain  = ["569B", "DSM 507", "ATCC 27562", "ATCC 17802", "ATCC 33913"]
IFL = [118, 126, 218, 143, 22]        # IntegronFinder --local-max
IFD = [45,  22,  158, 92,  22]        # IntegronFinder default
HAT = [133, 55,  153, 171, 50]        # HattCI (Viterbi >= 2)

# cohesive palette: two IntegronFinder settings share a hue family; HattCI = accent
C_LOCAL = "#17465c"   # IntegronFinder, local search (deep teal-blue)
C_DEF   = "#8fb7c7"   # IntegronFinder, default (same family, lighter)
C_HAT   = "#d98c3f"   # HattCI (warm accent)
INK, MUTE, FAINT = "#26303a", "#6b7680", "#e9edf0"

fig, ax = plt.subplots(figsize=(7.8, 4.5))
x = np.arange(len(species)); w = 0.25

ax.bar(x - w, IFL, w, color=C_LOCAL, zorder=3, label="IntegronFinder, local search")
ax.bar(x,     IFD, w, color=C_DEF,   zorder=3, label="IntegronFinder, default")
ax.bar(x + w, HAT, w, color=C_HAT,   zorder=3, label="HattCI")

# --- reported full-array cassette counts (only where a complete count exists) ---
C_REF = "#3a3a3a"
reported = {0: (150, u"≈150"), 2: (230, "230")}   # index -> (value, label)
for i, (val, lab) in reported.items():
    ax.plot([i - 1.7*w, i + 1.7*w], [val, val], color=C_REF, lw=1.6,
            ls=(0, (5, 2)), zorder=6)
    ax.text(i + 1.8*w, val, " " + lab, ha="left", va="center", fontsize=8.6, color=C_REF)

# --- axes ---
ax.set_ylim(0, 250); ax.set_yticks(range(0, 251, 50))
ax.set_ylabel("attC sites detected", fontsize=12, color=INK, labelpad=8)
ax.set_xticks(x)
ax.set_xticklabels([s + "\n" + t for s, t in zip(species, strain)], fontsize=9.5, color=INK)
ax.set_xlim(-0.6, len(species)-0.4)

for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.spines["left"].set_color("#c3ccd2"); ax.spines["bottom"].set_color("#c3ccd2")
ax.tick_params(axis="y", colors=MUTE, labelsize=9.5, length=4, width=0.9, color="#c3ccd2")
ax.tick_params(axis="x", length=0, pad=6)

handles, labels = ax.get_legend_handles_labels()
handles.append(Line2D([0], [0], color=C_REF, lw=1.6, ls=(0, (5, 2))))
labels.append("reported cassettes")
leg = ax.legend(handles, labels, fontsize=9.0, loc="upper left",
                bbox_to_anchor=(0.015, 0.985),
                handlelength=1.4, handletextpad=0.6, labelspacing=0.55,
                frameon=True, borderpad=0.8, edgecolor="#c3ccd2", facecolor="white",
                framealpha=1.0)
leg.get_frame().set_linewidth(0.8)
for t in leg.get_texts(): t.set_color(INK)

fig.tight_layout()
fig.savefig("fig1_detection.png", bbox_inches="tight")
fig.savefig("fig1_detection.pdf", bbox_inches="tight")
print("wrote fig1_detection.png / .pdf")
try:
    from IPython.display import Image, display
    display(Image(filename="fig1_detection.png"))
except Exception:
    pass
