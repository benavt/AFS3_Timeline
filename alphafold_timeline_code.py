import re
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import textwrap

BASE = Path(__file__).resolve().parent
CSV_PATH = BASE / "alphafold_timeline_data_full_with_labels.csv"
# Main overview (no descriptions) stays next to the script; detail PNGs go here
OUT_DIR = BASE
OUTPUTS_DIR = BASE / "outputs"


def normalize_date_cell(raw):
    s = str(raw).strip()
    if " / " in s:
        s = s.split(" / ")[0].strip()
    m = re.match(r"^(\d{4})\s*[\u2013\-]\s*(\d{4})$", s)
    if m:
        s = f"Jan {m.group(1)}"
    return s


# Load data
df = pd.read_csv(CSV_PATH)
df["Date"] = pd.to_datetime(df["Date"].map(normalize_date_cell), format="%b %Y")
df = df.sort_values("Date").reset_index(drop=True)

_x_num = np.asarray(mdates.date2num(df["Date"]), dtype=float)
_span = float(_x_num.max() - _x_num.min())
if _span <= 0:
    _span = 1e-6

# Event label column -> face color (case-insensitive; CASP rounds by Label or Event name)
FACE_CONFORMATIONAL = "lightblue"
FACE_PROTEIN_DESIGN = "#E5B4B9"  # pastel maroon / dusty rose
FACE_CASP = "#B8E0C8"  # pastel green


def category_facecolor(row):
    """Use CSV ``Label``; CASP rounds are green if Label is CASP or Event is ``CASP##``."""
    evt = str(row.get("Event", "") or "").strip()
    lbl = str(row.get("Label", "") or "").strip().lower()
    if re.match(r"^CASP\d+\s*$", evt, flags=re.I):
        lbl = "casp"
    if lbl == "protein design":
        return FACE_PROTEIN_DESIGN
    if lbl == "casp":
        return FACE_CASP
    # Conformational Modeling, unknown
    return FACE_CONFORMATIONAL


def event_slug(title):
    s = re.sub(r"[^\w\-]+", "_", title.strip(), flags=re.UNICODE)
    return re.sub(r"_+", "_", s).strip("_") or "event"


def estimate_label_halfwidth_data(event_text, fontsize, fig_w_in, subplot_w_frac=0.82):
    """Rough half-width on the date axis using character count and subplot width."""
    chars = max(len(event_text), 4)
    # ~0.6 * fontsize/72 inch per character, axis spans subplot_w_frac of figure
    inch_w = (chars * 0.65 * fontsize / 72.0) * 1.15
    data_w = inch_w / (fig_w_in * subplot_w_frac) * _span
    return 0.5 * data_w


def assign_label_y(
    fontsize,
    fig_w_in=14,
    subplot_w_frac=0.82,
    line_y=1.0,
    row_step=0.11,
    first_gap=0.08,
):
    """
    Alternate above/below the timeline with horizontal non-overlap in date space.
    Returns list of (y_center, above_line) for each row in df order.
    """
    n = len(df)
    halfw = [
        estimate_label_halfwidth_data(str(df.loc[i, "Event"]), fontsize, fig_w_in, subplot_w_frac)
        for i in range(n)
    ]
    x = _x_num
    above_row = [0] * n
    below_row = [0] * n

    for i in range(n):
        want_above = i % 2 == 0
        side_rows = above_row if want_above else below_row
        for lane in range(30):
            ok = True
            low_i = x[i] - halfw[i]
            high_i = x[i] + halfw[i]
            for j in range(i):
                same_side = (j % 2 == 0) == want_above
                if not same_side:
                    continue
                if side_rows[j] != lane:
                    continue
                low_j = x[j] - halfw[j]
                high_j = x[j] + halfw[j]
                if not (high_i < low_j or low_i > high_j):
                    ok = False
                    break
            if ok:
                side_rows[i] = lane
                break

    y_out = []
    for i in range(n):
        want_above = i % 2 == 0
        lane = above_row[i] if want_above else below_row[i]
        dist = first_gap + lane * row_step
        y = line_y + dist if want_above else line_y - dist
        y_out.append((y, want_above))
    return y_out


def plot_timeline_base(ax):
    ax.plot(df["Date"], [1] * len(df), "-o", color="blue", linewidth=3, markersize=8)
    ax.set_title("AlphaFold & AFsample Timeline", fontsize=16)
    ax.set_xlabel("Year", fontsize=12)
    ax.set_yticks([])
    ax.grid(True, axis="x", linestyle="--", alpha=0.7)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))


LABEL_FONTSIZE = 14
DESC_FONTSIZE = 11
TITLE_FONTSIZE = 16
AXIS_FONTSIZE = 12

ARROW_PROPS = dict(
    arrowstyle="->",
    lw=1.0,
    color="0.45",
    shrinkA=2,
    shrinkB=2,
    mutation_scale=12,
    connectionstyle="arc3,rad=0",
)


def draw_event_labels(ax, highlight_pos=None, label_positions=None):
    if label_positions is None:
        label_positions = assign_label_y(LABEL_FONTSIZE)
    for pos, row in df.iterrows():
        y_lab, _ = label_positions[pos]
        fs = LABEL_FONTSIZE + 2 if pos == highlight_pos else LABEL_FONTSIZE
        face = category_facecolor(row)
        edge = "black" if pos == highlight_pos else "0.35"
        lw = 1.2 if pos == highlight_pos else 0.6
        ax.annotate(
            row["Event"],
            xy=(row["Date"], 1.0),
            xytext=(row["Date"], y_lab),
            textcoords="data",
            ha="center",
            va="center",
            rotation=0,
            fontsize=fs,
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor=face,
                edgecolor=edge,
                linewidth=lw,
            ),
            arrowprops=ARROW_PROPS,
        )


def add_category_legend(ax):
    handles = [
        Patch(
            facecolor=FACE_CONFORMATIONAL,
            edgecolor="0.35",
            linewidth=0.6,
            label="Conformational modeling",
        ),
        Patch(
            facecolor=FACE_CASP,
            edgecolor="0.35",
            linewidth=0.6,
            label="CASP",
        ),
        Patch(
            facecolor=FACE_PROTEIN_DESIGN,
            edgecolor="0.35",
            linewidth=0.6,
            label="Protein design",
        ),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=11, framealpha=0.95)


def add_description_for_pos(ax, pos, row, label_positions):
    y_lab, above_line = label_positions[pos]
    desc = str(row.get("Description", "") or "")
    if not desc.strip():
        return
    wrapped = textwrap.fill(desc, width=52)
    gap = 0.1
    if above_line:
        y_desc = y_lab + gap
        va = "bottom"
    else:
        y_desc = y_lab - gap
        va = "top"
    ax.text(
        row["Date"],
        y_desc,
        wrapped,
        ha="center",
        va=va,
        fontsize=DESC_FONTSIZE,
        linespacing=1.15,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.92, edgecolor="gray"),
    )


def timeline_ylim(label_positions):
    ys = [y for y, _ in label_positions]
    top = max(ys) if ys else 1.2
    bot = min(ys) if ys else 0.8
    pad = 0.12
    return (min(bot - pad, 0.75), max(top + pad, 1.35))


# Precompute layout once so overview and per-event figures match
LABEL_POS = assign_label_y(LABEL_FONTSIZE)

OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

# Full overview (event names only)
fig, ax = plt.subplots(figsize=(14, 8))
plot_timeline_base(ax)
draw_event_labels(ax, highlight_pos=None, label_positions=LABEL_POS)
add_category_legend(ax)
ax.set_ylim(*timeline_ylim(LABEL_POS))
plt.tight_layout()
overview_path = OUT_DIR / "alphafold_timeline.png"
plt.savefig(overview_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Plot saved as {overview_path}")

# One PNG per event: same timeline context, description near that label only
for pos, row in df.iterrows():
    fig, ax = plt.subplots(figsize=(14, 9))
    plot_timeline_base(ax)
    draw_event_labels(ax, highlight_pos=pos, label_positions=LABEL_POS)
    add_category_legend(ax)
    add_description_for_pos(ax, pos, row, LABEL_POS)
    ax.set_ylim(*timeline_ylim(LABEL_POS))
    # Extra headroom / footroom for wrapped description
    y0, y1 = ax.get_ylim()
    _, above = LABEL_POS[pos]
    if above:
        ax.set_ylim(y0, y1 + 0.18)
    else:
        ax.set_ylim(y0 - 0.18, y1)
    plt.tight_layout()
    per_path = OUTPUTS_DIR / f"alphafold_timeline_{event_slug(row['Event'])}.png"
    plt.savefig(per_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {per_path}")
