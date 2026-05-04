import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import textwrap

BASE = Path(__file__).resolve().parent
CSV_PATH = BASE / "alphafold_timeline_data_full.csv"
OUT_DIR = BASE

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


def stagger_for_position(pos):
    """Match main timeline: alternating above / below the line."""
    y_offset = 0.15 if pos % 2 == 0 else -0.25
    rotation = 45 if pos % 2 == 0 else -45
    ha = "right" if pos % 2 == 0 else "left"
    above_line = pos % 2 == 0
    return y_offset, rotation, ha, above_line


def event_slug(title):
    s = re.sub(r"[^\w\-]+", "_", title.strip(), flags=re.UNICODE)
    return re.sub(r"_+", "_", s).strip("_") or "event"


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


def draw_event_labels(ax, highlight_pos=None):
    for pos, row in df.iterrows():
        y_off, rotation, ha, _ = stagger_for_position(pos)
        y_lab = 1 + y_off
        fs = LABEL_FONTSIZE + 2 if pos == highlight_pos else LABEL_FONTSIZE
        edge = "black" if pos == highlight_pos else None
        lw = 1.2 if pos == highlight_pos else 0
        ax.text(
            row["Date"],
            y_lab,
            row["Event"],
            rotation=rotation,
            ha=ha,
            va="center",
            fontsize=fs,
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor="lightblue",
                edgecolor=edge,
                linewidth=lw,
            ),
        )


def add_description_for_pos(ax, pos, row):
    y_off, _, _, above_line = stagger_for_position(pos)
    y_lab = 1 + y_off
    desc = str(row.get("Description", "") or "")
    if not desc.strip():
        return
    wrapped = textwrap.fill(desc, width=52)
    gap = 0.12
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


def timeline_ylim():
    """Room for rotated labels plus optional description bands."""
    return (0.2, 1.55)


# Full overview (event names only)
fig, ax = plt.subplots(figsize=(14, 6))
plot_timeline_base(ax)
draw_event_labels(ax, highlight_pos=None)
ax.set_ylim(*timeline_ylim())
plt.tight_layout()
overview_path = OUT_DIR / "alphafold_timeline.png"
plt.savefig(overview_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Plot saved as {overview_path}")

# One PNG per event: same timeline context, description near that label only
for pos, row in df.iterrows():
    fig, ax = plt.subplots(figsize=(14, 7))
    plot_timeline_base(ax)
    draw_event_labels(ax, highlight_pos=pos)
    add_description_for_pos(ax, pos, row)
    ax.set_ylim(*timeline_ylim())
    plt.tight_layout()
    per_path = OUT_DIR / f"alphafold_timeline_{event_slug(row['Event'])}.png"
    plt.savefig(per_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {per_path}")
