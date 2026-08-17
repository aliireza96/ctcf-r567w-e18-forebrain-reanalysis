#!/usr/bin/env python3
"""
Shared panel writer: every panel is emitted twice, captioned and bare.

WHY. Figure legends are written in the manuscript, not printed on the panel. But the grey
explanatory footnotes are genuinely useful while reviewing a panel on its own, so they are
kept in one copy and stripped from the other. Both come from the SAME figure object, so a
rebuild cannot leave the two out of step, which is the failure mode of maintaining a separate
"clean" version by hand.

  <stem>.pdf        for review: title, subtitle, footnote, colour keys
  <stem>_bare.pdf   for assembly: axes, tick labels, track names, colour keys only

WHAT COUNTS AS A CAPTION. Anything that explains rather than labels. Titles and footnotes go;
axis labels, tick labels, in-panel value annotations and colour keys stay, because without
them the panel cannot be read at all.

The caption text is not lost when it is stripped: each panel script keeps it in its docstring
and its source CSV header, which is where the Results text should be drawn from.

Date: 2026-08-16
"""
import os
import matplotlib.pyplot as plt


def save_twin(fig, out_dir, stem, captions=(), bare_size=None, bare_adjust=None,
              dpi=210, close=True, verbose=True):
    """Write <stem>.pdf/.png with captions, then <stem>_bare.pdf/.png without.

    captions     artists to hide in the bare copy (fig.text handles, suptitle, ax titles).
                 None entries are skipped, so callers need not guard optional artists.
    bare_size    (w, h) inches for the bare copy. Hidden text leaves its whitespace behind
                 unless the canvas shrinks with it.
    bare_adjust  dict passed to subplots_adjust for the bare copy, to reclaim the margins
                 the captions occupied.
    """
    p = os.path.join(out_dir, stem)
    fig.savefig(p + ".pdf")
    fig.savefig(p + ".png", dpi=dpi)

    for a in captions:
        if a is not None:
            a.set_visible(False)
    if bare_size:
        fig.set_size_inches(*bare_size)
    if bare_adjust:
        fig.subplots_adjust(**bare_adjust)
    fig.savefig(p + "_bare.pdf")
    fig.savefig(p + "_bare.png", dpi=dpi)
    if close:
        plt.close(fig)
    if verbose:
        print(f"      wrote {stem}.pdf and {stem}_bare.pdf")
