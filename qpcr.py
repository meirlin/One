#!/usr/bin/env python3
"""
qpcr.py — single-file LinRegPCR-style qPCR analysis pipeline for Bio-Rad CFX.

This is the whole toolchain in one file (reader, efficiency, quantification,
inter-run calibration, statistics, plots, config-driven batch mode). Run it with
--help for usage. See the README for the full description.
"""
from __future__ import annotations

import argparse
import csv
import difflib
import glob
import math
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ===========================================================================
# ==== qpcr_io.py =====================================================
# ===========================================================================

"""
qpcr_io.py
==========
Robust readers for Bio-Rad CFX Maestro / CFX Manager Excel exports.

CFX writes a slightly non-standard .xlsx (paths use back-slashes,
`[content_types].xml` is lower-case), so openpyxl/pandas often fail on it.
These functions read the workbook directly from the zip container and are
tolerant of that quirk.

Two kinds of file are handled:

* "Quantification Amplification Results" -> one column per well, one row per
  cycle, of fluorescence values. Two versions usually exist per run:
      - RAW (background-corrected, NOT baseline-corrected): early cycles sit on
        a high, roughly flat plateau (e.g. ~2400 a.u.). This is what LinRegPCR
        needs.
      - BASELINE-SUBTRACTED: early cycles fluctuate around 0.
  `read_amplification` auto-detects which is which.

* "Quantification Cq Results" -> one row per well with Well / Fluor / Target /
  Sample / Cq annotation columns.
"""


import re
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Any

_MAIN_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


# --------------------------------------------------------------------------- #
# low-level xlsx helpers
# --------------------------------------------------------------------------- #
def _col_letters_to_index(col: str) -> int:
    idx = 0
    for ch in col:
        if ch.isalpha():
            idx = idx * 26 + (ord(ch.upper()) - 64)
    return idx - 1


def _zip_key_map(z: zipfile.ZipFile) -> dict[str, str]:
    """Map normalised (lower-case, forward-slash) member names -> real names."""
    return {n.replace("\\", "/").lower(): n for n in z.namelist()}


def _read_shared_strings(z: zipfile.ZipFile, keymap: dict[str, str]) -> list[str]:
    shared: list[str] = []
    hits = [k for k in keymap if k.endswith("sharedstrings.xml")]
    if not hits:
        return shared
    root = ET.fromstring(z.read(keymap[hits[0]]))
    for si in root:
        text = "".join(t.text or "" for t in si.iter(f"{_MAIN_NS}t"))
        shared.append(text)
    return shared


def _sheet_members(z: zipfile.ZipFile, keymap: dict[str, str]) -> list[str]:
    """Return worksheet member keys, ordered sheet1, sheet2, ..."""
    hits = [k for k in keymap if "/worksheets/" in k and k.endswith(".xml")]

    def order(k: str) -> int:
        m = re.search(r"sheet(\d+)\.xml$", k)
        return int(m.group(1)) if m else 999

    return sorted(hits, key=order)


def read_xlsx_sheets(path: str) -> list[list[list[Any]]]:
    """Read every worksheet of a (possibly non-standard) xlsx into a list of
    row-major tables. Cell values are str for shared strings, str for numbers
    (kept as text; callers convert)."""
    z = zipfile.ZipFile(path)
    keymap = _zip_key_map(z)
    shared = _read_shared_strings(z, keymap)
    tables: list[list[list[Any]]] = []
    for member in _sheet_members(z, keymap):
        root = ET.fromstring(z.read(keymap[member]))
        rows: list[list[Any]] = []
        for row in root.iter(f"{_MAIN_NS}row"):
            cells: dict[int, Any] = {}
            max_c = -1
            for c in row.iter(f"{_MAIN_NS}c"):
                ref = c.get("r") or ""
                col_letters = "".join(ch for ch in ref if ch.isalpha())
                ci = _col_letters_to_index(col_letters) if col_letters else (max_c + 1)
                t = c.get("t")
                v = c.find(f"{_MAIN_NS}v")
                inline = c.find(f"{_MAIN_NS}is")
                val: Any = None
                if v is not None and v.text is not None:
                    if t == "s":
                        try:
                            val = shared[int(v.text)]
                        except (ValueError, IndexError):
                            val = v.text
                    else:
                        val = v.text
                elif inline is not None:
                    val = "".join(t.text or "" for t in inline.iter(f"{_MAIN_NS}t"))
                cells[ci] = val
                max_c = max(max_c, ci)
            rows.append([cells.get(i) for i in range(max_c + 1)])
        tables.append(rows)
    z.close()
    return tables


def _to_float(x: Any):
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip()
    if s == "" or s.lower() in {"nan", "na", "n/a", "-", "null", "none"}:
        return None
    s = s.replace(",", ".")  # tolerate European decimal comma
    try:
        return float(s)
    except ValueError:
        return None


# --------------------------------------------------------------------------- #
# amplification file
# --------------------------------------------------------------------------- #
@dataclass
class AmplificationData:
    cycles: list[float]
    wells: dict[str, list[float]]          # well id (e.g. 'A1') -> fluorescence
    kind: str                              # 'raw' | 'baseline_subtracted'
    source_path: str
    sheet_name: str = ""


def _normalise_well(w: str) -> str:
    """'A01' / 'a1' / 'A 1' -> 'A1' (letter + int)."""
    s = str(w).strip().upper().replace(" ", "")
    m = re.match(r"^([A-H])0*([0-9]{1,2})$", s)
    if m:
        return f"{m.group(1)}{int(m.group(2))}"
    return s


def read_amplification(path: str) -> AmplificationData:
    """Read a CFX 'Quantification Amplification Results' workbook.

    Layout: first sheet, header row = ['', 'Cycle', 'A1', 'A2', ...], then one
    row per cycle. Detects raw vs baseline-subtracted from the early-cycle level.
    """
    tables = read_xlsx_sheets(path)
    # pick the sheet that has a 'Cycle' header
    table = None
    for t in tables:
        if not t:
            continue
        header = [str(c).strip().lower() if c is not None else "" for c in t[0]]
        if "cycle" in header:
            table = t
            break
    if table is None:
        table = tables[0]

    header = [str(c).strip() if c is not None else "" for c in table[0]]
    lower = [h.lower() for h in header]
    cyc_col = lower.index("cycle") if "cycle" in lower else 0
    well_cols = [
        (i, header[i])
        for i in range(len(header))
        if i != cyc_col and header[i] != "" and re.match(r"^[A-Ha-h]\s?0*\d{1,2}$", header[i])
    ]

    cycles: list[float] = []
    wells: dict[str, list[float]] = {_normalise_well(name): [] for _, name in well_cols}
    for row in table[1:]:
        if cyc_col >= len(row):
            continue
        cyc = _to_float(row[cyc_col])
        if cyc is None:
            continue
        cycles.append(cyc)
        for i, name in well_cols:
            val = _to_float(row[i]) if i < len(row) else None
            wells[_normalise_well(name)].append(val if val is not None else float("nan"))

    kind = _detect_kind(wells)
    return AmplificationData(cycles=cycles, wells=wells, kind=kind, source_path=path)


def _detect_kind(wells: dict[str, list[float]]) -> str:
    """Raw data has a large, roughly constant positive offset in the first few
    cycles; baseline-subtracted data hovers near zero."""
    import statistics as _st

    early_levels = []
    spans = []
    for series in wells.values():
        vals = [v for v in series if v == v]  # drop NaN
        if len(vals) < 8:
            continue
        early = vals[:5]
        early_levels.append(_st.median(early))
        spans.append(max(vals) - min(vals))
    if not early_levels:
        return "unknown"
    med_early = _st.median(early_levels)
    med_span = _st.median(spans) if spans else 1.0
    # raw: early level is a big fraction of the total span (a real baseline offset)
    if med_span > 0 and abs(med_early) > 0.25 * med_span:
        return "raw"
    return "baseline_subtracted"


# --------------------------------------------------------------------------- #
# Cq / annotation file
# --------------------------------------------------------------------------- #
@dataclass
class WellAnnotation:
    well: str
    target: str = ""
    sample: str = ""
    fluor: str = ""
    content: str = ""
    cq: float | None = None
    extra: dict[str, Any] = field(default_factory=dict)


def _find_header_row(table: list[list[Any]], needed: list[str], max_scan: int = 30) -> int:
    needed_low = [n.lower() for n in needed]
    for i, row in enumerate(table[:max_scan]):
        low = [str(c).strip().lower() if c is not None else "" for c in row]
        if all(any(n == cell for cell in low) for n in needed_low):
            return i
    # fall back: any row containing 'well'
    for i, row in enumerate(table[:max_scan]):
        low = [str(c).strip().lower() if c is not None else "" for c in row]
        if "well" in low:
            return i
    return 0


def read_cq(path: str) -> list[WellAnnotation]:
    """Read a CFX 'Quantification Cq Results' workbook into annotations."""
    tables = read_xlsx_sheets(path)
    table = None
    for t in tables:
        if not t:
            continue
        flat = [str(c).strip().lower() for row in t[:5] for c in row if c is not None]
        if "well" in flat and ("target" in flat or "cq" in flat):
            table = t
            break
    if table is None:
        table = tables[0]

    hrow = _find_header_row(table, ["well"])
    header = [str(c).strip() if c is not None else "" for c in table[hrow]]
    low = [h.lower() for h in header]

    def col(*names: str):
        for nm in names:
            if nm.lower() in low:
                return low.index(nm.lower())
        return None

    c_well = col("well")
    c_target = col("target")
    c_sample = col("sample")
    c_fluor = col("fluor")
    c_content = col("content")
    c_cq = col("cq")

    anns: list[WellAnnotation] = []
    for row in table[hrow + 1:]:
        if c_well is None or c_well >= len(row) or row[c_well] in (None, ""):
            continue
        well = _normalise_well(row[c_well])
        if not re.match(r"^[A-H]\d{1,2}$", well):
            continue
        ann = WellAnnotation(
            well=well,
            target=str(row[c_target]).strip() if c_target is not None and c_target < len(row) and row[c_target] is not None else "",
            sample=str(row[c_sample]).strip() if c_sample is not None and c_sample < len(row) and row[c_sample] is not None else "",
            fluor=str(row[c_fluor]).strip() if c_fluor is not None and c_fluor < len(row) and row[c_fluor] is not None else "",
            content=str(row[c_content]).strip() if c_content is not None and c_content < len(row) and row[c_content] is not None else "",
            cq=_to_float(row[c_cq]) if c_cq is not None and c_cq < len(row) else None,
        )
        anns.append(ann)
    return anns

# ===========================================================================
# ==== linreg.py ======================================================
# ===========================================================================

"""
linreg.py
=========
A faithful, self-contained re-implementation of the LinRegPCR amplification-
efficiency algorithm (Ramakers et al. 2003; Ruijter et al., Nucleic Acids Res.
2009, 37:e45).

Per well (raw, non-baseline-corrected fluorescence):
  1. Ground noise from the earliest cycles.
  2. Exponential phase located from the 2nd derivative (take-off .. plateau).
  3. Fluorescence baseline estimated by a constrained search near the ground
     level: the baseline that makes log10(F - baseline) straightest over a
     window in the exponential phase (Ruijter 2009 idea), with the window gated
     to the lower/middle exponential (well above noise, well below plateau).
  4. Individual window-of-linearity by maximum R^2; E_ind = 10^slope.

Per amplicon (all wells of one target):
  5. A COMMON window-of-linearity is chosen in log-fluorescence space as the band
     that MINIMISES the coefficient of variation of per-well efficiencies (this is
     the defining LinRegPCR criterion). Each well is re-fitted there -> E_common.
  6. Mean PCR efficiency per amplicon = mean of E_common over wells that amplify,
     are not baseline errors / noisy, and are not efficiency outliers.

Quality flags mirror LinRegPCR: no_amplification, baseline_error, no_plateau,
noisy, efficiency_outlier.
"""


from dataclasses import dataclass, field

import numpy as np


# --------------------------------------------------------------------------- #
# numerical helpers
# --------------------------------------------------------------------------- #
def _linfit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    """OLS. Returns slope, intercept, R^2, residual SD."""
    n = len(x)
    if n < 2:
        return float("nan"), float("nan"), float("nan"), float("nan")
    x = x.astype(float)
    y = y.astype(float)
    xm, ym = x.mean(), y.mean()
    sxx = np.sum((x - xm) ** 2)
    if sxx <= 0:
        return float("nan"), float("nan"), float("nan"), float("nan")
    slope = np.sum((x - xm) * (y - ym)) / sxx
    intercept = ym - slope * xm
    yhat = slope * x + intercept
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - ym) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    resid_sd = np.sqrt(ss_res / (n - 2)) if n > 2 else 0.0
    return float(slope), float(intercept), float(r2), float(resid_sd)


def _smooth(y: np.ndarray, k: int = 3) -> np.ndarray:
    if k <= 1 or len(y) < k:
        return y
    return np.convolve(y, np.ones(k) / k, mode="same")


# --------------------------------------------------------------------------- #
# result container
# --------------------------------------------------------------------------- #
@dataclass
class WellFit:
    well: str
    baseline: float = float("nan")
    baseline_slope: float = float("nan")
    ground_noise: float = float("nan")
    plateau_fluor: float = float("nan")
    amplitude: float = float("nan")
    fold_increase: float = float("nan")
    takeoff_cycle: float = float("nan")
    plateau_cycle: float = float("nan")
    # individual window
    ind_cycle_start: float = float("nan")
    ind_cycle_end: float = float("nan")
    ind_n_points: int = 0
    ind_slope: float = float("nan")
    ind_intercept: float = float("nan")
    ind_r2: float = float("nan")
    ind_E: float = float("nan")
    # common window
    com_cycle_start: float = float("nan")
    com_cycle_end: float = float("nan")
    com_n_points: int = 0
    com_slope: float = float("nan")
    com_intercept: float = float("nan")
    com_r2: float = float("nan")
    com_E: float = float("nan")
    # flags
    flags: list[str] = field(default_factory=list)
    usable: bool = False
    # per-cycle diagnostics
    cycles: np.ndarray | None = None
    raw: np.ndarray | None = None
    corrected: np.ndarray | None = None

    @property
    def E(self) -> float:
        return self.com_E if np.isfinite(self.com_E) else self.ind_E

    @property
    def quality(self) -> str:
        if "no_amplification" in self.flags:
            return "no_amplification"
        if "baseline_error" in self.flags:
            return "baseline_error"
        if "noisy" in self.flags:
            return "noisy"
        if "no_plateau" in self.flags:
            return "no_plateau"
        if "efficiency_outlier" in self.flags:
            return "efficiency_outlier"
        return "OK"


# --------------------------------------------------------------------------- #
# phase detection
# --------------------------------------------------------------------------- #
def _detect_phase(cycles: np.ndarray, raw: np.ndarray) -> tuple[int, int]:
    """Return (takeoff_idx, plateau_idx) bounding the exponential phase."""
    n = len(raw)
    s = _smooth(raw, 3)
    rng = np.nanmax(s) - np.nanmin(s)
    if rng <= 0 or not np.isfinite(rng):
        return 1, n - 1
    d1 = np.gradient(s)
    d2 = np.gradient(d1)
    infl = int(np.nanargmax(d1))
    plateau_idx = infl + int(np.nanargmin(d2[infl:])) if infl < n - 2 else n - 1
    takeoff_idx = int(np.nanargmax(d2[: infl + 1])) if infl >= 2 else 0

    base = np.nanmedian(s[: max(3, n // 8)])
    top = np.nanmax(s)
    if top - base > 0:
        frac = (s - base) / (top - base)
        a5 = np.where(frac > 0.05)[0]
        a80 = np.where(frac > 0.80)[0]
        if len(a5):
            takeoff_idx = max(takeoff_idx, int(a5[0]) - 1)
        if len(a80):
            plateau_idx = min(plateau_idx if plateau_idx > 0 else n - 1, int(a80[0]))
    takeoff_idx = max(1, min(takeoff_idx, n - 3))
    plateau_idx = max(takeoff_idx + 2, min(plateau_idx, n - 1))
    return takeoff_idx, plateau_idx


# --------------------------------------------------------------------------- #
# baseline estimation (constrained near ground, gated exponential window)
# --------------------------------------------------------------------------- #
def _estimate_baseline(
    cycles: np.ndarray, raw: np.ndarray, ground: float, noise: float, amp: float,
    w: int, takeoff_idx: int
) -> tuple[np.ndarray, float, float]:
    """Per-cycle baseline modelling linear instrument drift, fitted to the flat
    ground region before take-off (mirrors the CFX baseline; avoids the LinReg
    pitfall where an R^2-maximising offset drifts the baseline upward and
    underestimates efficiency). Returns (baseline_array, baseline_at_takeoff,
    slope_per_cycle)."""
    # linear drift fitted to the ground region (cycles before take-off).
    gr = max(takeoff_idx - 1, 4)
    gr = min(gr, len(cycles) - 1)
    if gr >= 3:
        slope_b, intercept_b, _, _ = _linfit(cycles[:gr], raw[:gr])
        if not np.isfinite(slope_b) or abs(slope_b) > 0.05 * amp:
            slope_b, intercept_b = 0.0, ground
    else:
        slope_b, intercept_b = 0.0, ground
    baseline_array = intercept_b + slope_b * cycles
    b_at_takeoff = float(baseline_array[min(takeoff_idx, len(baseline_array) - 1)])
    return baseline_array, b_at_takeoff, float(slope_b)


# --------------------------------------------------------------------------- #
# per-well fit
# --------------------------------------------------------------------------- #
def fit_well(
    well: str,
    cycles_in: list[float],
    raw_in: list[float],
    *,
    already_baseline_subtracted: bool = False,
    window: int = 4,
    window_max: int = 6,
    min_fold: float = 7.0,
    noise_cycles: int = 5,
) -> WellFit:
    cycles = np.asarray(cycles_in, dtype=float)
    raw = np.asarray(raw_in, dtype=float)
    m = np.isfinite(cycles) & np.isfinite(raw)
    cycles, raw = cycles[m], raw[m]
    fit = WellFit(well=well, cycles=cycles, raw=raw)

    n = len(raw)
    if n < window + 3:
        fit.flags.append("too_few_points")
        return fit

    k = min(noise_cycles, max(3, n // 8))
    ground = float(np.median(raw[:k]))
    noise = float(np.std(raw[:k], ddof=1)) if k > 2 else float(np.std(raw[:k]))
    noise = max(noise, 1e-9)
    fit.ground_noise = noise

    takeoff, plateau = _detect_phase(cycles, raw)
    fit.takeoff_cycle = float(cycles[takeoff])
    fit.plateau_cycle = float(cycles[plateau])
    plateau_level = float(np.median(raw[-3:]))
    fit.plateau_fluor = float(np.nanmax(raw))
    amp = max(plateau_level - ground, 1e-9)
    fit.amplitude = amp

    if already_baseline_subtracted:
        baseline_array = np.zeros_like(raw)
        fit.baseline = 0.0
        fit.baseline_slope = 0.0
    else:
        baseline_array, b_scalar, b_slope = _estimate_baseline(
            cycles, raw, ground, noise, amp, window + 1, takeoff)
        fit.baseline = b_scalar
        fit.baseline_slope = b_slope

    corrected = raw - baseline_array
    fit.corrected = corrected

    fold = float(np.nanmax(corrected) / max(noise, 1e-9))
    fit.fold_increase = fold
    if np.nanmax(corrected) < min_fold * noise:
        fit.flags.append("no_amplification")

    # exponential region for the individual window (lower/middle exponential)
    lo_gate = max(3.0 * noise, 0.02 * amp)
    hi_gate = 0.50 * amp
    idx = np.where((corrected > lo_gate) & (corrected < hi_gate))[0]
    # keep the longest consecutive run
    if len(idx):
        runs, cur = [], [idx[0]]
        for j in idx[1:]:
            if j == cur[-1] + 1:
                cur.append(j)
            else:
                runs.append(cur)
                cur = [j]
        runs.append(cur)
        idx = np.array(max(runs, key=len))
    if len(idx) < window:
        if "no_amplification" not in fit.flags:
            fit.flags.append("baseline_error")
        return fit

    reg_c = cycles[idx]
    reg_f = corrected[idx]
    log_f = np.log10(reg_f)

    # Window selection: prefer the LOWEST straight window (closest to take-off),
    # because the true PCR efficiency is expressed at the bottom of the
    # exponential phase; higher windows are progressively biased towards the
    # plateau. Among widths window..window_max, take the lowest-positioned window
    # whose R^2 >= r2_ok; if none is that straight, take the highest-R^2 window.
    r2_ok = 0.990
    lowest = None       # (start_fluor, w, xs, ys, slope, intercept, r2)
    best_r2 = None      # fallback: highest R^2
    for w in range(window, window_max + 1):
        if len(log_f) < w:
            continue
        for s in range(0, len(log_f) - w + 1):
            xs = reg_c[s : s + w]
            ys = log_f[s : s + w]
            slope, intercept, r2, _ = _linfit(xs, ys)
            if slope <= 0 or not np.isfinite(r2):
                continue
            if best_r2 is None or r2 > best_r2[6]:
                best_r2 = (float(reg_f[s]), w, xs, ys, slope, intercept, r2)
            if r2 >= r2_ok:
                # position key = starting fluorescence (lower is better); tie-break wider window
                if lowest is None or (reg_f[s], -w) < (lowest[0], -lowest[1]):
                    lowest = (float(reg_f[s]), w, xs, ys, slope, intercept, r2)
    best = lowest if lowest is not None else best_r2
    if best is None:
        fit.flags.append("baseline_error")
        return fit

    _, w, xs, ys, slope, intercept, r2 = best
    # noisy flag: poor straightness of the best available window
    if r2 < 0.97:
        fit.flags.append("noisy")
    fit.ind_cycle_start = float(xs[0])
    fit.ind_cycle_end = float(xs[-1])
    fit.ind_n_points = int(w)
    fit.ind_slope = float(slope)
    fit.ind_intercept = float(intercept)
    fit.ind_r2 = float(r2)
    fit.ind_E = float(10 ** slope)

    if plateau >= n - 1 and plateau_level < 0.9 * np.nanmax(raw):
        fit.flags.append("no_plateau")

    fit.usable = "no_amplification" not in fit.flags and "baseline_error" not in fit.flags
    return fit


# --------------------------------------------------------------------------- #
# common window per amplicon (median low-exponential band; robust)
# --------------------------------------------------------------------------- #
@dataclass
class AmpliconWindow:
    target: str
    log_lower: float = float("nan")
    log_upper: float = float("nan")
    fluor_lower: float = float("nan")
    fluor_upper: float = float("nan")
    mean_E: float = float("nan")
    sd_E: float = float("nan")
    cv_E: float = float("nan")
    n_used: int = 0
    n_total: int = 0
    median_cycle_start: float = float("nan")
    median_cycle_end: float = float("nan")
    threshold_fluor: float = float("nan")


def _refit(fit: "WellFit", Llo: float, Lhi: float, min_pts: int) -> float | None:
    fc = fit.corrected
    if fc is None:
        return None
    m = fc > 0
    idx_all = np.where(m)[0]
    lf_all = np.log10(fc[m])
    in_band = (lf_all >= Llo) & (lf_all <= Lhi)
    sel_idx = idx_all[in_band]
    if len(sel_idx) < min_pts:
        return None
    # split into consecutive-cycle runs; a clean exponential crosses the band once
    runs, cur = [], [sel_idx[0]]
    for j in sel_idx[1:]:
        if j == cur[-1] + 1:
            cur.append(j)
        else:
            runs.append(cur)
            cur = [j]
    runs.append(cur)
    # choose the run giving the best positive-slope log-linear fit (>= min_pts)
    best = None
    for run in runs:
        if len(run) < min_pts:
            continue
        c = fit.cycles[run]
        y = np.log10(fc[run])
        slope, intercept, r2, _ = _linfit(c, y)
        if slope <= 0 or not np.isfinite(slope):
            continue
        key = (r2, len(run))
        if best is None or key > best[0]:
            best = (key, c, slope, intercept, r2)
    if best is None:
        return None
    _, c, slope, intercept, r2 = best
    fit.com_cycle_start = float(c[0])
    fit.com_cycle_end = float(c[-1])
    fit.com_n_points = int(len(c))
    fit.com_slope = float(slope)
    fit.com_intercept = float(intercept)
    fit.com_r2 = float(r2)
    fit.com_E = float(10 ** slope)
    return fit.com_E


def set_common_window(
    target: str,
    fits: list[WellFit],
    *,
    min_pts: int = 4,
    outlier_halfwidth: float = 0.10,
) -> AmpliconWindow:
    """Define the amplicon common window-of-linearity as the median
    low-exponential log-fluorescence band across usable wells, then re-fit every
    well within it so that all wells report an efficiency measured over the same
    fluorescence interval (LinRegPCR common window)."""
    aw = AmpliconWindow(target=target, n_total=len(fits))
    # Band search uses every well that has a valid corrected curve and amplified
    # (not only the strictly "usable" ones): wells that cannot be fitted inside a
    # candidate band drop out naturally, exactly as in LinRegPCR. Restricting to
    # "usable" here would let a narrow near-plateau band win on a small subset.
    good = [
        f for f in fits
        if f.corrected is not None
        and "no_amplification" not in f.flags
        and np.count_nonzero(f.corrected > 0) >= 5
    ]
    if len(good) < 3:
        for f in fits:
            f.com_E = f.ind_E
            f.com_cycle_start, f.com_cycle_end = f.ind_cycle_start, f.ind_cycle_end
            f.com_n_points, f.com_slope = f.ind_n_points, f.ind_slope
            f.com_intercept, f.com_r2 = f.ind_intercept, f.ind_r2
        es = [f.ind_E for f in good]
        if es:
            aw.mean_E = float(np.mean(es))
            aw.sd_E = float(np.std(es, ddof=1)) if len(es) > 1 else 0.0
            aw.n_used = len(es)
        return aw

    # ----------------------------------------------------------------- #
    # Common window-of-linearity (LinRegPCR criterion): the shared log-
    # fluorescence band [L, U] that MINIMISES the coefficient of variation of
    # the per-well efficiencies. The search runs over an upper edge U (kept
    # below the plateau) and a band width; L = U - width. For each candidate
    # band every well is fitted on the points whose corrected log-fluorescence
    # falls in [L, U]; wells within median +/- outlier_halfwidth define the CV.
    # ----------------------------------------------------------------- #
    def _fit_band(Llo: float, Lhi: float) -> list[tuple["WellFit", float]]:
        res = []
        for f in good:
            fc = f.corrected
            m = fc > 0
            lf = np.log10(fc[m])
            cc = f.cycles[m]
            sel = (lf >= Llo) & (lf <= Lhi)
            if np.count_nonzero(sel) < min_pts:
                continue
            slope, _, r2, _ = _linfit(cc[sel], lf[sel])
            if slope > 0 and np.isfinite(slope) and r2 > 0.90:
                res.append((f, float(10 ** slope)))
        return res

    # per-well log-fluorescence anchors: plateau (0.30 * max) and noise floor
    plateau_logs, noise_logs = [], []
    for f in good:
        fc = f.corrected
        if np.count_nonzero(fc > 0) < 5:
            continue
        mx = float(np.nanmax(fc))
        if mx <= 0:
            continue
        plateau_logs.append(np.log10(0.30 * mx))
        floor = max(3.0 * f.ground_noise, 0.03 * mx, 1.0)
        noise_logs.append(np.log10(floor))
    if not plateau_logs:
        return aw
    U0 = float(np.median(plateau_logs))
    L0 = float(np.median(noise_logs))

    n_need_fit = max(6, len(fits) // 3)
    # Keep the upper edge at or below ~0.30 * plateau (U0): above that the curve
    # bends into the plateau and per-well efficiencies read systematically low, so
    # a pure minimum-CV search would be biased toward a tight-but-low near-plateau
    # band. Among admissible exponential-phase bands we minimise CV, with a small
    # reward for bands supported by more wells (more robust common window).
    best = None  # (score, cv, L, U, meanE, sdE, n_kept, n_fit)
    for U in np.linspace(L0 + 0.5, U0, 16):
        for width in np.linspace(0.5, 1.3, 12):
            Llo = U - width
            res = _fit_band(Llo, U)
            if len(res) < n_need_fit:
                continue
            es = np.array([e for _, e in res])
            med = float(np.median(es))
            keep = es[np.abs(es - med) <= outlier_halfwidth]
            if len(keep) < 5:
                continue
            m_e = float(keep.mean())
            cv = float(keep.std(ddof=1) / m_e) if m_e > 0 else float("inf")
            # objective: low CV, mild reward for more supporting wells
            score = cv - 0.0015 * len(keep)
            if best is None or score < best[0]:
                best = (score, cv, float(Llo), float(U), m_e,
                        float(keep.std(ddof=1)), int(len(keep)), int(len(res)))
    if best is not None:  # drop the score field for the downstream unpack
        best = best[1:]

    # Fallback: if no band satisfied the quorum, relax the kept-count / fit-count
    if best is None:
        for U in np.linspace(L0 + 0.4, U0 + 0.3, 20):
            for width in np.linspace(0.4, 1.4, 14):
                Llo = U - width
                res = _fit_band(Llo, U)
                if len(res) < 3:
                    continue
                es = np.array([e for _, e in res])
                med = float(np.median(es))
                keep = es[np.abs(es - med) <= outlier_halfwidth]
                if len(keep) < 2:
                    keep = es
                m_e = float(keep.mean())
                cv = float(keep.std(ddof=1) / m_e) if (m_e > 0 and len(keep) > 1) else float("inf")
                if best is None or cv < best[0]:
                    best = (cv, float(Llo), float(U), m_e,
                            float(np.std(keep, ddof=1)) if len(keep) > 1 else 0.0,
                            int(len(keep)), int(len(res)))
    if best is None:
        # last resort: individual efficiencies
        for f in fits:
            f.com_E = f.ind_E
            f.com_cycle_start, f.com_cycle_end = f.ind_cycle_start, f.ind_cycle_end
            f.com_n_points, f.com_slope = f.ind_n_points, f.ind_slope
            f.com_intercept, f.com_r2 = f.ind_intercept, f.ind_r2
        es = [f.ind_E for f in good]
        aw.mean_E = float(np.mean(es))
        aw.sd_E = float(np.std(es, ddof=1)) if len(es) > 1 else 0.0
        aw.n_used = len(es)
        return aw

    _, L, U, _m_e, _sd_e, _nk, _nf = best
    aw.log_lower, aw.log_upper = float(L), float(U)
    aw.fluor_lower, aw.fluor_upper = float(10 ** L), float(10 ** U)
    aw.threshold_fluor = float(10 ** U)

    # Re-fit EVERY well in the chosen band so all report E over the same interval.
    for f in fits:
        _refit(f, L, U, min_pts)
        if not np.isfinite(f.com_E):  # band misses this well -> keep individual
            f.com_E = f.ind_E
            f.com_cycle_start, f.com_cycle_end = f.ind_cycle_start, f.ind_cycle_end
            f.com_n_points, f.com_slope = f.ind_n_points, f.ind_slope
            f.com_intercept, f.com_r2 = f.ind_intercept, f.ind_r2

    starts = [f.com_cycle_start for f in fits if np.isfinite(f.com_cycle_start)]
    ends = [f.com_cycle_end for f in fits if np.isfinite(f.com_cycle_end)]
    aw.median_cycle_start = float(np.median(starts)) if starts else float("nan")
    aw.median_cycle_end = float(np.median(ends)) if ends else float("nan")
    return aw


def finalise_amplicon(
    target: str, fits: list[WellFit], aw: AmpliconWindow, *, outlier_halfwidth: float = 0.10
) -> AmpliconWindow:
    """Mean PCR efficiency per amplicon from common-window E values, excluding
    no-amplification / baseline-error / noisy wells and efficiency outliers
    (default: outside median +/- 0.10)."""
    cand = [
        f for f in fits
        if f.usable and "noisy" not in f.flags and np.isfinite(f.com_E) and f.com_E > 1.0
    ]
    if not cand:
        cand = [f for f in fits if np.isfinite(f.com_E) and f.com_E > 1.0]
    if not cand:
        return aw
    es = np.array([f.com_E for f in cand])
    med = float(np.median(es))
    keep = np.abs(es - med) <= outlier_halfwidth
    for f, kk in zip(cand, keep):
        if not kk and "efficiency_outlier" not in f.flags:
            f.flags.append("efficiency_outlier")
    kept = es[keep] if keep.any() else es
    aw.mean_E = float(np.mean(kept))
    aw.sd_E = float(np.std(kept, ddof=1)) if len(kept) > 1 else 0.0
    aw.cv_E = aw.sd_E / aw.mean_E if aw.mean_E > 0 else float("nan")
    aw.n_used = int(len(kept))
    return aw

# ===========================================================================
# ==== metadata (ONE table) + quantification + stats + CLI ==================
# ===========================================================================
#
# Two files drive everything:
#
#   manifest.csv     Plate, Amp, Cq                      (unchanged)
#   sample_map.csv   ONE table with two kinds of row:
#       * SAMPLE rows  -> Plate, Sample, Group, Replicate, Type, Control
#       * GENE rows    -> Gene, Role, Calibrate  (Sample blank; Plate optional)
#
# Everything experiment-specific lives in these tables; NO gene or sample name
# is ever written inside the code, so the same code runs for any experiment.
#
# Rules:
#   * Only genes listed as GENE rows are computed. Anything not listed (e.g. a
#     gene you don't want, like an unused gene) is not computed at all.
#   * Role = target | reference. Any number of reference genes is allowed; the
#     normalisation factor is their GEOMETRIC MEAN, per biological sample.
#   * Calibrate = yes marks a gene for inter-run calibration; it is calibrated
#     only if it is on >=2 plates and a Type=calibrator sample carries it there.
#     Genes without Calibrate=yes are NEVER calibrated (nothing is automatic).
#   * Control group = the Group whose SAMPLE rows have Control=yes.
#   * Samples are matched across plates by (Group, Replicate), so inconsistent
#     sample names between plates do not matter.
#
import argparse
import csv
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np


def norm(name) -> str:
    import re
    s = str(name).strip().strip("'\"").lower()
    return re.sub(r"\s+", " ", s)


def norm_gene(name) -> str:
    import re
    return re.sub(r"[\s\-_.]+", "", norm(name))


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"1", "yes", "y", "true", "t", "on", "x", "+"}


def _sniff_delim(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        head = f.readline()
    if head.count("\t") > head.count(","):
        return "\t"
    if head.count(";") > head.count(",") and head.count(",") == 0:
        return ";"
    return ","


def read_table_rows(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xlsm", ".xltx"):
        try:
            import openpyxl
        except ModuleNotFoundError:
            raise SystemExit(f"{path} is .xlsx but openpyxl is missing; save as .csv or pip install openpyxl.")
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws = wb.active
        rows = [["" if c is None else str(c).strip() for c in r] for r in ws.iter_rows(values_only=True)]
        wb.close()
        return [r for r in rows if any(c for c in r)]
    delim = _sniff_delim(path)
    out = []
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        for r in csv.reader(f, delimiter=delim):
            cells = [("" if c is None else str(c).strip()) for c in r]
            if any(cells):
                out.append(cells)
    return out


def _header_index(header, aliases):
    import re
    def n(h): return re.sub(r"[\s\-_.]+", "", str(h).strip().lower())
    nh = [n(h) for h in header]
    idx = {}
    for canon, names in aliases.items():
        want = {n(a) for a in names}
        for i, h in enumerate(nh):
            if h in want:
                idx[canon] = i
                break
    return idx


# --------------------------------------------------------------------------- #
# dataclasses
# --------------------------------------------------------------------------- #
@dataclass
class PlateFiles:
    plate: str
    amp: list
    cq: str


@dataclass
class SampleRow:
    plate: str
    sample: str
    group: str
    replicate: str
    type: str = "sample"
    is_control: bool = False


_MANIFEST_ALIASES = {
    "plate": ["plate", "run", "batch", "id", "experiment"],
    "amp": ["amp", "amplification", "ampfile", "amplificationresults"],
    "cq": ["cq", "cqfile", "cqresults", "quantificationcq"],
}
_META_ALIASES = {
    "plate": ["plate", "run", "batch"],
    "sample": ["sample", "samplename", "name", "label"],
    "group": ["group", "condition", "treatment", "genotype", "cond"],
    "replicate": ["replicate", "rep", "biorep", "biologicalreplicate"],
    "type": ["type", "sampletype", "class", "kind"],
    "control": ["control", "iscontrol", "ctrl"],
    "gene": ["gene", "target", "amplicon", "assay"],
    "role": ["role", "generole"],
    "calibrate": ["calibrate", "calibration", "irc", "usecalibrator", "intercalibrate"],
}
_CAL_TYPES = {"calibrator", "irc", "cal"}
_NTC_TYPES = {"ntc", "notemplate", "blank", "water", "neg", "negative", "nrt", "standard", "std"}


def load_manifest(path):
    rows = read_table_rows(path)
    if not rows:
        raise SystemExit(f"{path}: empty manifest.")
    idx = _header_index(rows[0], _MANIFEST_ALIASES)
    for need in ("plate", "amp", "cq"):
        if need not in idx:
            raise SystemExit(f"{path}: a '{need}' column is required.")
    base = os.path.dirname(os.path.abspath(path))
    def resolve(p):
        p = p.strip()
        return "" if not p else (p if os.path.isabs(p) else os.path.normpath(os.path.join(base, p)))
    plates = {}
    for i, r in enumerate(rows[1:], start=2):
        def cell(k): return r[idx[k]] if idx[k] < len(r) else ""
        plate = cell("plate").strip()
        amp, cq = resolve(cell("amp")), resolve(cell("cq"))
        if not plate:
            continue
        if not amp or not cq:
            raise SystemExit(f"{path} row {i}: both Amp and Cq are required.")
        if plate in plates:
            plates[plate].amp.append(amp)
        else:
            plates[plate] = PlateFiles(plate=plate, amp=[amp], cq=cq)
    return list(plates.values())


def load_metadata(path):
    """Read the single sample_map (sample rows + gene rows).

    Returns: samples[(plate_norm, sample_norm)] -> SampleRow,
             control_groups (set of norm group),
             gene_role[gene_norm] -> 'target'|'reference',
             references (set of gene_norm),
             gene_display[gene_norm] -> display name."""
    rows = read_table_rows(path)
    if not rows:
        raise SystemExit(f"{path}: empty sample map.")
    idx = _header_index(rows[0], _META_ALIASES)
    if "sample" not in idx and "gene" not in idx:
        raise SystemExit(f"{path}: need at least a 'Sample' or 'Gene' column.")
    samples = {}
    control_groups = set()
    gene_role = {}
    gene_display = {}
    gene_calibrate = {}
    for r in rows[1:]:
        def cell(k): return r[idx[k]].strip() if k in idx and idx[k] < len(r) else ""
        gene = cell("gene")
        sample = cell("sample")
        if gene and not gene.startswith("#"):                     # GENE row
            gn = norm_gene(gene)
            role = norm(cell("role"))
            role = "reference" if role.startswith("ref") else "target"
            gene_role[gn] = role
            gene_display[gn] = gene
            gene_calibrate[gn] = _truthy(cell("calibrate"))       # explicit, per gene
            continue
        if not sample or sample.startswith("#"):
            continue
        plate = cell("plate")                                     # SAMPLE row
        typ = norm(cell("type")) or "sample"
        typ = "calibrator" if typ in _CAL_TYPES else "ntc" if typ in _NTC_TYPES else "sample"
        is_ctrl = _truthy(cell("control"))
        group = cell("group")
        if is_ctrl and group:
            control_groups.add(norm(group))
        samples[(norm(plate), norm(sample))] = SampleRow(
            plate=plate, sample=sample, group=group, replicate=cell("replicate"),
            type=typ, is_control=is_ctrl)
    references = {gn for gn, role in gene_role.items() if role == "reference"}
    if not gene_role:
        raise SystemExit(f"{path}: no GENE rows found (add rows with a Gene and Role).")
    return samples, control_groups, gene_role, references, gene_display, gene_calibrate


# --------------------------------------------------------------------------- #
# per-well record + per-plate processing
# --------------------------------------------------------------------------- #
@dataclass
class WellRecord:
    plate: str
    source_file: str
    well: str
    gene: str
    gene_norm: str
    role: str
    sample: str
    group: str
    replicate: str
    stype: str
    is_control: bool
    cq: float | None
    cq_source: str
    gene_mean_E: float
    ind_E: float
    com_E: float
    com_r2: float
    win_start: float
    win_end: float
    win_n: int
    quality: str
    flags: str
    used_for_mean: bool
    E_used: float
    N0: float = float("nan")
    N0_cal: float = float("nan")
    cf: float = 1.0


def load_exclude_wells(path):
    """Read a small table with Plate, Well columns -> set of (plate_norm, well_norm)."""
    if not path:
        return set()
    rows = read_table_rows(path)
    if not rows:
        return set()
    idx = _header_index(rows[0], {"plate": ["plate", "run"], "well": ["well", "position"]})
    if "plate" not in idx or "well" not in idx:
        raise SystemExit(f"{path}: needs Plate and Well columns.")
    out = set()
    for r in rows[1:]:
        pl = r[idx["plate"]].strip() if idx["plate"] < len(r) else ""
        wl = r[idx["well"]].strip() if idx["well"] < len(r) else ""
        if pl and wl and not pl.startswith("#"):
            out.add((norm(pl), _normalise_well(wl)))
    return out


def process_plate(pf, samples, gene_role, *, window=4, window_max=6, eff_outlier=0.10,
                  efficiency_mode="amplicon_mean", fixed_efficiency=2.0,
                  exclude_wells=None, force_baseline_subtracted=False):
    exclude_wells = exclude_wells or set()
    parsed = [(p, read_amplification(p)) for p in pf.amp]
    raws = [(p, d) for p, d in parsed if d.kind == "raw"]
    chosen_path, amp = (raws[0] if raws else parsed[0])
    already_bl = force_baseline_subtracted or amp.kind == "baseline_subtracted"
    source_file = os.path.basename(chosen_path)
    well2ann = {a.well: a for a in read_cq(pf.cq)}
    pkey = norm(pf.plate)

    wells_by_gene = defaultdict(list)
    display = {}
    for well in amp.wells:
        if (pkey, _normalise_well(well)) in exclude_wells:   # user-excluded well
            continue
        a = well2ann.get(well)
        if not a or not a.target:
            continue
        gn = norm_gene(a.target)
        if gn not in gene_role:                       # not declared -> skip entirely
            continue
        wells_by_gene[gn].append(well)
        display.setdefault(gn, a.target)

    records, summaries = [], []
    for gn, wells in sorted(wells_by_gene.items()):
        fits = [fit_well(w, amp.cycles, amp.wells[w], already_baseline_subtracted=already_bl,
                         window=window, window_max=window_max) for w in wells]
        aw = set_common_window(display[gn], fits)
        aw = finalise_amplicon(display[gn], fits, aw, outlier_halfwidth=eff_outlier)
        summaries.append({"plate": pf.plate, "source_file": source_file, "gene": display[gn],
                          "role": gene_role[gn], "mean_E": aw.mean_E,
                          "efficiency_percent": (aw.mean_E - 1) * 100 if math.isfinite(aw.mean_E) else float("nan"),
                          "sd_E": aw.sd_E, "n_wells_used": aw.n_used, "n_wells_total": aw.n_total})
        for f in fits:
            a = well2ann.get(f.well)
            sample = a.sample if a else ""
            sr = samples.get((pkey, norm(sample)))
            if sr:
                group, rep, stype, ctrl = sr.group, sr.replicate, sr.type, sr.is_control
            else:
                group, rep, stype, ctrl = "", "", "undeclared", False
            # which efficiency feeds N0 = E^(-Cq): gene mean / this well's own / fixed
            if efficiency_mode == "fixed":
                E = float(fixed_efficiency)
            elif efficiency_mode == "per_well":
                E = f.com_E if math.isfinite(f.com_E) and f.com_E > 1 else f.ind_E
            else:  # amplicon_mean (default)
                E = aw.mean_E
            N0 = float(E ** (-a.cq)) if (a and a.cq is not None and math.isfinite(a.cq)
                                         and math.isfinite(E) and E > 1) else float("nan")
            records.append(WellRecord(
                plate=pf.plate, source_file=source_file, well=f.well, gene=display[gn], gene_norm=gn,
                role=gene_role[gn], sample=sample, group=group, replicate=rep, stype=stype, is_control=ctrl,
                cq=(a.cq if a else None), cq_source=("file" if (a and a.cq is not None) else "none"),
                gene_mean_E=aw.mean_E, ind_E=f.ind_E, com_E=f.com_E, com_r2=f.com_r2,
                E_used=E,
                win_start=f.com_cycle_start, win_end=f.com_cycle_end, win_n=f.com_n_points,
                quality=f.quality, flags=";".join(f.flags),
                used_for_mean=(f.usable and "noisy" not in f.flags and "efficiency_outlier" not in f.flags),
                N0=N0, N0_cal=N0))
    return records, summaries, amp.kind


# --------------------------------------------------------------------------- #
# automatic per-gene inter-run calibration
# --------------------------------------------------------------------------- #
def _geomean(xs):
    xs = [x for x in xs if x and math.isfinite(x) and x > 0]
    return float(np.exp(np.mean(np.log(xs)))) if xs else float("nan")


def apply_calibration(records, calibrate_genes):
    """Inter-run calibrate ONLY the genes explicitly flagged Calibrate=yes
    (calibrate_genes = set of gene_norm). A flagged gene is calibrated only if
    it is on >=2 plates AND calibrator wells carry it on >=2 of them; otherwise
    it is left uncalibrated. Genes not flagged are never touched."""
    for r in records:
        r.N0_cal = r.N0
        r.cf = 1.0
    factors, warns = [], []
    by_gene = defaultdict(list)
    for r in records:
        by_gene[r.gene_norm].append(r)
    for gn, recs in by_gene.items():
        if gn not in calibrate_genes:                          # not flagged -> never calibrate
            continue
        plates = sorted({r.plate for r in recs})
        if len(plates) < 2:                                    # single plate -> nothing to reconcile
            warns.append(f"gene '{recs[0].gene}' flagged Calibrate but on 1 plate; left uncalibrated.")
            continue
        cal_by_plate = {}
        for p in plates:
            cal_by_plate[p] = _geomean([r.N0 for r in recs if r.plate == p and r.stype == "calibrator"])
        cal_by_plate = {p: v for p, v in cal_by_plate.items() if math.isfinite(v)}
        if len(cal_by_plate) < 2:                              # no calibrator -> no calibration
            warns.append(f"gene '{recs[0].gene}' flagged Calibrate but calibrator found on "
                         f"{len(cal_by_plate)}/{len(plates)} plates; left uncalibrated.")
            continue
        geo = _geomean(list(cal_by_plate.values()))
        for p in plates:
            cf = cal_by_plate.get(p)
            if not cf or not geo:
                continue
            factor = cf / geo
            factors.append({"gene": recs[0].gene, "plate": p, "calibrator_N0": cf,
                            "cross_plate_geomean": geo, "calibration_factor": factor,
                            "n_calibrator_wells": sum(1 for r in recs if r.plate == p and r.stype == "calibrator")})
            for r in recs:
                if r.plate == p and math.isfinite(r.N0):
                    r.cf = factor
                    r.N0_cal = r.N0 / factor
    return factors, warns


# --------------------------------------------------------------------------- #
# aggregate to (Group, Replicate) x gene -> NRQ, fold change
# --------------------------------------------------------------------------- #
@dataclass
class ExprValue:
    group: str
    replicate: str
    gene: str
    gene_norm: str
    role: str
    n_tech: int
    mean_cq: float
    mean_N0: float
    mean_N0_cal: float
    nrq: float = float("nan")
    fold: float = float("nan")
    log2_fold: float = float("nan")
    tech_cv: float = float("nan")     # CV of technical wells (QC)
    n_tech_flag: str = ""             # 'check' if well count != gene's usual
    outlier: str = ""                 # 'outlier' if biological replicate is a MAD outlier


def _mad_outliers(x):
    """Return boolean list: modified z-score (MAD-based) > 3.5. Needs n>=3."""
    x = np.asarray(x, float)
    if len(x) < 3:
        return [False] * len(x)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:
        return [False] * len(x)
    mz = 0.6745 * (x - med) / mad
    return list(np.abs(mz) > 3.5)


def aggregate(records, references_norm, control_group_norm, control_agg="median", drop_outliers=False):
    warns = []
    use = [r for r in records if r.stype == "sample" and r.group and math.isfinite(r.N0_cal)]
    buck = defaultdict(list)
    for r in use:
        buck[(norm(r.group), str(r.replicate), r.gene_norm)].append(r)
    vals = {}
    for (g, rep, gn), rs in buck.items():
        cqs = [x.cq for x in rs if x.cq is not None and math.isfinite(x.cq)]
        n0c = [x.N0_cal for x in rs]
        cv = float(np.std(n0c, ddof=1) / np.mean(n0c)) if len(n0c) > 1 and np.mean(n0c) > 0 else float("nan")
        vals[(g, rep, gn)] = ExprValue(
            group=rs[0].group, replicate=rep, gene=rs[0].gene, gene_norm=gn, role=rs[0].role,
            n_tech=len(rs), mean_cq=float(np.mean(cqs)) if cqs else float("nan"),
            mean_N0=float(np.mean([x.N0 for x in rs])), mean_N0_cal=float(np.mean(n0c)), tech_cv=cv)

    # QC 1: technical-well count anomaly (compare to the gene's usual count = mode)
    per_gene_counts = defaultdict(list)
    for (g, rep, gn), v in vals.items():
        per_gene_counts[gn].append(v.n_tech)
    usual = {gn: max(set(c), key=c.count) for gn, c in per_gene_counts.items()}
    for v in vals.values():
        if v.n_tech != usual.get(v.gene_norm):
            v.n_tech_flag = "check"
            warns.append(f"{v.gene} / {v.group} rep{v.replicate}: {v.n_tech} technical wells "
                         f"(usual {usual.get(v.gene_norm)}) - possible sample-name collision.")

    # QC 2: high technical CV
    for v in vals.values():
        if math.isfinite(v.tech_cv) and v.tech_cv > 0.5:
            warns.append(f"{v.gene} / {v.group} rep{v.replicate}: high technical CV = {v.tech_cv:.0%}.")

    # NRQ = calibrated N0 / geomean(reference N0), per (group, replicate)
    bios = sorted({(g, rep) for (g, rep, gn) in vals})
    if references_norm:
        ref_geo = {}
        for (g, rep) in bios:
            rv = [vals[(g, rep, rn)].mean_N0_cal for rn in references_norm
                  if (g, rep, rn) in vals and vals[(g, rep, rn)].mean_N0_cal > 0]
            if rv:
                ref_geo[(g, rep)] = _geomean(rv)
        for (g, rep, gn), v in vals.items():
            rg = ref_geo.get((g, rep))
            if rg and rg > 0:
                v.nrq = v.mean_N0_cal / rg
    else:
        warns.append("no reference gene: NRQ = calibrated N0 (not normalised).")
        for v in vals.values():
            v.nrq = v.mean_N0_cal

    # QC 3: biological-replicate outliers (MAD on log2 NRQ), per (gene, group)
    by_gg = defaultdict(list)
    for v in vals.values():
        if math.isfinite(v.nrq) and v.nrq > 0:
            by_gg[(v.gene_norm, norm(v.group))].append(v)
    for grp_vals in by_gg.values():
        logs = [math.log2(v.nrq) for v in grp_vals]
        for v, is_out in zip(grp_vals, _mad_outliers(logs)):
            if is_out:
                v.outlier = "outlier"
                warns.append(f"{v.gene} / {v.group} rep{v.replicate}: biological OUTLIER "
                             f"(NRQ={v.nrq:.3g}); check before trusting its fold change.")

    # fold change vs control, using a robust central value of the control group
    def central(xs):
        xs = [x for x in xs if math.isfinite(x)]
        if not xs:
            return None
        if control_agg == "mean":
            return float(np.mean(xs))
        if control_agg == "geomean":
            pos = [x for x in xs if x > 0]
            return float(np.exp(np.mean(np.log(pos)))) if pos else None
        return float(np.median(xs))            # default: median (robust to one outlier)

    if control_group_norm:
        ctrl_ref = {}
        tmp = defaultdict(list)
        for v in vals.values():
            if norm(v.group) == control_group_norm and math.isfinite(v.nrq):
                if drop_outliers and v.outlier == "outlier":
                    continue
                tmp[v.gene_norm].append(v.nrq)
        for gn, xs in tmp.items():
            c = central(xs)
            if c and c > 0:
                ctrl_ref[gn] = c
        for v in vals.values():
            cr = ctrl_ref.get(v.gene_norm)
            if cr and cr > 0 and math.isfinite(v.nrq):
                v.fold = v.nrq / cr
                v.log2_fold = float(np.log2(v.fold)) if v.fold > 0 else float("nan")
    else:
        warns.append("no control group (Control=yes): fold change not computed.")
    return list(vals.values()), warns


def reference_stability(records, references_norm, gene_display):
    """geNorm-style stability M per reference gene (mean pairwise SD of log2 ratios
    across biological samples). Lower M = more stable. M<0.5 good, <1.0 acceptable."""
    if len(references_norm) < 2:
        return []
    # calibrated N0 per (group,replicate,gene) averaged over technical wells
    buck = defaultdict(list)
    for r in records:
        if r.stype == "sample" and r.group and r.gene_norm in references_norm and math.isfinite(r.N0_cal) and r.N0_cal > 0:
            buck[(norm(r.group), str(r.replicate), r.gene_norm)].append(r.N0_cal)
    n0 = {k: float(np.mean(v)) for k, v in buck.items()}
    samples = sorted({(g, rep) for (g, rep, gn) in n0})
    refs = sorted(references_norm)
    rows = []
    for j in refs:
        As = []
        for k in refs:
            if k == j:
                continue
            diffs = [math.log2(n0[(g, rep, j)] / n0[(g, rep, k)])
                     for (g, rep) in samples if (g, rep, j) in n0 and (g, rep, k) in n0]
            if len(diffs) >= 2:
                As.append(float(np.std(diffs, ddof=1)))
        M = float(np.mean(As)) if As else float("nan")
        verdict = "good" if M < 0.5 else "acceptable" if M < 1.0 else "UNSTABLE"
        rows.append({"reference_gene": gene_display.get(j, j), "geNorm_M": M, "stability": verdict})
    return rows


def bootstrap_fold_ci(values, control_group_norm, control_agg, n_boot=2000, seed=0, drop_outliers=False):
    """Percentile bootstrap 95% CI on fold change per (gene, group)."""
    if not control_group_norm or n_boot <= 0:
        return {}
    if drop_outliers:
        values = [v for v in values if v.outlier != "outlier"]
    rng = np.random.default_rng(seed)
    by = defaultdict(lambda: defaultdict(list))   # gene_norm -> group_norm -> [nrq]
    disp = {}
    for v in values:
        if v.role == "target" and math.isfinite(v.nrq) and v.nrq > 0:
            by[v.gene_norm][norm(v.group)].append(v.nrq)
            disp[(v.gene_norm, norm(v.group))] = (v.gene, v.group)

    def central(xs):
        if control_agg == "mean":
            return np.mean(xs)
        if control_agg == "geomean":
            return np.exp(np.mean(np.log(xs)))
        return np.median(xs)

    out = {}
    for gn, groups in by.items():
        ctrl = groups.get(control_group_norm)
        if not ctrl or len(ctrl) < 2:
            continue
        ctrl = np.array(ctrl)
        for gg, arr in groups.items():
            if gg == control_group_norm or len(arr) < 2:
                continue
            arr = np.array(arr)
            folds = np.empty(n_boot)
            for b in range(n_boot):
                cs = central(rng.choice(ctrl, len(ctrl), replace=True))
                ts = np.mean(rng.choice(arr, len(arr), replace=True))
                folds[b] = ts / cs if cs > 0 else np.nan
            folds = folds[np.isfinite(folds)]
            if len(folds):
                lo, hi = np.percentile(folds, [2.5, 97.5])
                out[disp[(gn, gg)]] = (float(lo), float(hi))
    return out


# --------------------------------------------------------------------------- #
# statistics: choice of many-to-one test vs control + omnibus
#
#   qPCR relative quantities are ratio data (log-normal), so PARAMETRIC tests
#   are run on log2(NRQ). Because you compare each treatment to ONE control:
#     * dunnett     (default) Dunnett's test, parametric many-to-one, FWER-controlled
#     * ttest       Welch t-test each-vs-control, Benjamini-Hochberg adjusted
#     * mannwhitney Mann-Whitney U each-vs-control, BH adjusted (non-parametric)
#     * dunn        Dunn's test vs control after Kruskal-Wallis, BH (non-parametric)
#   Omnibus per gene is always reported: one-way ANOVA (on log2 NRQ) and
#   Kruskal-Wallis (non-parametric).
# --------------------------------------------------------------------------- #
@dataclass
class StatRow:
    gene: str
    group: str
    is_control: bool
    n: int
    mean_nrq: float
    sd_nrq: float
    sem_nrq: float
    mean_fold: float
    log2_fold: float
    test: str = ""
    statistic: float | None = None
    p_value: float | None = None
    p_adj: float | None = None
    stars: str = ""


def _bh(pvals):
    idx = [i for i, p in enumerate(pvals) if p is not None and math.isfinite(p)]
    m = len(idx)
    out = [None] * len(pvals)
    if not m:
        return out
    order = sorted(idx, key=lambda i: pvals[i])
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        prev = min(prev, pvals[i] * m / (m - rank + 1))
        out[i] = min(prev, 1.0)
    return out


def _stars(p):
    if p is None or not math.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def _dunn_vs_control(groups, control_key):
    """Dunn's test comparing each group to the control (rank-based, tie-corrected)."""
    from scipy.stats import rankdata, norm as _norm
    names = list(groups)
    allv = np.concatenate([groups[k] for k in names])
    ranks = rankdata(allv)
    N = len(allv)
    # mean rank per group
    mr, sizes, pos = {}, {}, 0
    for k in names:
        n = len(groups[k])
        mr[k] = float(np.mean(ranks[pos:pos + n]))
        sizes[k] = n
        pos += n
    # tie correction
    _, counts = np.unique(allv, return_counts=True)
    ties = np.sum(counts ** 3 - counts)
    sigma2 = (N * (N + 1) / 12.0) - ties / (12.0 * (N - 1)) if N > 1 else float("nan")
    out = {}
    for k in names:
        if k == control_key:
            continue
        se = math.sqrt(sigma2 * (1.0 / sizes[k] + 1.0 / sizes[control_key])) if sigma2 > 0 else float("nan")
        if not se or not math.isfinite(se):
            out[k] = (float("nan"), None)
            continue
        z = (mr[k] - mr[control_key]) / se
        p = float(2 * _norm.sf(abs(z)))
        out[k] = (float(z), p)
    return out


def compute_stats(values, control_group_norm, test="dunnett", drop_outliers=False):
    if drop_outliers:
        values = [v for v in values if v.outlier != "outlier"]
    try:
        from scipy import stats as ss
        have = True
    except ModuleNotFoundError:
        have = False
    rows, omni = [], []
    targets = sorted({v.gene for v in values if v.role == "target"}, key=norm_gene)
    for gene in targets:
        vv = [v for v in values if v.gene == gene]
        groups = sorted({v.group for v in vv}, key=norm)
        nrq_by = {g: np.array([v.nrq for v in vv if v.group == g and math.isfinite(v.nrq)]) for g in groups}
        log_by = {g: np.log2(nrq_by[g][nrq_by[g] > 0]) for g in groups}
        ctrl_name = next((g for g in groups if control_group_norm and norm(g) == control_group_norm), None)

        # ---- omnibus (always) ----
        F = Pf = H = Ph = None
        arrs_log = [log_by[g] for g in groups if len(log_by[g]) >= 2]
        arrs_raw = [nrq_by[g] for g in groups if len(nrq_by[g]) >= 2]
        if have and len(arrs_log) >= 2:
            try:
                F, Pf = map(float, ss.f_oneway(*arrs_log))
            except Exception:
                F = Pf = None
        if have and len(arrs_raw) >= 2:
            try:
                H, Ph = map(float, ss.kruskal(*arrs_raw))
            except Exception:
                H = Ph = None
        omni.append({"gene": gene, "n_groups": len(arrs_log),
                     "ANOVA_F": F, "ANOVA_p": Pf, "ANOVA_sig": _stars(Pf),
                     "KruskalWallis_H": H, "KruskalWallis_p": Ph, "KruskalWallis_sig": _stars(Ph)})

        # ---- post-hoc vs control ----
        pmap = {}
        smap = {}
        if have and ctrl_name is not None:
            treats = [g for g in groups if g != ctrl_name]
            if test == "dunnett":
                good = [g for g in treats if len(log_by[g]) >= 2]
                if len(log_by[ctrl_name]) >= 2 and good:
                    try:
                        res = ss.dunnett(*[log_by[g] for g in good], control=log_by[ctrl_name],
                                         random_state=np.random.default_rng(20240101))
                        for g, st, p in zip(good, np.atleast_1d(res.statistic), np.atleast_1d(res.pvalue)):
                            pmap[g] = float(p); smap[g] = float(st)
                    except Exception:
                        pass
            elif test == "dunn":
                gd = {g: nrq_by[g] for g in groups if len(nrq_by[g]) >= 1}
                if ctrl_name in gd and len(gd) >= 2:
                    for g, (z, p) in _dunn_vs_control(gd, ctrl_name).items():
                        raw = p
                        pmap[g] = raw; smap[g] = z
                    for g, adj in zip(list(pmap), _bh([pmap[g] for g in pmap])):
                        pmap[g] = adj                       # BH-adjust Dunn
            else:
                for g in treats:
                    a, c = (log_by[g], log_by[ctrl_name]) if test == "ttest" else (nrq_by[g], nrq_by[ctrl_name])
                    if len(a) >= 2 and len(c) >= 2:
                        try:
                            if test == "ttest":
                                st, p = ss.ttest_ind(a, c, equal_var=False)
                            else:
                                st, p = ss.mannwhitneyu(a, c, alternative="two-sided")
                            pmap[g] = float(p); smap[g] = float(st)
                        except Exception:
                            pass
                for g, adj in zip(list(pmap), _bh([pmap[g] for g in pmap])):
                    pmap[g] = adj                            # BH for t-test / MWU

        for g in groups:
            nrq = nrq_by[g]
            fold = np.array([v.fold for v in vv if v.group == g and math.isfinite(v.fold)])
            is_ctrl = bool(ctrl_name and g == ctrl_name)
            n = len(nrq)
            padj = pmap.get(g)
            rows.append(StatRow(
                gene=gene, group=g, is_control=is_ctrl, n=n,
                mean_nrq=float(np.mean(nrq)) if n else float("nan"),
                sd_nrq=float(np.std(nrq, ddof=1)) if n > 1 else float("nan"),
                sem_nrq=float(np.std(nrq, ddof=1) / math.sqrt(n)) if n > 1 else float("nan"),
                mean_fold=float(np.mean(fold)) if len(fold) else float("nan"),
                log2_fold=float(np.log2(np.mean(fold))) if len(fold) and np.mean(fold) > 0 else float("nan"),
                test=("-" if is_ctrl else test), statistic=smap.get(g),
                p_value=padj, p_adj=padj, stars=("" if is_ctrl else _stars(padj))))
    return rows, omni, have


# --------------------------------------------------------------------------- #
# plots
# --------------------------------------------------------------------------- #
_PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
            "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"]


def plot_expression(values, stat_rows, control_group_norm, out_path, metric="fold"):
    targets = sorted({v.gene for v in values if v.role == "target"}, key=norm_gene)
    if not targets:
        return False
    groups = sorted({v.group for v in values if v.role == "target"}, key=norm)
    if control_group_norm:
        groups = [g for g in groups if norm(g) == control_group_norm] + \
                 [g for g in groups if norm(g) != control_group_norm]
    look = {(s.gene, s.group): s for s in stat_rows}
    ncol = min(len(targets), 2)
    nrow = math.ceil(len(targets) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 4.2 * nrow), squeeze=False)
    ykey = "fold" if metric == "fold" else "nrq"
    ylabel = "Fold change vs control" if metric == "fold" else "NRQ"
    for ti, gene in enumerate(targets):
        ax = axes[ti // ncol][ti % ncol]
        for ci, grp in enumerate(groups):
            pts = [getattr(v, ykey) for v in values
                   if v.gene == gene and v.group == grp and math.isfinite(getattr(v, ykey))]
            if not pts:
                continue
            m = float(np.mean(pts))
            err = float(np.std(pts, ddof=1) / math.sqrt(len(pts))) if len(pts) > 1 else 0.0
            ax.bar(ci, m, width=0.62, color=_PALETTE[ci % len(_PALETTE)], alpha=0.85,
                   edgecolor="black", linewidth=0.6, zorder=2)
            if err > 0:
                ax.errorbar(ci, m, yerr=err, fmt="none", ecolor="black", elinewidth=1, capsize=4, zorder=3)
            jit = np.random.default_rng(ci).uniform(-0.12, 0.12, size=len(pts))
            ax.scatter(np.full(len(pts), ci) + jit, pts, s=26, color="black", alpha=0.7, zorder=4)
            st = look.get((gene, grp))
            if st and st.stars and st.stars not in ("", "ns"):
                ax.text(ci, m + err + 0.03 * max(1.0, m), st.stars, ha="center", va="bottom",
                        fontsize=12, fontweight="bold")
        if metric == "fold":
            ax.axhline(1.0, color="grey", lw=0.8, ls="--", zorder=1)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=35, ha="right")
        ax.set_title(gene, fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.margins(y=0.18)
    for k in range(len(targets), nrow * ncol):
        axes[k // ncol][k % ncol].set_visible(False)
    fig.suptitle("Relative gene expression", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, bbox_inches="tight")
    fig.savefig(out_path.rsplit(".", 1)[0] + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return True


def plot_efficiency(summaries, records, out_path):
    keys = [(s["plate"], s["gene"]) for s in summaries]
    if not keys:
        return False
    labels = [f"{g}\n{p}" for (p, g) in keys]
    fig, ax = plt.subplots(figsize=(1.5 * len(keys) + 2, 4.6))
    for i, s in enumerate(summaries):
        es = [r.com_E for r in records if r.plate == s["plate"] and r.gene == s["gene"] and math.isfinite(r.com_E)]
        used = [r.com_E for r in records if r.plate == s["plate"] and r.gene == s["gene"]
                and r.used_for_mean and math.isfinite(r.com_E)]
        c = _PALETTE[i % len(_PALETTE)]
        if es:
            ax.scatter(np.full(len(es), i) + np.random.default_rng(i).uniform(-0.14, 0.14, len(es)),
                       es, s=22, color=c, alpha=0.3, linewidths=0)
        if used:
            ax.scatter(np.full(len(used), i) + np.random.default_rng(i + 9).uniform(-0.14, 0.14, len(used)),
                       used, s=26, color=c, alpha=0.9, edgecolors="black", linewidths=0.4)
        me = s["mean_E"]
        if isinstance(me, float) and math.isfinite(me):
            ax.hlines(me, i - 0.32, i + 0.32, color="black", lw=2.2)
            ax.text(i, me, f" {me:.3f}", va="center", ha="left", fontsize=9, fontweight="bold")
    ax.axhline(2.0, color="grey", ls=":", lw=0.8)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylabel("Amplification factor E")
    ax.set_title("PCR efficiency per gene per plate", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def _fmt(x, nd=6):
    if x is None:
        return ""
    if isinstance(x, bool):
        return "yes" if x else ""
    if isinstance(x, float):
        if not math.isfinite(x):
            return ""
        if x != 0 and abs(x) < 1e-4:
            return f"{x:.4e}"
        return round(x, nd)
    return x


def write_csv(path, rows, fields):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(r.get(k)) for k in fields})


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def run(manifest_path, sample_map_path, out_prefix, *, test="dunnett", make_plots=True,
        window=4, window_max=6, eff_outlier=0.10, efficiency_mode="amplicon_mean",
        fixed_efficiency=2.0, control_agg="mean", n_boot=2000, drop_outliers=False,
        exclude_wells_path=None, force_baseline_subtracted=False):
    plates = load_manifest(manifest_path)
    samples, control_groups, gene_role, references, gene_display, gene_calibrate = load_metadata(sample_map_path)
    exclude_wells = load_exclude_wells(exclude_wells_path)
    if exclude_wells:
        print(f"Excluding {len(exclude_wells)} well(s) listed in {os.path.basename(exclude_wells_path)}.")
    control_group_norm = next(iter(control_groups)) if control_groups else None
    if len(control_groups) > 1:
        print(f"  NOTE: multiple control groups {sorted(control_groups)}; using all.")
    emode = f"{efficiency_mode}" + (f" (E={fixed_efficiency})" if efficiency_mode == "fixed" else "")
    print(f"Efficiency mode for N0: {emode}")

    all_records, all_summaries = [], []
    print("Plates:")
    for pf in plates:
        recs, summ, kind = process_plate(pf, samples, gene_role, window=window, window_max=window_max,
                                         eff_outlier=eff_outlier, efficiency_mode=efficiency_mode,
                                         fixed_efficiency=fixed_efficiency, exclude_wells=exclude_wells,
                                         force_baseline_subtracted=force_baseline_subtracted)
        all_records += recs
        all_summaries += summ
        print(f"  [{pf.plate}] kind={kind} genes={sorted({s['gene'] for s in summ})}")
        for s in summ:
            if isinstance(s["mean_E"], float) and math.isfinite(s["mean_E"]):
                print(f"       {s['gene']:>10} ({s['role']}): E={s['mean_E']:.3f} "
                      f"({s['efficiency_percent']:.1f}%) n={s['n_wells_used']}/{s['n_wells_total']}")

    calibrate_genes = {gn for gn, yes in gene_calibrate.items() if yes}
    ref_names = sorted(gene_display[gn] for gn in references)
    print(f"  references (geomean): {ref_names or '— none —'}")
    print(f"  calibrated genes: {sorted(gene_display[gn] for gn in calibrate_genes) or '— none —'}")
    cal_factors, cal_warns = apply_calibration(all_records, calibrate_genes)
    for w in cal_warns:
        print("  NOTE:", w)
    values, warns = aggregate(all_records, sorted(references), control_group_norm,
                              control_agg=control_agg, drop_outliers=drop_outliers)
    for w in warns:
        print("  NOTE:", w)
    print(f"  fold uses control center = {control_agg}"
          + ("; flagged biological outliers DROPPED from stats/fold." if drop_outliers
             else "; outliers flagged but kept (use --drop-outliers to exclude)."))
    ref_stab = reference_stability(all_records, sorted(references), gene_display)
    for r in ref_stab:
        if r["stability"] != "good":
            print(f"  NOTE: reference {r['reference_gene']} geNorm M={r['geNorm_M']:.2f} -> {r['stability']}.")
    ci = bootstrap_fold_ci(values, control_group_norm, control_agg, n_boot=n_boot, drop_outliers=drop_outliers)
    stat_rows, omni, have = compute_stats(values, control_group_norm, test=test, drop_outliers=drop_outliers)
    if not have:
        print("  NOTE: scipy missing -> stats skipped. pip install scipy.")

    written = []

    eff_rows = [{
        "source_file": r.source_file, "plate": r.plate, "gene": r.gene, "role": r.role, "well": r.well,
        "sample": r.sample, "group": r.group, "replicate": r.replicate, "sample_type": r.stype,
        "Cq": r.cq, "Cq_source": r.cq_source, "window_start_cycle": r.win_start,
        "window_end_cycle": r.win_end, "window_points": r.win_n, "individual_E": r.ind_E,
        "common_window_E": r.com_E, "efficiency_percent": (r.com_E - 1) * 100 if math.isfinite(r.com_E) else None,
        "R2": r.com_r2, "gene_mean_E": r.gene_mean_E, "E_used_for_N0": r.E_used,
        "used_for_gene_mean": r.used_for_mean,
        "N0_relative": r.N0, "quality": r.quality, "flags": r.flags,
    } for r in all_records]
    eff_fields = ["source_file", "plate", "gene", "role", "well", "sample", "group", "replicate",
                  "sample_type", "Cq", "Cq_source", "window_start_cycle", "window_end_cycle",
                  "window_points", "individual_E", "common_window_E", "efficiency_percent", "R2",
                  "gene_mean_E", "E_used_for_N0", "used_for_gene_mean", "N0_relative", "quality", "flags"]
    p = f"{out_prefix}_efficiencies.csv"; write_csv(p, eff_rows, eff_fields); written.append(p)

    expr_rows = [{
        "group": v.group, "replicate": v.replicate, "gene": v.gene, "role": v.role, "n_technical": v.n_tech,
        "technical_CV": v.tech_cv, "n_tech_check": v.n_tech_flag, "biological_outlier": v.outlier,
        "mean_Cq": v.mean_cq, "mean_N0": v.mean_N0, "mean_N0_calibrated": v.mean_N0_cal, "NRQ": v.nrq,
        "fold_vs_control": v.fold, "log2_fold": v.log2_fold,
    } for v in sorted(values, key=lambda v: (norm_gene(v.gene), norm(v.group), v.replicate))]
    expr_fields = ["group", "replicate", "gene", "role", "n_technical", "technical_CV", "n_tech_check",
                   "biological_outlier", "mean_Cq", "mean_N0", "mean_N0_calibrated", "NRQ",
                   "fold_vs_control", "log2_fold"]
    p = f"{out_prefix}_expression.csv"; write_csv(p, expr_rows, expr_fields); written.append(p)

    st_rows = [{
        "gene": s.gene, "group": s.group, "is_control": s.is_control, "n": s.n, "mean_NRQ": s.mean_nrq,
        "sd_NRQ": s.sd_nrq, "sem_NRQ": s.sem_nrq, "mean_fold": s.mean_fold, "log2_fold": s.log2_fold,
        "fold_CI95_low": (ci.get((s.gene, s.group)) or (None, None))[0],
        "fold_CI95_high": (ci.get((s.gene, s.group)) or (None, None))[1],
        "test": s.test, "statistic": s.statistic, "p_value": s.p_value, "p_adjusted": s.p_adj,
        "significance": s.stars,
    } for s in stat_rows]
    st_fields = ["gene", "group", "is_control", "n", "mean_NRQ", "sd_NRQ", "sem_NRQ", "mean_fold",
                 "log2_fold", "fold_CI95_low", "fold_CI95_high", "test", "statistic", "p_value",
                 "p_adjusted", "significance"]
    p = f"{out_prefix}_statistics.csv"; write_csv(p, st_rows, st_fields); written.append(p)

    if ref_stab:
        p = f"{out_prefix}_reference_stability.csv"
        write_csv(p, ref_stab, ["reference_gene", "geNorm_M", "stability"]); written.append(p)

    p = f"{out_prefix}_omnibus.csv"
    write_csv(p, omni, ["gene", "n_groups", "ANOVA_F", "ANOVA_p", "ANOVA_sig",
                        "KruskalWallis_H", "KruskalWallis_p", "KruskalWallis_sig"]); written.append(p)

    long_rows = [{
        "plate": r.plate, "source_file": r.source_file, "well": r.well, "gene": r.gene, "role": r.role,
        "sample": r.sample, "group": r.group, "replicate": r.replicate, "sample_type": r.stype,
        "is_control": r.is_control, "Cq": r.cq, "gene_mean_E": r.gene_mean_E, "N0_relative": r.N0,
        "calibration_factor": r.cf, "N0_calibrated": r.N0_cal, "quality": r.quality,
    } for r in all_records]
    long_fields = ["plate", "source_file", "well", "gene", "role", "sample", "group", "replicate",
                   "sample_type", "is_control", "Cq", "gene_mean_E", "N0_relative", "calibration_factor",
                   "N0_calibrated", "quality"]
    p = f"{out_prefix}_long.csv"; write_csv(p, long_rows, long_fields); written.append(p)

    if cal_factors:
        p = f"{out_prefix}_calibration_factors.csv"
        write_csv(p, cal_factors, ["gene", "plate", "calibrator_N0", "cross_plate_geomean",
                                   "calibration_factor", "n_calibrator_wells"]); written.append(p)

    if make_plots:
        metric = "fold" if any(math.isfinite(v.fold) for v in values) else "nrq"
        if any(math.isfinite(getattr(v, "fold" if metric == "fold" else "nrq")) for v in values):
            if plot_expression(values, stat_rows, control_group_norm, f"{out_prefix}_expression.png", metric=metric):
                written += [f"{out_prefix}_expression.png", f"{out_prefix}_expression.pdf"]
        if plot_efficiency(all_summaries, all_records, f"{out_prefix}_efficiency.png"):
            written.append(f"{out_prefix}_efficiency.png")

    print("\nWritten files:")
    for w in written:
        print("  ", w)
    return written


def main() -> int:
    ap = argparse.ArgumentParser(description="qPCR efficiency + relative expression (manifest + one sample_map).")
    ap.add_argument("--manifest", required=True, help="Plate, Amp, Cq")
    ap.add_argument("--sample-map", required=True,
                    help="ONE table: SAMPLE rows (Plate,Sample,Group,Replicate,Type,Control) + GENE rows (Gene,Role)")
    ap.add_argument("--out", default="qpcr_out/RESULT")
    ap.add_argument("--test", choices=["dunnett", "ttest", "mannwhitney", "dunn"], default="dunnett",
                    help="Post-hoc each-group-vs-control test (default dunnett).")
    ap.add_argument("--efficiency-mode", choices=["amplicon_mean", "per_well", "fixed"],
                    default="amplicon_mean",
                    help="Which E feeds N0=E^(-Cq): amplicon_mean (trimmed mean per gene, default) | "
                         "per_well (each well's own E) | fixed (use --fixed-efficiency).")
    ap.add_argument("--fixed-efficiency", type=float, default=2.0,
                    help="E used when --efficiency-mode fixed (2.0 = perfect doubling).")
    ap.add_argument("--window", type=int, default=4,
                    help="Min points in the log-linear window fit.")
    ap.add_argument("--window-max", type=int, default=6, help="Max points in the window fit.")
    ap.add_argument("--eff-outlier", type=float, default=0.10,
                    help="Half-width around median E for outlier well exclusion.")
    ap.add_argument("--fold-center", choices=["mean", "median", "geomean"], default="mean",
                    help="Central value of the control group for fold change (default median, robust to outliers).")
    ap.add_argument("--drop-outliers", action="store_true",
                    help="Exclude MAD-flagged biological-replicate outliers from fold and statistics.")
    ap.add_argument("--n-boot", type=int, default=2000,
                    help="Bootstrap resamples for fold 95%% CI (0 to disable).")
    ap.add_argument("--exclude-wells", default=None,
                    help="Optional CSV (Plate, Well) of individual wells to drop before fitting.")
    ap.add_argument("--baseline-subtracted", action="store_true",
                    help="Force treating amplification files as already baseline-subtracted.")
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args()
    run(a.manifest, a.sample_map, a.out, test=a.test, make_plots=not a.no_plots,
        window=a.window, window_max=a.window_max, eff_outlier=a.eff_outlier,
        efficiency_mode=a.efficiency_mode, fixed_efficiency=a.fixed_efficiency,
        control_agg=a.fold_center, n_boot=a.n_boot, drop_outliers=a.drop_outliers,
        exclude_wells_path=a.exclude_wells, force_baseline_subtracted=a.baseline_subtracted)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
