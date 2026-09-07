#!/usr/bin/env python3
"""
tes_screening.py - Monte Carlo screening of high-temperature thermal energy
storage media with uncertainty propagation and a rank-stability test.

What it does
------------
1. Reads a long-format property table (one row per material+property+source).
   It prefers a curated table named tes_property_table.csv in the same
   directory; if that is missing it falls back to tes_screening_demo_inputs.csv
   and marks the run "demo_placeholder" in the output. Override with --input.
2. Builds three figures of merit (FOMs) per material:
     energy_density_kWh_per_m3  volumetric stored energy over a temperature
                                swing, sensible plus latent if the melting
                                point falls inside the swing. Higher is better.
     cost_per_kWh_th_USD        media cost only, indicative. Lower is better.
     diffusivity_mm2_per_s      thermal diffusivity k/(rho*cp), used only as a
                                crude charge/discharge-rate proxy for
                                conduction-limited media. Higher is better.
3. Propagates uncertainty by Monte Carlo: each draw picks one competing source
   row per material+property uniformly at random (this is how source
   disagreement enters), then samples uniformly inside that row's
   [value_low, value_high]; if the interval is missing it uses
   value * (1 +/- DEFAULT_REL_UNCERTAINTY).
4. Reports the distribution of ranks per material per FOM, pairwise
   "A beats B" probabilities, and flags pairs that are indistinguishable
   given the input intervals.
5. Separately enumerates the discrete source-disagreement scenarios (every
   combination of competing rows, nominal values, no sampling) and reports
   how many distinct orderings source choice alone can produce.
6. Writes everything to tes_screening_output.json.

Temperature windows
-------------------
Property tables in this field usually give a maximum operating temperature and
a melting point but no cold-end temperature, because the cold end is a system
design choice. The script therefore derives a window and records the basis:
  upper end  max_operating_temperature if present; else, for a solid sensible
             medium, melting_point - SOLID_HEADROOM_C; else the material is
             excluded.
  lower end  min_operating_temperature if present; else for a latent/PCM
             family melting_point +/- PCM_HALF_SWING_C (capped above by the
             upper end); else for a liquid medium melting_point +
             FREEZE_MARGIN_C; else the assumed DEFAULT_COLD_END_C.
These fallbacks are assumptions, not data. Every material's basis string is in
the output, and materials relying on an assumed cold end are listed.

Two swing scenarios are computed:
  own_window     each material over its own derived window. Materials are then
                 compared across different swings, which favours wide-window
                 media.
  common_window  a fixed swing (default 500-700 C). Materials whose derived
                 window does not cover it are excluded and listed with the
                 reason. This is the like-for-like comparison.

Determinism
-----------
Fixed seed (default 20260907), sorted iteration everywhere, no timestamps in
the output. Two runs with the same input file and arguments produce
byte-identical JSON.

Usage
-----
  python tes_screening.py
  python tes_screening.py --input tes_screening_demo_inputs.csv
  python tes_screening.py --draws 2000 --seed 20260907 --cold 500 --hot 700
"""

import argparse
import csv
import hashlib
import itertools
import json
import os
import platform
import sys
from collections import defaultdict

import numpy as np

try:
    import scipy
    SCIPY_VERSION = scipy.__version__
except Exception:  # pragma: no cover
    SCIPY_VERSION = "not available"

# ----------------------------------------------------------------------------
# Configuration constants (all echoed into the output for provenance)
# ----------------------------------------------------------------------------
DEFAULT_SEED = 20260907
DEFAULT_DRAWS = 2000
DEFAULT_COLD_C = 500.0
DEFAULT_HOT_C = 700.0
DEFAULT_REL_UNCERTAINTY = 0.10   # used when a row has no low/high interval
FREEZE_MARGIN_C = 20.0           # cold-end margin above melting point (liquid)
PCM_HALF_SWING_C = 50.0          # half-width of the swing around a PCM melt
SOLID_HEADROOM_C = 100.0         # headroom below melting point for a solid
DEFAULT_COLD_END_C = 200.0       # assumed cold end when nothing else applies
INDISTINGUISHABLE_BAND = (0.40, 0.60)   # P(A better than B) inside -> a tie
MAX_SCENARIO_COMBINATIONS = 256  # cap on enumerated source-choice scenarios

# Canonical property names and accepted synonyms (lower case, underscored).
PROPERTY_SYNONYMS = {
    "density": "density",
    "bulk_density": "density",
    "mass_density": "density",
    "heat_capacity": "heat_capacity",
    "specific_heat": "heat_capacity",
    "specific_heat_capacity": "heat_capacity",
    "cp": "heat_capacity",
    "thermal_conductivity": "thermal_conductivity",
    "conductivity": "thermal_conductivity",
    "k": "thermal_conductivity",
    "melting_point": "melting_point",
    "melting_temperature": "melting_point",
    "melt_temperature": "melting_point",
    "freezing_point": "melting_point",
    "liquidus": "melting_point",
    "latent_heat": "latent_heat",
    "latent_heat_fusion": "latent_heat",
    "latent_heat_of_fusion": "latent_heat",
    "enthalpy_of_fusion": "latent_heat",
    "heat_of_fusion": "latent_heat",
    "max_operating_temperature": "max_operating_temperature",
    "maximum_operating_temperature": "max_operating_temperature",
    "max_operating_temp": "max_operating_temperature",
    "max_temperature": "max_operating_temperature",
    "upper_operating_temperature": "max_operating_temperature",
    "min_operating_temperature": "min_operating_temperature",
    "minimum_operating_temperature": "min_operating_temperature",
    "min_operating_temp": "min_operating_temperature",
    "min_temperature": "min_operating_temperature",
    "lower_operating_temperature": "min_operating_temperature",
    "cost_per_tonne": "cost_per_tonne",
    "cost_indicative": "cost_per_tonne",
    "indicative_cost": "cost_per_tonne",
    "cost": "cost_per_tonne",
    "material_cost": "cost_per_tonne",
    "media_cost": "cost_per_tonne",
    "cost_per_ton": "cost_per_tonne",
}

TEMPERATURE_PROPERTIES = {
    "melting_point", "max_operating_temperature", "min_operating_temperature"}

# Readable unit keys; both these and the file's units go through norm_unit().
RAW_UNIT_FACTORS = {
    ("density", "kg/m3"): 1.0,
    ("density", "kg m-3"): 1.0,
    ("density", "kg/m^3"): 1.0,
    ("density", "g/cm3"): 1000.0,
    ("density", "g/cm^3"): 1000.0,
    ("density", "kg/L"): 1000.0,
    ("heat_capacity", "J/(kg K)"): 1.0,
    ("heat_capacity", "J/kg-K"): 1.0,
    ("heat_capacity", "J/kg/K"): 1.0,
    ("heat_capacity", "J kg-1 K-1"): 1.0,
    ("heat_capacity", "kJ/(kg K)"): 1000.0,
    ("heat_capacity", "kJ/kg-K"): 1000.0,
    ("heat_capacity", "kJ/kg/K"): 1000.0,
    ("heat_capacity", "J/(g K)"): 1000.0,
    ("thermal_conductivity", "W/(m K)"): 1.0,
    ("thermal_conductivity", "W/m-K"): 1.0,
    ("thermal_conductivity", "W/m/K"): 1.0,
    ("thermal_conductivity", "W m-1 K-1"): 1.0,
    ("latent_heat", "J/kg"): 1.0,
    ("latent_heat", "kJ/kg"): 1000.0,
    ("latent_heat", "J/g"): 1000.0,
    ("latent_heat", "MJ/kg"): 1.0e6,
    ("cost_per_tonne", "USD/tonne"): 1.0,
    ("cost_per_tonne", "USD_per_tonne"): 1.0,
    ("cost_per_tonne", "USD/t"): 1.0,
    ("cost_per_tonne", "$/tonne"): 1.0,
    ("cost_per_tonne", "$/t"): 1.0,
    ("cost_per_tonne", "USD/kg"): 1000.0,
    ("cost_per_tonne", "USD_per_kg"): 1000.0,
    ("cost_per_tonne", "$/kg"): 1000.0,
    ("cost_per_tonne", "EUR/tonne"): 1.0,   # treated as USD 1:1; logged
    ("cost_per_tonne", "EUR/t"): 1.0,
}

CANONICAL_UNITS = {
    "density": "kg/m3",
    "heat_capacity": "J/(kg K)",
    "thermal_conductivity": "W/(m K)",
    "latent_heat": "J/kg",
    "cost_per_tonne": "USD/tonne",
    "melting_point": "C",
    "max_operating_temperature": "C",
    "min_operating_temperature": "C",
}

REQUIRED_FOR_ENERGY = ("density", "heat_capacity")
FOM_DIRECTION = {                      # +1: higher is better, -1: lower better
    "energy_density_kWh_per_m3": +1,
    "cost_per_kWh_th_USD": -1,
    "diffusivity_mm2_per_s": +1,
}

PCM_FAMILY_MARKERS = ("pcm", "latent")
SOLID_FAMILY_MARKERS = ("solid", "concrete", "cementitious", "ceramic", "rock",
                        "sand", "silica", "graphite", "carbon", "refractory",
                        "basalt", "oxide", "brick")

WARNINGS = []


def warn(message):
    if message not in WARNINGS:
        WARNINGS.append(message)


# ----------------------------------------------------------------------------
# Parsing helpers
# ----------------------------------------------------------------------------
def norm_key(text):
    return (text or "").strip().lower().replace(" ", "_").replace("-", "_")


def norm_unit(text):
    text = (text or "").strip().lower()
    for junk in (" ", "(", ")", "*", "·", "⋅", "_", "-", "^"):
        text = text.replace(junk, "")
    return text


UNIT_FACTORS = {(prop, norm_unit(unit)): factor
                for (prop, unit), factor in RAW_UNIT_FACTORS.items()}


def to_float(text):
    if text is None:
        return None
    text = str(text).strip().replace(",", "")
    if text == "" or text.lower() in {"na", "n/a", "none", "null", "-", "?"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def convert(prop, unit, value):
    """Convert value to the canonical unit for prop. Returns (value, ok)."""
    if value is None:
        return None, True
    u = norm_unit(unit)
    if prop in TEMPERATURE_PROPERTIES:
        if u in {"c", "degc", "celsius", "oc", "°c", ""}:
            return value, True
        if u in {"k", "kelvin"}:
            return value - 273.15, True
        if u in {"f", "degf"}:
            return (value - 32.0) * 5.0 / 9.0, True
        return value, False
    if u == "":
        return value, False
    factor = UNIT_FACTORS.get((prop, u))
    if factor is None:
        return value, False
    if prop == "cost_per_tonne" and "eur" in u:
        warn("cost given in EUR was treated as USD 1:1")
    return value * factor, True


def read_table(path):
    """Return rows as dicts with canonical property names and canonical units."""
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit("input CSV has no header row: %s" % path)
        headers = {norm_key(h): h for h in reader.fieldnames if h}
        for column in ("material", "property", "value"):
            if column not in headers:
                raise SystemExit(
                    "input CSV is missing required column '%s'; found %s"
                    % (column, sorted(headers)))

        def get(record, name):
            key = headers.get(name)
            return record.get(key) if key else None

        for record in reader:
            material = (get(record, "material") or "").strip()
            raw_prop = norm_key(get(record, "property"))
            if not material or not raw_prop:
                continue
            prop = PROPERTY_SYNONYMS.get(raw_prop)
            if prop is None:
                warn("ignored unrecognised property '%s' (material: %s)"
                     % (raw_prop, material))
                continue
            unit = get(record, "unit")
            value = to_float(get(record, "value"))
            low = to_float(get(record, "value_low"))
            high = to_float(get(record, "value_high"))
            if value is None and low is not None and high is not None:
                value = 0.5 * (low + high)
            if value is None:
                warn("skipped row with no usable value: %s / %s"
                     % (material, prop))
                continue
            value, ok_v = convert(prop, unit, value)
            low, ok_l = convert(prop, unit, low)
            high, ok_h = convert(prop, unit, high)
            if not (ok_v and ok_l and ok_h):
                warn("unrecognised unit '%s' for %s / %s; value used unchanged"
                     % ((unit or "").strip(), material, prop))
            if low is not None and high is not None and low > high:
                low, high = high, low
            rows.append({
                "material": material,
                "family": (get(record, "family") or "").strip(),
                "property": prop,
                "value": value,
                "low": low,
                "high": high,
                "temperature_C": to_float(get(record, "temperature_c")),
                "source_citation": (get(record, "source_citation") or "").strip(),
                "source_locator": (get(record, "source_locator") or "").strip(),
                "disagreement_flag": (
                    get(record, "disagreement_flag") or "").strip(),
            })
    if not rows:
        raise SystemExit("no usable rows parsed from %s" % path)
    return rows


def interval(row):
    """Sampling interval for one row."""
    low, high = row["low"], row["high"]
    span = abs(row["value"]) * DEFAULT_REL_UNCERTAINTY
    if low is None:
        low = row["value"] - span
    if high is None:
        high = row["value"] + span
    if low > high:
        low, high = high, low
    return float(low), float(high)


# ----------------------------------------------------------------------------
# Temperature window policy (works elementwise on scalars or numpy arrays)
# ----------------------------------------------------------------------------
def derive_window(min_op, melt, max_op, family):
    """Return (t_lo, t_hi, basis, assumed_cold_end) or (None, None, reason, _)."""
    fam = (family or "").lower()
    is_pcm = any(marker in fam for marker in PCM_FAMILY_MARKERS)
    is_solid = any(marker in fam for marker in SOLID_FAMILY_MARKERS)
    basis = []

    if max_op is not None:
        t_hi = max_op
        basis.append("upper: max_operating_temperature from table")
    elif melt is not None and is_pcm:
        t_hi = melt + PCM_HALF_SWING_C
        basis.append("upper: melting_point + %.0f C (PCM stage, no "
                     "max_operating_temperature in table)" % PCM_HALF_SWING_C)
    elif melt is not None and is_solid:
        t_hi = melt - SOLID_HEADROOM_C
        basis.append("upper: melting_point - %.0f C headroom (solid medium, "
                     "no max_operating_temperature in table)" % SOLID_HEADROOM_C)
    else:
        return None, None, ("no upper temperature: neither "
                            "max_operating_temperature nor a usable "
                            "melting_point for a solid medium"), False

    assumed_cold = False
    if min_op is not None:
        t_lo = min_op
        basis.append("lower: min_operating_temperature from table")
    elif is_pcm and melt is not None:
        t_lo = melt - PCM_HALF_SWING_C
        t_hi = np.minimum(t_hi, melt + PCM_HALF_SWING_C)
        basis.append("lower: PCM stage, melting_point +/- %.0f C capped by the "
                     "upper limit" % PCM_HALF_SWING_C)
    elif melt is not None and not is_solid:
        t_lo = melt + FREEZE_MARGIN_C
        basis.append("lower: melting_point + %.0f C freeze margin (liquid "
                     "medium)" % FREEZE_MARGIN_C)
    else:
        t_lo = DEFAULT_COLD_END_C
        assumed_cold = True
        basis.append("lower: ASSUMED cold end %.0f C, no data in table"
                     % DEFAULT_COLD_END_C)
    return t_lo, t_hi, "; ".join(basis), assumed_cold


def build_materials(rows):
    """Group rows into {material: {...}} and derive nominal windows."""
    grouped = defaultdict(lambda: defaultdict(list))
    families = {}
    for row in rows:
        grouped[row["material"]][row["property"]].append(row)
        if row["family"]:
            families.setdefault(row["material"], row["family"])

    materials = {}
    for material in sorted(grouped):
        props = grouped[material]
        family = families.get(material, "")

        def nominal(prop):
            if prop not in props:
                return None
            return float(np.mean([r["value"] for r in props[prop]]))

        t_lo_nom, t_hi_nom, basis, assumed_cold = derive_window(
            nominal("min_operating_temperature"), nominal("melting_point"),
            nominal("max_operating_temperature"), family)
        midpoint = (0.5 * (t_lo_nom + t_hi_nom)
                    if (t_lo_nom is not None and t_hi_nom is not None) else None)

        # Where several rows exist for one property at different temperatures,
        # keep the temperature nearest the window midpoint; the rows remaining
        # at that temperature are treated as competing sources.
        selected = {}
        for prop in sorted(props):
            candidates = props[prop]
            temps = sorted({r["temperature_C"] for r in candidates
                            if r["temperature_C"] is not None})
            if len(temps) > 1 and midpoint is not None:
                best = min(temps, key=lambda t: abs(t - midpoint))
                keep = [r for r in candidates if r["temperature_C"] == best]
                warn("%s / %s: kept %d row(s) at %.0f C nearest the window "
                     "midpoint, dropped %d row(s) at other temperatures"
                     % (material, prop, len(keep), best,
                        len(candidates) - len(keep)))
                candidates = keep
            selected[prop] = candidates

        materials[material] = {
            "family": family,
            "properties": selected,
            "window_basis": basis,
            "assumed_cold_end": assumed_cold,
            "nominal_window_C": [t_lo_nom, t_hi_nom],
        }
    return materials


# ----------------------------------------------------------------------------
# Sampling and figures of merit
# ----------------------------------------------------------------------------
def sample_property(rng, rows, draws):
    """Pick a competing source row per draw, then sample inside its interval."""
    lows = np.array([interval(r)[0] for r in rows], dtype=float)
    highs = np.array([interval(r)[1] for r in rows], dtype=float)
    if len(rows) == 1:
        choice = np.zeros(draws, dtype=int)
    else:
        choice = rng.integers(0, len(rows), size=draws)
    unit_draws = rng.random(draws)
    return lows[choice] + unit_draws * (highs[choice] - lows[choice])


def compute_foms(density, cp, cond, latent, melt, t_lo, t_hi, cost_per_tonne):
    """Vectorised figures of merit in canonical units."""
    delta_t = np.maximum(t_hi - t_lo, 0.0)
    specific = cp * delta_t                       # J/kg sensible
    latent_fraction = None
    if latent is not None and melt is not None:
        inside = (melt >= t_lo) & (melt <= t_hi)
        specific = specific + np.where(inside, latent, 0.0)
        latent_fraction = float(np.mean(inside))
    out = {"energy_density_kWh_per_m3": density * specific / 3.6e6,
           "_delta_T_C": delta_t,
           "_latent_counted_fraction": latent_fraction}
    energy_mass = 1000.0 * specific / 3.6e6       # kWh_th per tonne
    if cost_per_tonne is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            out["cost_per_kWh_th_USD"] = np.where(
                energy_mass > 0, cost_per_tonne / np.where(
                    energy_mass > 0, energy_mass, 1.0), np.inf)
    if cond is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            out["diffusivity_mm2_per_s"] = cond / (density * cp) * 1.0e6
    return out


def rank_matrix(values_by_material, direction):
    """Rank 1 is best. Ties broken deterministically by material name order."""
    names = sorted(values_by_material)
    stack = np.vstack([values_by_material[n] for n in names])
    keys = stack if direction < 0 else -stack
    keys = np.where(np.isfinite(keys), keys, np.inf)
    order = np.argsort(keys, axis=0, kind="stable")
    ranks = np.empty_like(order)
    rows = np.broadcast_to(np.arange(len(names))[:, None], order.shape).copy()
    np.put_along_axis(ranks, order, rows, axis=0)
    return {n: ranks[i] + 1 for i, n in enumerate(names)}


def summarise_ranks(ranks_by_material, n_materials):
    summary = {}
    for material in sorted(ranks_by_material):
        r = ranks_by_material[material]
        counts = np.bincount(r, minlength=n_materials + 2)[1:n_materials + 1]
        summary[material] = {
            "mean_rank": round(float(np.mean(r)), 3),
            "median_rank": int(np.median(r)),
            "modal_rank": int(np.argmax(counts) + 1),
            "rank_p05": int(np.percentile(r, 5)),
            "rank_p95": int(np.percentile(r, 95)),
            "prob_rank_1": round(float(np.mean(r == 1)), 4),
            "prob_top_3": round(float(np.mean(r <= 3)), 4),
            "rank_histogram": [int(c) for c in counts],
        }
    return summary


def summarise_values(values_by_material):
    out = {}
    for material in sorted(values_by_material):
        v = values_by_material[material]
        finite = v[np.isfinite(v)]
        if finite.size == 0:
            out[material] = {"p05": None, "p50": None, "p95": None}
            continue
        out[material] = {
            "p05": round(float(np.percentile(finite, 5)), 4),
            "p50": round(float(np.percentile(finite, 50)), 4),
            "p95": round(float(np.percentile(finite, 95)), 4),
        }
    return out


def pairwise(values_by_material, direction):
    names = sorted(values_by_material)
    pairs = []
    for a, b in itertools.combinations(names, 2):
        va, vb = values_by_material[a], values_by_material[b]
        ok = np.isfinite(va) & np.isfinite(vb)
        if not ok.any():
            continue
        better = (va[ok] > vb[ok]) if direction > 0 else (va[ok] < vb[ok])
        p = float(np.mean(better))
        pairs.append({
            "a": a, "b": b,
            "prob_a_better": round(p, 4),
            "verdict": ("indistinguishable"
                        if INDISTINGUISHABLE_BAND[0] <= p <= INDISTINGUISHABLE_BAND[1]
                        else ("a_better" if p > 0.5 else "b_better")),
            "decisive": bool(p >= 0.95 or p <= 0.05),
        })
    return pairs


# ----------------------------------------------------------------------------
# Deterministic scenario enumeration over source disagreements
# ----------------------------------------------------------------------------
def disagreement_groups(materials):
    groups = []
    for material in sorted(materials):
        for prop in sorted(materials[material]["properties"]):
            rows = materials[material]["properties"][prop]
            if len(rows) > 1:
                groups.append((material, prop, rows))
    return groups


def deterministic_fom(entry, material, overrides, cold, hot, scenario_kind):
    """Nominal volumetric energy density for one material under overrides."""
    props = entry["properties"]

    def val(prop):
        if (material, prop) in overrides and prop in props:
            return overrides[(material, prop)]
        if prop not in props:
            return None
        return float(np.mean([r["value"] for r in props[prop]]))

    density, cp = val("density"), val("heat_capacity")
    if density is None or cp is None:
        return None
    t_lo, t_hi, _, _ = derive_window(val("min_operating_temperature"),
                                     val("melting_point"),
                                     val("max_operating_temperature"),
                                     entry["family"])
    if t_lo is None:
        return None
    if scenario_kind == "common_window":
        if t_lo > cold + 1e-9 or t_hi < hot - 1e-9:
            return None
        t_lo, t_hi = cold, hot
    latent, melt = val("latent_heat"), val("melting_point")
    specific = cp * max(float(t_hi) - float(t_lo), 0.0)
    if latent is not None and melt is not None and t_lo <= melt <= t_hi:
        specific += latent
    return density * specific / 3.6e6


def scenario_rankings(materials, groups, cold, hot, scenario_kind):
    if not groups:
        return {"n_source_choice_groups": 0,
                "note": "no competing source rows in the input table"}
    sizes = [len(rows) for _, _, rows in groups]
    total = int(np.prod([float(s) for s in sizes]))
    orderings = defaultdict(list)
    seen = 0
    for combo in itertools.product(*[range(s) for s in sizes]):
        if seen >= MAX_SCENARIO_COMBINATIONS:
            break
        seen += 1
        overrides, label_bits = {}, []
        for (material, prop, rows), idx in zip(groups, combo):
            overrides[(material, prop)] = rows[idx]["value"]
            label_bits.append("%s/%s=%g" % (material, prop, rows[idx]["value"]))
        values = {}
        for material in sorted(materials):
            fom = deterministic_fom(materials[material], material, overrides,
                                    cold, hot, scenario_kind)
            if fom is not None:
                values[material] = fom
        ordered = tuple(m for m, _ in sorted(values.items(),
                                             key=lambda kv: (-kv[1], kv[0])))
        orderings[ordered].append(" ; ".join(label_bits))
    return {
        "n_source_choice_groups": len(groups),
        "source_choice_groups": ["%s / %s (%d competing rows)" % (m, p, len(r))
                                 for m, p, r in groups],
        "n_scenarios_possible": total,
        "n_scenarios_evaluated": seen,
        "truncated": total > MAX_SCENARIO_COMBINATIONS,
        "n_distinct_orderings": len(orderings),
        "distinct_orderings": [
            {"order_best_to_worst": list(order),
             "n_scenarios": len(examples),
             "example_source_choice": examples[0]}
            for order, examples in sorted(orderings.items(),
                                          key=lambda kv: (-len(kv[1]), kv[0]))],
    }


# ----------------------------------------------------------------------------
# Scenario driver
# ----------------------------------------------------------------------------
def run_scenario(materials, rng, draws, scenario_kind, cold, hot):
    included, excluded = {}, {}
    fom_samples = defaultdict(dict)
    material_detail = {}

    for material in sorted(materials):
        entry = materials[material]
        props = entry["properties"]
        missing = [p for p in REQUIRED_FOR_ENERGY if p not in props]
        if missing:
            excluded[material] = ("missing required property/properties: %s"
                                  % ", ".join(missing))
            continue
        lo_nom, hi_nom = entry["nominal_window_C"]
        if lo_nom is None or hi_nom is None:
            excluded[material] = entry["window_basis"]
            continue
        if hi_nom - lo_nom <= 0:
            excluded[material] = ("derived window has zero or negative swing "
                                  "(%.0f to %.0f C)" % (lo_nom, hi_nom))
            continue

        samples = {prop: sample_property(rng, props[prop], draws)
                   for prop in sorted(props)}
        t_lo, t_hi, _, _ = derive_window(
            samples.get("min_operating_temperature"),
            samples.get("melting_point"),
            samples.get("max_operating_temperature"),
            entry["family"])

        if scenario_kind == "common_window":
            if lo_nom > cold + 1e-9 or hi_nom < hot - 1e-9:
                excluded[material] = (
                    "derived window %.0f-%.0f C does not cover the common "
                    "swing %.0f-%.0f C" % (lo_nom, hi_nom, cold, hot))
                continue
            t_lo = np.full(draws, cold)
            t_hi = np.full(draws, hot)

        foms = compute_foms(
            density=samples["density"], cp=samples["heat_capacity"],
            cond=samples.get("thermal_conductivity"),
            latent=samples.get("latent_heat"),
            melt=samples.get("melting_point"),
            t_lo=t_lo, t_hi=t_hi,
            cost_per_tonne=samples.get("cost_per_tonne"))

        included[material] = entry
        for name in FOM_DIRECTION:
            if name in foms:
                fom_samples[name][material] = foms[name]

        material_detail[material] = {
            "family": entry["family"],
            "window_C_mean": [round(float(np.mean(t_lo)), 1),
                              round(float(np.mean(t_hi)), 1)],
            "delta_T_C_mean": round(float(np.mean(foms["_delta_T_C"])), 1),
            "window_basis": entry["window_basis"],
            "cold_end_is_an_assumption": entry["assumed_cold_end"],
            "properties_used": sorted(props),
            "n_competing_source_rows": {p: len(props[p]) for p in sorted(props)
                                        if len(props[p]) > 1},
            "latent_heat_counted_fraction_of_draws":
                foms["_latent_counted_fraction"],
        }

    results = {
        "scenario": scenario_kind,
        "swing_C": ([cold, hot] if scenario_kind == "common_window"
                    else "per material, own derived window"),
        "n_materials_included": len(included),
        "materials_included": sorted(included),
        "materials_excluded": excluded,
        "materials_with_assumed_cold_end": sorted(
            m for m in included if materials[m]["assumed_cold_end"]),
        "material_detail": material_detail,
        "figures_of_merit": {},
    }

    for name, direction in FOM_DIRECTION.items():
        values = fom_samples.get(name, {})
        if len(values) < 2:
            results["figures_of_merit"][name] = {
                "note": "fewer than two materials carry this figure of merit",
                "materials": sorted(values)}
            continue
        ranks = rank_matrix(values, direction)
        medians = {m: float(np.median(values[m][np.isfinite(values[m])]))
                   if np.isfinite(values[m]).any() else np.inf
                   for m in values}
        ordering = sorted(values, key=lambda m: (-direction * medians[m], m))
        results["figures_of_merit"][name] = {
            "direction": ("higher_is_better" if direction > 0
                          else "lower_is_better"),
            "n_materials": len(values),
            "median_value_ordering_best_to_worst": ordering,
            "value_percentiles": summarise_values(values),
            "rank_distribution": summarise_ranks(ranks, len(values)),
            "pairwise": pairwise(values, direction),
        }
    return results


def stability_summary(scenario_result):
    out = {}
    for name, block in scenario_result["figures_of_merit"].items():
        if "rank_distribution" not in block:
            continue
        firm, soft = [], []
        for material, stats in block["rank_distribution"].items():
            if stats["rank_p05"] == stats["rank_p95"]:
                firm.append({"material": material, "rank": stats["modal_rank"]})
            else:
                soft.append({"material": material,
                             "modal_rank": stats["modal_rank"],
                             "rank_p05_p95": [stats["rank_p05"],
                                              stats["rank_p95"]]})
        ties = [p for p in block["pairwise"]
                if p["verdict"] == "indistinguishable"]
        out[name] = {
            "materials_holding_one_rank_in_90pct_of_draws":
                sorted(firm, key=lambda d: d["rank"]),
            "materials_whose_rank_moves":
                sorted(soft, key=lambda d: d["modal_rank"]),
            "indistinguishable_pairs": ties,
            "n_indistinguishable_pairs": len(ties),
            "n_decisive_pairs": sum(1 for p in block["pairwise"]
                                    if p["decisive"]),
            "n_pairs_total": len(block["pairwise"]),
        }
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(
        description="Monte Carlo screening of thermal energy storage media.")
    parser.add_argument("--input", default=None,
                        help="property CSV; default: tes_property_table.csv if "
                             "present, else tes_screening_demo_inputs.csv")
    parser.add_argument("--output", default=os.path.join(
        here, "tes_screening_output.json"))
    parser.add_argument("--draws", type=int, default=DEFAULT_DRAWS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--cold", type=float, default=DEFAULT_COLD_C)
    parser.add_argument("--hot", type=float, default=DEFAULT_HOT_C)
    args = parser.parse_args()

    curated = os.path.join(here, "tes_property_table.csv")
    demo = os.path.join(here, "tes_screening_demo_inputs.csv")
    if args.input:
        path = (args.input if os.path.isabs(args.input)
                else os.path.join(here, args.input))
        base = os.path.basename(path)
        kind = ("curated" if base == "tes_property_table.csv"
                else ("demo_placeholder"
                      if base == "tes_screening_demo_inputs.csv"
                      else "user_supplied"))
    elif os.path.exists(curated):
        path, kind = curated, "curated"
    elif os.path.exists(demo):
        path, kind = demo, "demo_placeholder"
    else:
        raise SystemExit("no input table found in %s" % here)
    if not os.path.exists(path):
        raise SystemExit("input table not found: %s" % path)

    with open(path, "rb") as handle:
        digest = hashlib.sha256(handle.read()).hexdigest()

    rows = read_table(path)
    materials = build_materials(rows)
    groups = disagreement_groups(materials)

    rng = np.random.default_rng(args.seed)
    own = run_scenario(materials, rng, args.draws, "own_window",
                       args.cold, args.hot)
    common = run_scenario(materials, rng, args.draws, "common_window",
                          args.cold, args.hot)

    output = {
        "meta": {
            "script": "tes_screening.py",
            "script_version": "1.2",
            "input_file": os.path.basename(path),
            "input_kind": kind,
            "input_sha256": digest,
            "input_is_demo_placeholder": kind == "demo_placeholder",
            "seed": args.seed,
            "monte_carlo_draws": args.draws,
            "default_relative_uncertainty_when_no_interval":
                DEFAULT_REL_UNCERTAINTY,
            "freeze_margin_C": FREEZE_MARGIN_C,
            "pcm_half_swing_C": PCM_HALF_SWING_C,
            "solid_headroom_below_melt_C": SOLID_HEADROOM_C,
            "assumed_cold_end_C_when_no_data": DEFAULT_COLD_END_C,
            "common_window_C": [args.cold, args.hot],
            "indistinguishable_band_on_pairwise_probability":
                list(INDISTINGUISHABLE_BAND),
            "sampling_model": ("per draw: pick one competing source row "
                               "uniformly per material+property, then sample "
                               "uniformly inside that row's [low, high]; if no "
                               "interval is given, value*(1 +/- 0.10)"),
            "deterministic": ("fixed seed, sorted iteration, no timestamps; "
                              "identical inputs give byte-identical output"),
            "environment": {
                "python": sys.version.split()[0],
                "platform": platform.platform(),
                "numpy": np.__version__,
                "scipy": SCIPY_VERSION,
            },
            "canonical_units": CANONICAL_UNITS,
            "n_input_rows_parsed": len(rows),
            "n_materials_in_table": len(materials),
        },
        "figure_of_merit_definitions": {
            "energy_density_kWh_per_m3":
                "rho * (cp * dT + L if the melting point falls inside dT) "
                "/ 3.6e6. Volumetric, media only, no void fraction or vessel.",
            "cost_per_kWh_th_USD":
                "cost_per_tonne / (1000 kg * (cp*dT + L) / 3.6e6). Media cost "
                "only, mixed and mostly unstated price-basis years.",
            "diffusivity_mm2_per_s":
                "1e6 * k / (rho * cp). A conduction-limited charge-rate proxy "
                "only; meaningless where convection or particle transport sets "
                "the rate.",
        },
        "limitations": [
            "Media-only screening. Containment, vessel, insulation, pumps, "
            "heat exchangers and system integration are excluded, and they "
            "usually dominate installed cost per kWh_th.",
            "No exergy or temperature-quality weighting: a kWh stored at 700 C "
            "counts the same as a kWh stored at 300 C.",
            "No degradation term. Cycle life, corrosion allowance, salt "
            "make-up and particle attrition are not represented, so media "
            "known to degrade are not penalised here.",
            "Properties are constant over the whole swing; real cp, rho and k "
            "vary with temperature, and the script picks the tabulated "
            "temperature nearest the window midpoint.",
            "Cold-end temperatures are system design choices. Where the table "
            "gives none, the script assumes one; those materials are listed "
            "and their energy densities are assumption-driven.",
            "Uniform sampling inside [low, high] is a deliberately weak "
            "uncertainty model, not a calibrated posterior; its width is "
            "whatever the table compiler entered.",
            "Competing sources are drawn with equal probability; no source is "
            "judged more reliable than another.",
            "Volumetric figures mix bulk packed-bed density with monolith and "
            "melt density, so cross-family volumetric comparison is coarse.",
            "Ranks are computed within the set of materials in the input "
            "table; adding or removing a material changes every rank.",
        ],
        "scenarios": {"own_window": own, "common_window": common},
        "stability_summary": {"own_window": stability_summary(own),
                              "common_window": stability_summary(common)},
        "source_disagreement_scenarios": {
            "note": ("deterministic rankings on volumetric energy density for "
                     "every combination of competing source rows, nominal "
                     "values only, no sampling"),
            "own_window": scenario_rankings(materials, groups, args.cold,
                                            args.hot, "own_window"),
            "common_window": scenario_rankings(materials, groups, args.cold,
                                               args.hot, "common_window"),
        },
        "warnings": WARNINGS,
    }

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2, sort_keys=False)
        handle.write("\n")

    # ---------------- console summary ----------------
    print("input: %s (%s)" % (os.path.basename(path), kind))
    print("sha256: %s" % digest)
    print("rows parsed: %d   materials in table: %d"
          % (len(rows), len(materials)))
    print("seed=%d draws=%d  python=%s numpy=%s scipy=%s"
          % (args.seed, args.draws, sys.version.split()[0], np.__version__,
             SCIPY_VERSION))
    for scenario_name in ("own_window", "common_window"):
        block = output["scenarios"][scenario_name]
        print("\n=== scenario: %s (swing: %s) ==="
              % (scenario_name, block["swing_C"]))
        print("included: %d   excluded: %d"
              % (block["n_materials_included"],
                 len(block["materials_excluded"])))
        for material, reason in sorted(block["materials_excluded"].items()):
            print("  excluded: %s -- %s" % (material, reason))
        if block["materials_with_assumed_cold_end"]:
            print("  cold end assumed (%.0f C) for: %s"
                  % (DEFAULT_COLD_END_C,
                     ", ".join(block["materials_with_assumed_cold_end"])))
        for fom_name, fom in block["figures_of_merit"].items():
            if "rank_distribution" not in fom:
                print("  -- %s: %s" % (fom_name, fom["note"]))
                continue
            print("  -- %s (%s, %d materials)"
                  % (fom_name, fom["direction"], fom["n_materials"]))
            for material in fom["median_value_ordering_best_to_worst"]:
                stats = fom["rank_distribution"][material]
                pct = fom["value_percentiles"][material]
                print("     %-44s median=%10.3f  p05-p95=%9.3f-%10.3f  "
                      "modal_rank=%2d  rank p05-p95=%2d-%2d  P(1st)=%.3f"
                      % (material[:44], pct["p50"], pct["p05"], pct["p95"],
                         stats["modal_rank"], stats["rank_p05"],
                         stats["rank_p95"], stats["prob_rank_1"]))
            ties = [p for p in fom["pairwise"]
                    if p["verdict"] == "indistinguishable"]
            print("     indistinguishable pairs (0.40 <= P(a better) <= 0.60):"
                  " %d of %d" % (len(ties), len(fom["pairwise"])))
            for pair in ties:
                print("        %s  vs  %s   P(a better)=%.3f"
                      % (pair["a"], pair["b"], pair["prob_a_better"]))
    for scenario_name in ("own_window", "common_window"):
        block = output["source_disagreement_scenarios"][scenario_name]
        if "n_distinct_orderings" not in block:
            print("\nsource-choice scenarios (%s): %s"
                  % (scenario_name, block.get("note")))
            continue
        print("\nsource-choice scenarios (%s): %d groups, %d of %d "
              "combinations evaluated, %d distinct orderings"
              % (scenario_name, block["n_source_choice_groups"],
                 block["n_scenarios_evaluated"], block["n_scenarios_possible"],
                 block["n_distinct_orderings"]))
    if WARNINGS:
        print("\nwarnings (%d):" % len(WARNINGS))
        for message in WARNINGS:
            print("  - %s" % message)
    print("\nwrote %s" % args.output)


if __name__ == "__main__":
    main()
