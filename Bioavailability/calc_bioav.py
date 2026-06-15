"""
calc_bioav.py
-------------
Given a CSV with a molecular formula column (and optionally
pre-existing CHONPS + Mass columns), compute per-compound:

  Class         -- MSCC compound class (Rivas-Ubach et al. 2018)
  stoichMet_OC  -- metabolic stoichiometry of the organic carbon donor
  stoichMet_O2  -- metabolic stoichiometry of O2
  lambda_O2_pH7     -- lambda at pH 7 (fixed reference)
  lambda0_O2        -- lambda at pH 0 (standard state)
  lambda_O2_pH<x>   -- lambda at user-specified pH (column name encodes pH)
  ne                -- electrons transferred
  TER           -- threshold elemental ratio
  f_C           -- 1 / ne
  f_CN          -- C:N nutrient limitation factor
  bioav_pH<p>_VhOC_<x>          -- exp(-|sOC| / (Vh*[OC])) * f_CN * f_C
                                  column name encodes pH and Vh×[OC]
  rO2_pH<p>_VhOC_<x>_VhO2_<y>  -- -bioav * exp(-|sO2| / (Vh_O2*[O2])) * sO2
                                  column name encodes pH, Vh×[OC], and Vh_O2×[O2]

Usage
-----
    python calc_bioav.py --input compounds.csv [options]

Required argument
-----------------
  --input       Path to input CSV.

Optional arguments
------------------
  --output      Output CSV path (default: <input>_bioav.csv).
  --molform-col Column name containing molecular formulas (default: MolForm).
  --pH          pH for lambda_pH calculation (default: 7.0).
  --Vh          Harvest volume for OC, L (default: 1.0).
  --OC          Organic carbon concentration, mol/L (default: 1.0).
                Ignored per-row when the input CSV contains an
                'OM_Concentration' column.
  --VhO2        Harvest volume for O2, L (default: 1.0).
  --O2          O2 concentration, mol/L (default: abundant → exp term = 1).

CHONPS detection
----------------
If columns C, H, O, N, P, S are ALL present they are used directly.
Otherwise they are parsed from the molecular formula column.
Mass is re-computed from monoisotopic weights when parsed or when absent.

OM_Concentration column
-----------------------
If the input CSV contains a column named 'OM_Concentration' it is used as a
compound-specific [OC]_i (mol/L) for each row. Rows where the value is missing
or zero fall back to the scalar --OC argument. When the column is used, the
bioav/rO2 output column names contain 'VhOC_colOC' instead of a fixed number.

Attributions
------------
  MSCC — Metabolomics-based Compound Classification
    Compound classification is performed using the MSCC scheme
    (Rivas-Ubach et al. 2018) developed at Pacific Northwest National
    Laboratory.
    https://github.com/PNNL-Comp-Mass-Spec/MSCC
    License: BSD 2-Clause
    https://github.com/PNNL-Comp-Mass-Spec/MSCC?tab=BSD-2-Clause-1-ov-file
"""

import argparse
import math
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Monoisotopic masses and element list
# ---------------------------------------------------------------------------
_MONO = {
    "C": 12.0000000,
    "H":  1.0078250,
    "O": 15.9949146,
    "N": 14.0030740,
    "P": 30.9737634,
    "S": 31.9720707,
}
_CHONPS = ["C", "H", "O", "N", "P", "S"]
_FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def _parse_formula(formula: str) -> dict:
    """Return element-count dict for a molecular formula string."""
    if not isinstance(formula, str) or not formula.strip():
        return {}
    counts: dict = {}
    for el, n in _FORMULA_RE.findall(formula):
        if el:
            counts[el] = counts.get(el, 0) + (int(n) if n else 1)
    return counts


def _chonps_from_formula(formula: str) -> dict:
    """Parse a formula string and return CHONPS floats + monoisotopic Mass."""
    counts = _parse_formula(formula)
    result = {el: float(counts.get(el, 0)) for el in _CHONPS}
    result["Mass"] = sum(counts.get(el, 0) * _MONO[el]
                         for el in _MONO if el in counts)
    return result


# ---------------------------------------------------------------------------
# MSCC Classification  (Rivas-Ubach et al. 2018, Anal. Chem.)
# ---------------------------------------------------------------------------
def mscc_classify(C, H, O, N, P, S, Mass) -> str:
    """Return MSCC compound-class label for a single compound."""
    try:
        C = float(C); H = float(H); O = float(O)
        N = float(N); P = float(P); S = float(S)
        Mass = 0.0 if pd.isna(Mass) else float(Mass)
    except (TypeError, ValueError):
        return "Unclassified"
    if math.isnan(C) or C <= 0:
        return "Unclassified"

    OC = O / C
    HC = H / C
    NC = N / C
    PC = P / C if C > 0 else 0.0
    NP = N / P if P > 0 else 0.0  # MSCC convention: NP = 0 when P = 0

    lipid = (OC <= 0.6   and HC >= 1.32  and NC <= 0.126
             and PC < 0.35  and NP <= 5)
    carb  = (OC >= 0.8   and HC >= 1.65  and HC < 2.7   and N == 0)
    amino = (OC >= 0.61  and HC >= 1.45
             and NC <= 0.2   and NC > 0.07
             and PC < 0.3    and NP <= 2
             and O >= 3      and N >= 1)
    phyto = (OC <= 1.15  and HC < 1.32   and NC < 0.126
             and PC <= 0.2   and NP <= 3)
    prot1 = (OC > 0.12   and OC <= 0.6
             and HC > 0.9    and HC < 2.5
             and NC >= 0.126 and NC <= 0.7
             and PC < 0.17   and N >= 1)
    prot2 = (OC > 0.6    and OC <= 1.0
             and HC > 1.2    and HC < 2.5
             and NC > 0.2    and NC <= 0.7
             and PC < 0.17   and N >= 1)
    nuc   = (OC >= 0.5   and OC < 1.7
             and HC > 1.0    and HC < 1.8
             and NC >= 0.2   and NC <= 0.5
             and PC >= 0.1   and PC <= 0.35
             and NP > 0.6    and NP <= 5
             and N >= 2      and P >= 1
             and S == 0
             and Mass > 305  and Mass < 523)

    # Priority: Nucleotide > Carbohydrate > Lipid > AminoSugar >
    #           Phytochemical > Protein
    matches = []
    if nuc:             matches.append("Nucleotides")
    if carb:            matches.append("Carbohydrates")
    if lipid:           matches.append("Lipid")
    if amino:           matches.append("AminoSugars")
    if phyto:           matches.append("Phytochemicals")
    if prot1 or prot2:  matches.append("Protein")

    if not matches:
        return "Unclassified"
    if len(matches) == 1:
        return matches[0]
    # Nucleotide trumps Protein / AminoSugar double-matches
    if "Nucleotides" in matches and (
        "Protein" in matches or "AminoSugars" in matches
    ):
        return "Nucleotides"
    return matches[0]


# ---------------------------------------------------------------------------
# Core thermodynamic function
# ---------------------------------------------------------------------------
def getThermoStoich(chemForm, pH_user: float = 7.0):
    """
    Compute thermodynamic stoichiometry and lambda values for a compound.

    Parameters
    ----------
    chemForm : list  [C, H, N, O, P, S, charge]
    pH_user  : float  pH for lambda_pH output (default 7.0)

    Returns
    -------
    lambda_pH7 : lambda at pH 7
    lambda0    : lambda at pH 0 (standard state)
    lambda_pH  : lambda at pH_user
    stoichMet_OC  : metabolic stoich of OC donor
    stoichMet_O2  : metabolic stoich of O2
    ne         : electrons transferred
    TER        : threshold elemental ratio
    """
    a, b, c, d, e, f, z = (chemForm[0], chemForm[1], chemForm[2],
                            chemForm[3], chemForm[4], chemForm[5],
                            chemForm[6])

    # Step 1a – donor half-reaction
    yH2o = -(3*a + 4*e - d)
    yH   =   5*a + b - 4*c - 2*d + 7*e - f
    yE   =  -z   + 4*a + b - 3*c - 2*d + 5*e - 2*f
    stoichD = [-1, yH2o, a, c, e, f, yH, yE] + [0.0] * 20

    # Step 1b – O2 acceptor half-reaction (acceptor 14)
    stoichA = [0.0] * 28
    stoichA[8] = -1.0   # O2
    stoichA[6] = -4.0   # H+
    stoichA[7] = -4.0   # e-
    stoichA[1] =  2.0   # H2O

    # Step 1c – catabolic reaction
    yEd = stoichD[7]
    yEa = stoichA[7]
    stoichCat = [stoichD[i] - (yEd / yEa) * stoichA[i] for i in range(28)]

    # Step 2a – anabolic biomass target: CH1.8 O0.5 N0.2
    aB, bB, cB, dB, eB, fB, zB = 1, 1.8, 0.2, 0.5, 0, 0, 0
    stoichAnStarB = [
        -1,
        -(3*aB + 4*eB - dB),
         aB, cB, eB, fB,
         5*aB + bB - 4*cB - 2*dB + 7*eB - fB,
        -zB   + 4*aB + bB - 3*cB - 2*dB + 5*eB - 2*fB,
    ] + [0.0] * 20
    stoichAnStarB = [-x for x in stoichAnStarB]
    stoichAnStarB[27] = stoichAnStarB[0]
    stoichAnStarB[0]  = 0.0

    # Step 2b – overall anabolic reaction
    stoichAnStar = [stoichAnStarB[i] + (1.0 / a) * stoichD[i]
                    for i in range(28)]
    yEana = stoichAnStar[7]
    if yEana > 0:
        stoichAn = [stoichAnStar[i] - (yEana / yEa) * stoichA[i]
                    for i in range(28)]
    elif yEana < 0:
        stoichAn = [stoichAnStar[i] - (yEana / yEd) * stoichD[i]
                    for i in range(28)]
    else:
        stoichAn = list(stoichAnStar)

    # Step 3 – Gibbs energies
    ne   = -z + 4*a + b - 3*c - 2*d + 5*e - 2*f
    nosc = -ne / a + 4
    delGcox0 = 60.3 - 28.5 * nosc
    delGd0   = delGcox0 * a * abs(stoichD[0])

    delGf0 = [0, -237.2, -586.9, -79.37, -1089.1, 12.05, 0, 0,
              16.5, -111.3, -32.2, 18.19, -4.6, -78.87, 0,
              -744.63, -33.4, -486.6, 0, 522.5, -690, -489.8,
              -1012.6, -34.06, -392, -465.14, -228, -67]
    delGcox0_zero = float(np.dot(delGf0, stoichD))
    delGf0[0] = (delGd0 - delGcox0_zero) / stoichD[0]

    delGcat0 = float(np.dot(delGf0, stoichCat))
    delGan0  = float(np.dot(delGf0, stoichAn))

    if math.isnan(delGan0) or math.isnan(delGcat0):
        nan6 = (float("nan"),) * 7
        return nan6

    R = 0.008314   # kJ / (mol·K)
    T = 298.15     # K
    eta     = 0.43
    delGsyn = 200.0

    def _lambda_at(H_activity: float):
        """Compute lambda given [H+] activity."""
        dGcat = delGcat0 + R * T * stoichCat[6] * math.log(H_activity)
        dGan  = delGan0  + R * T * stoichAn[6]  * math.log(H_activity)
        if math.isnan(dGan):
            return float("nan")
        m = 1 if dGan < 0 else -1
        lam = (dGan * (eta ** m) + delGsyn) / (-dGcat * eta)
        return lam

    lambda0   = _lambda_at(1.0)          # pH 0  (H+ activity = 1)
    lambda_7  = _lambda_at(1e-7)         # pH 7
    lambda_pH = _lambda_at(10.0 ** (-pH_user))  # user pH

    # Metabolic stoichiometry at user-specified pH
    lam_for_stoich = lambda_pH
    if lam_for_stoich > 0:
        stoichMet = [lam_for_stoich * stoichCat[i] + stoichAn[i]
                     for i in range(28)]
    else:
        stoichMet = list(stoichAn)

    stoichMet_OC   = stoichMet[0]
    stoichMet_O2   = stoichMet[8]
    stoichMet_biom = stoichMet[27]

    # TER
    try:
        CUE = stoichMet_biom * 1.0 / (abs(stoichMet_OC) * a)
        NUE = stoichMet_biom * 0.2 / (abs(stoichMet_OC) * c
                                       + abs(stoichMet[3]) * 1.0)
        TER = (NUE / CUE) * (1.0 / 0.2) if CUE != 0 else float("nan")
    except (ZeroDivisionError, ValueError):
        TER = float("nan")

    return lambda_7, lambda0, lambda_pH, stoichMet_OC, stoichMet_O2, ne, TER


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------
def process(input_path: str, output_path: str, molform_col: str,
            pH: float, Vh: float, OC: float,
            VhO2: float = 1.0, O2: float = None) -> pd.DataFrame:
    """Load CSV, compute all columns, return and save result.

    O2=None means O2 is abundant: exp(-|sO2| / (VhO2 * [O2])) = 1.
    """

    df = pd.read_csv(input_path, low_memory=False)
    print(f"Loaded {len(df):,} rows from {input_path}")

    # ---- CHONPS detection / parsing ---------------------------------------
    has_chonps = all(col in df.columns for col in _CHONPS)
    if has_chonps:
        print("CHONPS columns detected — using existing values.")
        for col in _CHONPS:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        if "Mass" not in df.columns:
            df["Mass"] = df[molform_col].apply(
                lambda f: (_chonps_from_formula(str(f)) or {}).get("Mass", float("nan"))
            )
            print("  Mass column absent — computed from formula.")
        else:
            df["Mass"] = pd.to_numeric(df["Mass"], errors="coerce")
    else:
        print(f"CHONPS columns not found — parsing from '{molform_col}'.")
        if molform_col not in df.columns:
            sys.exit(f"ERROR: column '{molform_col}' not found in input CSV.")
        for col in _CHONPS + ["Mass"]:
            df[col] = float("nan")
        for idx, formula in df[molform_col].items():
            parsed = _chonps_from_formula(str(formula))
            if parsed and parsed.get("Mass", 0) > 0:
                for col in _CHONPS + ["Mass"]:
                    df.at[idx, col] = parsed[col]

    # Fill any remaining NaN CHONPS with 0 for elements other than C
    for col in ["H", "O", "N", "P", "S"]:
        df[col] = df[col].fillna(0.0)

    # ---- MSCC Class --------------------------------------------------------
    print("Assigning MSCC compound classes …")
    if molform_col not in df.columns:
        df["Class"] = "Unclassified"
    else:
        df["Class"] = df.apply(
            lambda r: mscc_classify(r["C"], r["H"], r["O"],
                                    r["N"], r["P"], r["S"], r.get("Mass", 0)),
            axis=1,
        )
    print("  Class counts:")
    for cls, cnt in df["Class"].value_counts().items():
        print(f"    {cls:20s}  {cnt:6,}")

    # ---- OC concentration: per-row column or scalar ----------------------
    OM_COL = "OM_Concentration"
    if OM_COL in df.columns:
        df[OM_COL] = pd.to_numeric(df[OM_COL], errors="coerce")
        # Fall back to scalar --OC where column is missing or zero
        oc_series = df[OM_COL].where(
            df[OM_COL].notna() & (df[OM_COL] > 0), other=OC
        )
        per_row_oc = True
        n_col = (df[OM_COL].notna() & (df[OM_COL] > 0)).sum()
        print(f"OM_Concentration column found: {n_col:,} rows use column value, "
              f"{len(df)-n_col:,} fall back to --OC={OC}")
    else:
        oc_series  = pd.Series(OC, index=df.index)
        per_row_oc = False
        print(f"No OM_Concentration column — using scalar [OC]={OC} mol/L")

    # ---- Thermodynamics ----------------------------------------------------
    VhOC_scalar = Vh * OC   # used only for column name when scalar
    VhO2OC = None if O2 is None else VhO2 * O2
    o2_desc   = "abundant" if VhO2OC is None else f"{VhO2OC:.4g}"
    vhoc_desc = f"Vh×colOC" if per_row_oc else f"{VhOC_scalar:.4g}"
    print(f"Computing thermodynamics (pH={pH}, Vh×[OC]={vhoc_desc}, "
          f"Vh_O2×[O2]={o2_desc}) …")

    # Column names encode pH and the parameter products used
    # When per-row OC is used, encode 'colOC' instead of a fixed number.
    vhoc_tag  = "colOC" if per_row_oc else f"{VhOC_scalar:.4g}"
    col_lambda_pH = f"lambda_O2_pH{pH}"
    col_bioav     = f"bioav_pH{pH}_VhOC_{vhoc_tag}"
    col_rO2       = (f"rO2_pH{pH}_VhOC_{vhoc_tag}_VhO2_abund"
                     if VhO2OC is None else
                     f"rO2_pH{pH}_VhOC_{vhoc_tag}_VhO2_{VhO2OC:.4g}")
    print(f"  Output columns: '{col_lambda_pH}', '{col_bioav}', '{col_rO2}'")

    out_cols = ["lambda_O2_pH7", "lambda0_O2", col_lambda_pH,
                "stoichMet_OC", "stoichMet_O2", "ne", "TER",
                "f_C", "f_CN", col_bioav, col_rO2]
    for col in out_cols:
        df[col] = float("nan")

    eligible = (
        df["C"].notna() & (df["C"] > 0) &
        df[_CHONPS].notna().all(axis=1)
    )
    print(f"  Eligible rows (C > 0, CHONPS complete): {eligible.sum():,}")

    n_ok = n_err = 0
    for idx in df.index[eligible]:
        try:
            chemForm = [df.at[idx, "C"], df.at[idx, "H"], df.at[idx, "N"],
                        df.at[idx, "O"], df.at[idx, "P"], df.at[idx, "S"],
                        0.0]  # neutral charge
            lam7, lam0, lam_pH, sOC, sO2, ne, TER = getThermoStoich(
                chemForm, pH_user=pH
            )

            # f_C = 1 / ne
            f_C = (1.0 / ne) if (math.isfinite(ne) and ne != 0) else float("nan")

            # f_CN = exp(-alpha * (x - 1)^2),  x = (C / (N+1)) / TER
            C_val = df.at[idx, "C"]
            N_val = df.at[idx, "N"]
            if math.isfinite(TER) and TER != 0:
                x    = (C_val / (N_val + 1.0)) / TER
                f_CN = math.exp(-0.780057283995686 * (x - 1.0) ** 2)
            else:
                f_CN = 0.0

            # rO2_OClim and bioav use per-row or scalar Vh * [OC]
            fC = f_C if (isinstance(f_C, float) and math.isfinite(f_C)) \
                else float("nan")
            VhOC_i = Vh * float(oc_series.at[idx])
            exp_O2 = (1.0 if VhO2OC is None
                      else math.exp(-abs(sO2) / VhO2OC))
            bioav  = math.exp(-abs(sOC) / VhOC_i) * f_CN * fC
            # rO2 = -bioav * exp(-|sO2| / (Vh_O2*[O2])) * sO2
            rO2    = -bioav * exp_O2 * sO2

            df.at[idx, "lambda_O2_pH7"] = lam7
            df.at[idx, "lambda0_O2"]    = lam0
            df.at[idx, col_lambda_pH]   = lam_pH
            df.at[idx, "stoichMet_OC"]  = sOC
            df.at[idx, "stoichMet_O2"] = sO2
            df.at[idx, "ne"]           = ne
            df.at[idx, "TER"]          = TER
            df.at[idx, "f_C"]          = fC
            df.at[idx, "f_CN"]         = f_CN
            df.at[idx, col_bioav] = bioav
            df.at[idx, col_rO2]  = rO2
            n_ok += 1

        except Exception as exc:
            n_err += 1
            if n_err <= 3:
                print(f"  [err sample idx={idx}]: {type(exc).__name__}: {exc}")

    print(f"  Computed: {n_ok:,}  |  errors skipped: {n_err}")

    for col in out_cols:
        vals = df[col].dropna()
        if len(vals):
            print(f"  {col:20s}  range: [{vals.min():.5g}, {vals.max():.5g}]")

    # ---- Save --------------------------------------------------------------
    df.to_csv(output_path, index=False)
    print(f"\nSaved -> {output_path}")
    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Calculate MSCC class, lambda, stoichiometry, and "
                    "bioavailability for compounds in a CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",  required=True,
                   help="Path to input CSV file.")
    p.add_argument("--output", default=None,
                   help="Path for output CSV (default: <input>_bioav.csv).")
    p.add_argument("--molform-col", default="MolForm",
                   dest="molform_col",
                   help="Column name containing molecular formulas.")
    p.add_argument("--pH",  type=float, default=7.0,
                   help="pH for lambda_pH calculation.")
    p.add_argument("--Vh",  type=float, default=1.0,
                   help="Harvest volume (L).")
    p.add_argument("--OC",  type=float, default=1.0,
                   help="Organic carbon concentration (mol/L).")
    p.add_argument("--VhO2", type=float, default=1.0,
                   help="Harvest volume for O2 (L).")
    p.add_argument("--O2",  type=float, default=None,
                   help="O2 concentration (mol/L). Omit for O2-abundant "
                        "(exp term = 1).")
    return p


if __name__ == "__main__":
    args = _build_parser().parse_args()

    out = args.output
    if out is None:
        from pathlib import Path
        p = Path(args.input)
        out = str(p.parent / (p.stem + "_bioav.csv"))

    o2_str = (f"[O2]={args.O2} mol/L  →  Vh_O2×[O2]={args.VhO2 * args.O2:.4g}"
              if args.O2 is not None else "[O2]=abundant")
    print(f"Parameters: pH={args.pH}, Vh={args.Vh} L, "
          f"[OC]={args.OC} mol/L  →  Vh×[OC]={args.Vh * args.OC:.4g}  |  "
          + o2_str)

    process(
        input_path=args.input,
        output_path=out,
        molform_col=args.molform_col,
        pH=args.pH,
        Vh=args.Vh,
        OC=args.OC,
        VhO2=args.VhO2,
        O2=args.O2,
    )
