#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
extract_lvm_target_spectra.py

Reads the match file produced by the drpall cross-match (e.g., candidatos_lvm.matched.txt), 
groups entries by exposure, temporarily downloads the corresponding lvmSFrame file from the SAS, 
identifies the fibers closest to each target, saves only those spectra in small .npz files, 
and removes the temporary SFrame file afterward.

Basado en la estructura de lvmSFrame usada en el notebook tutorial:
- SLITMAP
- WAVE
- FLUX
- IVAR

Usage:

python extract_lvm_target_spectra.py \
    --matched candidatos_lvm.matched.txt \
    --outdir extracted_spectra \
    --radius-arcsec 35 \
    --nclosest 5


Optional argument to specify the redux version:

--redux-version 1.2.1


Outputs:

.npz files containing spectra for each target and exposure

a CSV summary with all selected fibers
"""

import os
import re
import json
import shutil
import argparse
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from astropy.io import fits, ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

# ============================================================
# CREDENCIALES SAS
# ============================================================

USER = "user for SAS"
PASSWORD = "password for SAS"

# ============================================================
# UTILIDADES
# ============================================================

def clean_str(x):
    """Convierte bytes/strings variados a str limpia."""
    if isinstance(x, bytes):
        x = x.decode("utf-8", errors="ignore")
    return str(x).strip()


def find_column(columns, candidates, required=True):
    """
    Search for a column in a list of names, ignores caps
    """
    lowmap = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lowmap:
            return lowmap[cand.lower()]
    if required:
        raise KeyError(f"Could not find any of these columns: {candidates}")
    return None


def normalize_tileid(tileid):
    """
    Convert tileid to int if possible
    """
    if pd.isna(tileid):
        return None
    s = str(tileid).strip()
    s = re.sub(r"\.0$", "", s)
    try:
        return int(s)
    except Exception:
        return None


def tile_group_from_tileid(tileid):
    """
    Builds a group tipe 1044XX from tile ID .e.g 10444086
    """
    tileid = normalize_tileid(tileid)
    if tileid is None:
        return None
    s = f"{tileid:d}"
    if len(s) < 4:
        return None
    return s[:4] + "XX"


def build_sframe_url(expnum, mjd, tileid, redux_version="1.2.1"):
    """
    Builds and URL from the lvmSFrame from expnum, mjd, tileid
    """
    tileid = normalize_tileid(tileid)
    if tileid is None:
        raise ValueError(f"tileid inválido: {tileid}")
    group = tile_group_from_tileid(tileid)
    expnum = int(expnum)
    mjd = int(float(mjd))

    fname = f"lvmSFrame-{expnum:08d}.fits"
    base = "https://data.sdss5.org/sas/sdsswork/lvm/spectro/redux"
    url = f"{base}/{redux_version}/{group}/{tileid}/{mjd}/{fname}"
    return url, fname


session = requests.Session()
session.auth = (USER, PASSWORD)


def download_file(url, dest_path, chunk_size=1024 * 1024):

    r = session.get(
        url,
        stream=True,
        timeout=120
    )

    if r.status_code == 401:
        raise RuntimeError(
            "Failed authentification, check USER and PASSWORD"
        )

    r.raise_for_status()

    with open(dest_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)

def read_matched_table(path):
    """
    Reads the candidate.txt file produced by find_obs.py
    """
    tab = ascii.read(path, format="fixed_width_two_line")
    df = tab.to_pandas()

    # Normalización suave de nombres esperados
    # Target catalog:
    source_name_col = find_column(df.columns, ["Source_name"], required=False)
    ra_col = find_column(df.columns, ["RA"], required=False)
    dec_col = find_column(df.columns, ["Dec"], required=False)

    # Obs catalog:
    expnum_col = find_column(df.columns, ["expnum", "expnum_obs"], required=False)
    mjd_col = find_column(df.columns, ["mjd", "mjd_obs"], required=False)
    tileid_col = find_column(df.columns, ["tileid", "tileid_obs"], required=False)
    location_col = find_column(df.columns, ["location", "location_obs"], required=False)
    sep_col = find_column(df.columns, ["separation", "separation_arcsec"], required=False)

    rename_map = {}
    if source_name_col and source_name_col != "Source_name":
        rename_map[source_name_col] = "Source_name"
    if ra_col and ra_col != "RA":
        rename_map[ra_col] = "RA"
    if dec_col and dec_col != "Dec":
        rename_map[dec_col] = "Dec"
    if expnum_col and expnum_col != "expnum":
        rename_map[expnum_col] = "expnum"
    if mjd_col and mjd_col != "mjd":
        rename_map[mjd_col] = "mjd"
    if tileid_col and tileid_col != "tileid":
        rename_map[tileid_col] = "tileid"
    if location_col and location_col != "location":
        rename_map[location_col] = "location"
    if sep_col and sep_col != "separation":
        rename_map[sep_col] = "separation"

    if rename_map:
        df = df.rename(columns=rename_map)

    required = ["Source_name", "RA", "Dec", "expnum", "mjd", "tileid"]
    missing = [x for x in required if x not in df.columns]
    if missing:
        raise KeyError(f"Faltan columnas requeridas en matched.txt: {missing}")

    # numerical cleaning
    df["RA"] = pd.to_numeric(df["RA"], errors="coerce")
    df["Dec"] = pd.to_numeric(df["Dec"], errors="coerce")
    df["expnum"] = pd.to_numeric(df["expnum"], errors="coerce").astype("Int64")
    df["mjd"] = pd.to_numeric(df["mjd"], errors="coerce").astype("Int64")
    df["tileid"] = pd.to_numeric(df["tileid"], errors="coerce").astype("Int64")

    if "separation" in df.columns:
        try:
            df["separation"] = pd.to_numeric(df["separation"], errors="coerce")
        except Exception:
            pass

    df = df.dropna(subset=["RA", "Dec", "expnum", "mjd", "tileid"]).copy()
    return df


def find_fiber_coord_columns(df):
    """
    Search RA/DEC columns in SLITMAP
    """
    ra_col = find_column(df.columns, ["ra", "RA", "racen", "raobs"], required=False)
    dec_col = find_column(df.columns, ["dec", "DEC", "deccen", "decobs"], required=False)

    if ra_col is None or dec_col is None:
        raise KeyError(
            f"Could not find column coordinate in SLITMAP."
            f"Available columns: {list(df.columns)}"
        )
    return ra_col, dec_col


def load_good_science_fibers(sframe_path):
    """
    Opens lvmSFrame and returns
    - good science fibers
    - wave
    - flux
    - error
    """
    with fits.open(sframe_path, memmap=True) as hdul:
        # HDUs esperados según el notebook
        slit = Table(hdul["SLITMAP"].data)
        wave = np.array(hdul["WAVE"].data, dtype=float)
        flux_all = np.array(hdul["FLUX"].data, dtype=float)
        ivar_all = np.array(hdul["IVAR"].data, dtype=float)

        # Filtrado como en el tutorial
        targettype = np.array([clean_str(x).lower() for x in slit["targettype"]])
        fibstatus = np.array(slit["fibstatus"])

        good = (targettype == "science") & (fibstatus == 0)

        slit_good = slit[good].to_pandas().reset_index(drop=True)
        flux = flux_all[good]
        ivar = ivar_all[good]

        error = np.full_like(flux, np.nan, dtype=float)
        m = ivar > 0
        error[m] = 1.0 / np.sqrt(ivar[m])

    return slit_good, wave, flux, error


def extract_target_spectra_from_exposure(
    sframe_path,
    expnum,
    tileid,
    mjd,
    targets_df,
    outdir,
    radius_arcsec=35.0,
    nclosest=5,
):
    """
    Extracts nearby spectra for all targets associated with a given exposure.
    """
    slit, wave, flux, error = load_good_science_fibers(sframe_path)

    # Columnas esperadas en SLITMAP
    fiberid_col = find_column(slit.columns, ["fiberid", "fiber", "orig_ifulabel"], required=False)
    if fiberid_col is None:
        raise KeyError(f"No fiberid column found in SLITMAP. Available columns: {list(slit.columns)}")

    ra_fib_col, dec_fib_col = find_fiber_coord_columns(slit)

    c_fib = SkyCoord(slit[ra_fib_col].values * u.deg, slit[dec_fib_col].values * u.deg)

    rows_summary = []

    for _, row in targets_df.iterrows():
        src = str(row["Source_name"])
        ra = float(row["RA"])
        dec = float(row["Dec"])

        c_tar = SkyCoord(ra * u.deg, dec * u.deg)
        sep = c_tar.separation(c_fib).arcsec

        order = np.argsort(sep)
        nearest = order[:nclosest]
        within = nearest[sep[nearest] <= radius_arcsec]

        # Si ninguna cae dentro del radio, al menos guarda la más cercana
        if len(within) == 0 and len(nearest) > 0:
            within = nearest[:1]

        fiber_ids = slit.iloc[within][fiberid_col].values
        ra_sel = slit.iloc[within][ra_fib_col].values
        dec_sel = slit.iloc[within][dec_fib_col].values
        sep_sel = sep[within]
        flux_sel = flux[within]
        err_sel = error[within]

        outname = f"{src}_exp{int(expnum):05d}.npz"
        outpath = os.path.join(outdir, outname)

        meta = {
            "Source_name": src,
            "target_ra_deg": ra,
            "target_dec_deg": dec,
            "expnum": int(expnum),
            "tileid": int(tileid),
            "mjd": int(mjd),
            "n_fibers_saved": int(len(within)),
            "radius_arcsec": float(radius_arcsec),
            "nclosest": int(nclosest),
            "sframe_name": os.path.basename(sframe_path),
        }

        np.savez_compressed(
            outpath,
            wave=wave,
            flux=flux_sel,
            error=err_sel,
            fiberid=np.array(fiber_ids),
            ra_fiber_deg=np.array(ra_sel),
            dec_fiber_deg=np.array(dec_sel),
            sep_arcsec=np.array(sep_sel),
            meta=json.dumps(meta),
        )

        for j, fib in enumerate(fiber_ids):
            rows_summary.append(
                {
                    "Source_name": src,
                    "target_ra_deg": ra,
                    "target_dec_deg": dec,
                    "expnum": int(expnum),
                    "tileid": int(tileid),
                    "mjd": int(mjd),
                    "fiberid": int(fib) if str(fib).isdigit() else fib,
                    "sep_arcsec": float(sep_sel[j]),
                    "outfile": outname,
                }
            )

    return rows_summary


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matched",
        required=True,
        help="Matches file, e.g candidates_lvm.matched.txt",
    )
    parser.add_argument(
        "--outdir",
        default="extracted_spectra",
        help="output directory",
    )
    parser.add_argument(
        "--redux-version",
        default="1.2.1",
        help="SAS redux version, e.g. 1.2.1",
    )
    parser.add_argument(
        "--radius-arcsec",
        type=float,
        default=35.0,
        help="Maximum radius to store nearby fibers.",
    )
    parser.add_argument(
        "--nclosest",
        type=int,
        default=5,
        help="Maximum number of nearby fibers to save per target.",
    )
    parser.add_argument(
        "--keep-sframe",
        action="store_true",
        help="Maximum number of nearby fibers to save per target.",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading matches from: {args.matched}")
    matches = read_matched_table(args.matched)
    print(f"Valid matches: {len(matches)}")

    # Agrupar por exposición para no descargar repetidamente el mismo archivo
    grouped = matches.groupby("expnum", sort=True)

    all_rows = []

    for expnum, g in grouped:
        expnum = int(expnum)

        # Tomamos la primera fila para mjd/tileid
        mjd = int(g.iloc[0]["mjd"])
        tileid = int(g.iloc[0]["tileid"])

        try:
            url, fname = build_sframe_url(
                expnum=expnum,
                mjd=mjd,
                tileid=tileid,
                redux_version=args.redux_version,
            )
        except Exception as e:
            print(f"[{expnum}] Could not build URL: {e}")
            continue

        print(f"\n[{expnum}] Processing exposure")
        print(f"  URL: {url}")
        print(f"  Targets in this exposure: {len(g)}")

        tmpdir = tempfile.mkdtemp(prefix=f"lvm_{expnum:05d}_")
        local_path = os.path.join(tmpdir, fname)

        try:
            print("  Dowloading SFrame...")
            download_file(url, local_path)

            print("  Extracting spectrum...")
            rows = extract_target_spectra_from_exposure(
                sframe_path=local_path,
                expnum=expnum,
                tileid=tileid,
                mjd=mjd,
                targets_df=g,
                outdir=str(outdir),
                radius_arcsec=args.radius_arcsec,
                nclosest=args.nclosest,
            )
            all_rows.extend(rows)
            print(f"  Done. Saved spectra: {len(set([r['outfile'] for r in rows]))}")

        except requests.HTTPError as e:
            print(f"  HTTP error when downloading {url}: {e}")
        except Exception as e:
            print(f"  Error processing expnum={expnum}: {e}")
        finally:
            if args.keep_sframe:
                print(f"  Keeping temporal in: {local_path}")
            else:
                shutil.rmtree(tmpdir, ignore_errors=True)
                print("  Temporal file deleted.")

    # Guardar resumen
    if all_rows:
        summary = pd.DataFrame(all_rows)
        summary_path = outdir / "extracted_spectra_summary.csv"
        summary.to_csv(summary_path, index=False)
        print(f"\nSummary saved in: {summary_path}")
    else:
        print("\nNo spectra were saved..")


if __name__ == "__main__":
    main()