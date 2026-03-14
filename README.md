
# Extract LVM Target Spectra

This repository contains a small Python tool to automatically extract spectra from **SDSS-V Local Volume Mapper (LVM)** observations for a list of targets.

The script takes as input the **cross-match table produced by Knox Long's `find_obs.py` routine**, identifies the corresponding LVM exposures, and extracts spectra from the **nearest fibers** to each target.

Only the relevant spectra are saved, avoiding the need to download or store entire LVM datasets.

---

# Overview

The script performs the following steps:

1. Reads the match file produced by the `drpall` cross-match (e.g. `candidatos_lvm.matched.txt`).
2. Groups targets by **exposure number (`expnum`)**.
3. Downloads the corresponding **`lvmSFrame`** file from the SDSS Science Archive Server (SAS).
4. Identifies the **nearest N fibers** to each target position.
5. Extracts the spectra from those fibers.
6. Saves the extracted spectra in **compressed `.npz` files**.
7. Deletes the temporary `lvmSFrame` file to minimize disk usage.

This allows efficient extraction of spectra for **large target lists** without downloading full LVM datasets.

---

# Requirements

Python ≥ 3.8

Required packages:

```
numpy
pandas
astropy
requests
```

Install them with:

```
pip install numpy pandas astropy requests
```

---

# Input

The script expects the **cross-match output produced by Knox Long's `find_obs.py` script**.

Example file:

```
candidatos_lvm.matched.txt
```

This file contains the list of targets and the LVM exposures that overlap their positions.

---

# Usage

Example command:

```
python extract_lvm_target_spectra.py   --matched candidatos_lvm.matched.txt   --outdir extracted_spectra   --redux-version 1.2.1   --radius-arcsec 15   --nclosest 2
```

---

# Parameters

### `--matched`

Cross-match file produced by `find_obs.py`.

Example:

```
candidatos_lvm.matched.txt
```

---

### `--outdir`

Directory where the extracted spectra will be saved.

Example:

```
extracted_spectra/
```

---

### `--redux-version`

LVM redux version stored on the SDSS SAS.

Example:

```
--redux-version 1.2.1
```

---

### `--radius-arcsec`

Maximum radius (arcseconds) used to select nearby fibers.

Example:

```
--radius-arcsec 15
```

---

### `--nclosest`

Maximum number of closest fibers to extract per target.

Example:

```
--nclosest 2
```

---

# Output

The script produces:

### Extracted spectra

Compressed `.npz` files containing:

- wavelength array
- flux
- flux error
- fiber IDs
- fiber coordinates
- angular separation from the target

Example file:

```
cand_001_exp15294.npz
```

---

### Summary table

A CSV file listing all extracted fibers:

```
extracted_spectra_summary.csv
```

Columns include:

- target name
- target coordinates
- exposure number
- tile ID
- fiber ID
- angular separation
- output file name

---

# Data Access

The script downloads **`lvmSFrame` files directly from the SDSS Science Archive Server (SAS)**.

Access requires valid **SDSS credentials**.

Credentials must be provided in the script using:

```
USER = "your_sdss_username"
PASSWORD = "your_sdss_password"
```

---

# Notes

- The script downloads only the necessary exposures and deletes them after processing to reduce disk usage.
- Multiple targets falling within the same exposure are processed simultaneously.
- The number of fibers extracted per target can be controlled using `--nclosest`.

---

# Acknowledgments

This tool relies on the cross-matching routine developed by **Knox Long** (`find_obs.py`) to identify LVM exposures overlapping the target positions.

LVM data are part of the **SDSS-V Local Volume Mapper survey**.

---

# License

MIT License (or specify your preferred license).
