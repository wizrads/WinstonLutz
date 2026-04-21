# MP‑566 Lab 8 — Winston‑Lutz Test

Materials for **MP‑566: Physics of Radiotherapy, Lab 8** at UW‑Madison. The lab walks students through acquiring Winston‑Lutz (WL) MV portal images on a Varian TrueBeam, then analyzing them to check whether the **mechanical** and **radiation** isocenters coincide within AAPM TG‑142 / TG‑198 tolerance for SRS/SBRT machines.

Both the original MATLAB analysis and a modern Python/Colab version are included. The Python version uses [**pylinac**](https://pylinac.readthedocs.io/en/latest/winston_lutz.html), which implements the same Low *et al.* (1995) formalism as the MATLAB code.

---

## Repository contents

| File | Purpose |
|---|---|
| `Lab8_WinstonLutz.ipynb` | **Student notebook** (Google Colab–ready). Walks through DICOM familiarization, pylinac analysis, and pass/fail against TG‑142. Contains four `TODO` cells. |
| `Analyze_Fields.m` | Original MATLAB analysis (hand‑written Low *et al.* formalism with `imfindcircles` + `pinv`). |
| `WL_plotter.m` | MATLAB helper to overlay the radiation field and BB on an image. |
| `Coding_Workspace.m` | MATLAB entry point / scratchpad. |
| `Lab_Data_Final.zip` | 14 `.dcm` portal images (one per WL field in Table 1 of the handout). |

The 14 portal images correspond to the gantry/collimator/couch combinations in **Table 1** of the handout:

| Fields | Gantry | Collimator | Couch |
|---|---|---|---|
| 1–3 | 270° | 90 / 0 / 270 | 0 |
| 4–6 | 0° | 90 / 0 / 270 | 0 |
| 7–9 | 90° | 90 / 0 / 270 | 0 |
| 10–12 | 180° | 90 / 0 / 270 | 0 |
| 13–14 | 0° | 0 | 90 / 270 |

---

## Quick start — Google Colab (recommended for students)

1. Download `Lab8_WinstonLutz.ipynb` and the 14 `.dcm` files from `Lab_Data_Final/` (unzip `Lab_Data_Final.zip` locally first).
2. Open [Google Colab](https://colab.research.google.com/) → **Upload notebook** → select `Lab8_WinstonLutz.ipynb`.
3. Run the cells top‑to‑bottom. The upload cell will prompt you to pick the 14 `.dcm` files; drag‑select them all in the dialog.
4. Fill in the four `TODO` cells as you go:
   - **TODO #1** — the BB diameter in mm.
   - **TODO #2** — the `Δu` / `Δv` / total per image.
   - **TODO #3** — the magnitude of the 3D shift vector.
   - **TODO #4** — the TG‑142/TG‑198 tolerance and pass/fail check.

## Quick start — local Python

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python analyze_fields.py            # end‑to‑end analysis, no TODOs
# or
jupyter lab Lab8_WinstonLutz.ipynb  # work through the student version
```

## Quick start — MATLAB (original reference)

```matlab
[images, shift_3D] = Analyze_Fields('Lab_Data_Final');
WL_plotter(images(1))       % overlay field CAX + BB on field #1
```

---

## What the analysis produces

Running either the MATLAB code or the Python notebook on the provided data yields (within segmentation precision):

**Per‑image 2D shifts** (Δu, Δv, total — 14 rows, one per field).

**3D shift vector** from the radiation isocenter to the BB center:

---

## Method

For each image, locate the **radiation field centroid** and the **BB centroid** on the EPID, convert pixel coordinates to mm at the isocenter plane (multiply by `ImagePlanePixelSpacing / (RTImageSID/RadiationMachineSAD)`), and record the in‑plane shift `(Δu, Δv)`.

Stack the 14 in‑plane shifts with the first two rows of each field's `R_gantry · R_couch` rotation matrix and solve, by least‑squares (Moore–Penrose pseudoinverse), for the single 3D offset between the radiation isocenter and the BB. This is the formalism from Low *et al.* (1995) and is what both the MATLAB script and pylinac implement.

---

## References

1. Lutz, W., Winston, K. R., & Maleki, N. (1988). *A system for stereotactic radiosurgery with a linear accelerator.*
2. Low, D. A., Li, Z., & Drzymala, R. E. (1995). *Minimization of target positioning error in accelerator-based radiosurgery.* Med. Phys. 22(4): 443–448.
3. AAPM TG‑142. *Quality assurance of medical accelerators.*
4. AAPM TG‑198. *An Implementation Guide for TG‑142 Quality Assurance of Medical Accelerators.*
5. AAPM MPPG Report 9.a. *Guidelines for SRS‑SBRT.*
6. Pylinac documentation: <https://pylinac.readthedocs.io/en/latest/winston_lutz.html>

---

## License / attribution

Lab materials © the UW‑Madison MP‑566 course staff. Code and notebook provided for educational use. Pylinac is licensed under the MIT license.
