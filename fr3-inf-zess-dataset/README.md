# Wideband Channel Measurements in the Lower FR3 Band (6–12 GHz) — ZESS Production Hall, RUB

Wideband VNA channel measurements in an industrial production hall at Ruhr University Bochum,
with boresight-aligned horn antennas, for line-of-sight (LOS) and non-line-of-sight (NLOS, blockage
by an industrial robot) conditions. This repository contains the raw S-parameters, the receiver
positions, the antenna gain data and the processing code used in:

> F. Rauf, A. Umra, O. Alsoloh and A. Sezgin, "Wideband Channel Measurements in the Lower FR3 Band
> (6–12 GHz) for Industrial Environments," submitted to IEEE ICC 2027.

## Measurement setup

| Parameter | Value |
|---|---|
| Environment | ZESS production hall, Ruhr University Bochum, 32.73 m × 23.29 m, ceiling height 10 m |
| Equipment | Rohde & Schwarz ZNA vector network analyser |
| Measured range | 6–17 GHz, 11 001 points (1 MHz step) |
| IF bandwidth / power / averaging | 1 kHz / 10 dBm / 10 sweeps (vector) |
| Calibration | 2-port UOSM, 2.92 mm kit, reference plane at the antenna connectors |
| Antennas | Schwarzbeck BBHA 9120 L double-ridged horns (Tx and Rx), vertical polarisation |
| Alignment | boresight to boresight at every position |
| Tx / Rx height | 1.5 m / 1.5 m |
| LOS campaign | 135 Rx positions on a 1-m grid, 1–18 m (Tx 1); positions behind the transparent robot-cell fence are treated as LOS |
| NLOS campaign | 30 Rx positions behind an industrial robot, 2.5–16.5 m (Tx 2) |

## Repository structure

```
data/
  LOS/            pt1.s2p ... pt135.s2p      raw S-parameters, LOS campaign (Touchstone, 6-17 GHz)
                  positions.csv              pt, x_m, y_m, distance_m, los, tx
  NLOS/           pt1.s2p ... pt30.s2p       raw S-parameters, NLOS campaign
                  positions_nlos.csv         pt, x_m, y_m, distance_m, los, tx
  calibration/    thru / reference checks    through and 1-m reference measurements
  antenna/        gain.csv                   BBHA 9120 L gain (freq_GHz, gain_dBi), datasheet values
  band_validity.csv                          validity of each 1-GHz sub-band (see below)
code/
  fr3_inf_analysis.m                         full processing: path loss, CI/FI fits, delay spread, figures, tables
  export_plot_csv.m                          plot-ready CSV export
  check_1m.m                                 in-lab check of a 1-m reference measurement against free space
results/                                     fitted parameters and per-point results as CSV
```

LOS and NLOS point numbers are independent (LOS `pt1` and NLOS `pt1` are different positions).
The coordinates in the two position files refer to the hall frame described in the paper.

## Data format

* `.s2p`: Touchstone 2-port files written by the VNA (`# HZ S RI R 50`), columns
  `freq, re S11, im S11, re S21, im S21, re S12, im S12, re S22, im S22`.
  The channel is in S21. The antenna gains are **not** removed (see `data/antenna/gain.csv`).
* `positions*.csv`: `pt` point number, `x_m`, `y_m` hall coordinates, `distance_m` Tx–Rx distance,
  `los` (1 = LOS, 0 = NLOS), `tx` transmitter position (1 = LOS campaign, 2 = NLOS campaign).
  The hall origin is the bottom-left corner of the floor plan in Fig. 1 of the paper; Tx 1 is at (26.737, 11.147) m.

## Sub-band validity — please read

Each 1-GHz sub-band was checked against free-space propagation at a 1-m reference position.
**Only the sub-bands from 6 to 12 GHz passed this check and are used in the paper.**
Above 12 GHz, the measurements are affected by a calibration/cable error and should not be used
for propagation analysis. See `data/band_validity.csv`.

## Processing

Run `code/fr3_inf_analysis.m` (MATLAB, no toolboxes; also runs in GNU Octave). All paths are relative
to the repository, so no configuration is needed after cloning:

```matlab
cd code
fr3_inf_analysis
```

It writes the figures, the LaTeX parameter table and CSV files to `results/`.
`check_1m.m` repeats the 1-m free-space check used to select the valid sub-bands.


## Citation

If you use this dataset, please cite the paper above (see `CITATION.cff`).

## Acknowledgement

This work was supported in part by the German Federal Ministry of Research, Technology and Space
(BMFTR) in the course of the 6GEM+ Transfer Hub under grant 16KIS2411.

## Contact

Faizan Rauf, Chair of Digital Communication Systems, Ruhr University Bochum — faizan.rauf@rub.de
