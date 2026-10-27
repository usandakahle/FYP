# SHARC: Directional wave measurement across the ice–ocean transition
 
Code and results for the EEE4022S final-year project by Sandakahle Bhengu, University of Cape Town, 2026.
 
The report quotes only the results that gave significant insight or led to a design change. All other results, and the tests that led to them, are in this repository, organised by the experiments E1–E7 of the report (Chapter 4).
 
> **Status:** the folders below are being populated.
 
---
 
## Where to find the results
 
| Results | Experiment | Hypothesis | Folder |
|---|---|---|---|
| Heave spectrum and H_m0 recovery on synthetic seas | E1 | H1 | `01_synthetic/` |
| Direction recovery on synthetic seas, estimator and fusion error | E2 | H1 | `01_synthetic/` |
| Swell direction with a second (wind sea) system | E3 | H1 | `01_synthetic/` |
| IMU noise, bias and six-position calibration | E4 | H1 | `02_imu_bench/` |
| Packet splitting, screening metrics and platform states | E5 | H3 | `03_screening/` |
| Anti-alias filter designs and alias margins | – | – | `04_filter_design/` |
| SB04 and SB06 directions against LP8, by platform phase | E6 | H2–H4 | `05_deployment/` |

 
---
 
## Repository layout
 
```
.
├── README.md
├── 00_pipeline/                 final processing pipeline
│   ├── run_sharc_jack.m         
│   ├── batchBias.m              whole deployment: heading proxy + direction tracking
│   └── build_yaw_checkpoint.m   heading proxy from GNSS drift bearing and yaw gyro
│
├── 01_synthetic/                E1–E3
│   ├── ocean_wave_spectrum.m
│   ├── synth_packet.m
│   ├── synth_verification.m
│   ├── initial_chain/           first iteration: sharc_meas.m, directional_validation.m
│   ├── results/
│   └── figures/
│
├── 02_imu_bench/                E4
│   ├── results/
│   └── figures/
│
├── 03_screening/                E5
│   ├── split_sb04_packets.m
│   ├── sharc_explore.m
│   ├── results/
│   └── figures/
│
├── 04_filter_design/
│   ├── sharc_alias_analysis.m
│   ├── results/                 alias_designs.csv, alias_summary_current.csv
│   └── figures/
│
├── 05_deployment/               E6
│   ├── SB04/
│   │   ├── results/
│   │   └── figures/
│   └── SB06/
│       ├── results/
│       └── figures/
│
│
└── data/                        input data, not stored on GitHub (see data/README.md)
```
 
In every folder, `results/` holds the tables (`.csv`) and saved workspaces (`.mat`), and `figures/` holds the figures produced by the scripts in that folder.
