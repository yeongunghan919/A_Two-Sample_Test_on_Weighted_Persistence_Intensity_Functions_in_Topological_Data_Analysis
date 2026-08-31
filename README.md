# A Two-Sample Testing for Random Persistence Diagrams

This repository contains the code used in the paper

**“A Two-Sample Test for Random Persistence Diagrams Against Alternatives Characterized by Persistence Intensity Differences.”**

The proposed method performs two-sample testing for persistence
intensity functions using kernel-based test statistics and permutation
calibration. The repository also contains implementations of competing
methods and scripts for the simulation and real-data experiments presented
in the paper.

## Repository Structure

```text
.
├── source/                     
├── simulations/
│   ├── Circles/               
│   ├── Instrument/         
│   ├── orbit5k/              
│   └── Torus/
├── Per_Image_testing/
│   ├── PI_orbit5k_data/
│   ├── PI_circle_data/
│   └── PI_Torus_data/
├── results/
│   ├── figures/
│   │   └── figure_file/
│   └── simulation_results/
│       ├── Circle_results/
│       ├── Torus_results/
│       └── Orbit5k_results/
├── data/
│   ├── clarinet_diagrams/                
│   ├── flute_diagrams/          
│   ├── json/ 
│   ├── orbit_data/
│   └── Rdata
├── Figures.ipynb              
└── README.md
```

## Folder Description

### `source/`

This directory contains Python modules implementing the proposed bandwidth-aggregated two-sample test and other competing methods. 
The implementation of the aggregation procedure is adapted from the MMDAgg code of Schrab et al. (2023), with modifications for persistence intensity functions.


### `simulations/`

This directory contains scripts for the simulation studies.

- `Circles/` – Circles simulation experiment.
- `Instrument/` – Musical instrument experiment (flute vs. clarinet).
- `orbit5k/` – ORBIT5K simulation experiment.
- `Torus/` – Torus simulation experiment.

### `Per_Image_testing/`

This directory contains the implementation of the persistence image two-sample test.
The implementation of the two-stage two-sample testing method is adapted from the code of Moon and Lazar (2023).

### `results/`

This directory contains processed experimental outputs.

### `data/`

This directory contains source code for generating the data used in the paper rather than the generated data themselves

### `Figures.ipynb`

This is a jupyter notebook used to visualize the experimental results.