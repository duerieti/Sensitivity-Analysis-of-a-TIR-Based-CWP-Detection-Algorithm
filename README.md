# Sensitivity Analysis of a TIR-Based Cold-Water Patch Detection Algorithm

**Author:** Etienne Dürig  
**Institution:** ZHAW Institute of Computational Life Sciences (ICLS)  
**Date:** July 2026

---

## 1. Project Overview

TODO

---
> **Note:** This project was developed and run on a Debian-based Linux system
> and is best suited for execution on a distributed HPC cluster running a
> Debian-like Linux OS with SLURM as the job scheduler and Conda as the
> environment manager. The parallel execution scripts rely on GNU Parallel,
> which is a Linux/macOS tool and is not natively available on Windows (but can be deployed using Conda). Windows
> users who wish to run the full pipeline locally are advised to install WSL 2
> with Ubuntu:  
> https://learn.microsoft.com/en-us/windows/wsl/install
>
> When running on Windows, file paths and line endings may need to be adjusted.
> Any hardcoded paths using `/` may need to be updated, and shell scripts should

---
## 2. Prerequisites

### 2.1 ZHAW HPC Account (recommended)

If you do not already have an account on the ZHAW HPC, submit a request at:  
https://wiki.hpc.zhaw.ch/hpcuserwiki/index.php/Access:New%3F_Get_an_account  
*(Only accessible within the ZHAW network.)*

After account creation, complete the first login and connect to the cluster:
- First login: https://wiki.hpc.zhaw.ch/hpcuserwiki/index.php/Getting_started:First_Login  
- Connect: https://wiki.hpc.zhaw.ch/hpcuserwiki/index.php/Getting_started:Connect_to_the_cluster

The cluster runs Linux. The following terminal commands will be useful for navigation:

| Command | Purpose |
|---------|---------|
| `cd` | Change directory |
| `ls` | List directory contents |
| `mkdir` | Create a directory |
| `nano` | View and edit files in the terminal |

### 2.2 Conda

If not already installed, download Conda from the official documentation:  
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

### 2.3 Git

**Windows**  
https://git-scm.com/install/windows

**macOS** (requires Homebrew)
```bash
brew install git
```

**Linux (Debian-based)**
```bash
sudo apt install git
```

### 2.4 GNU Parallel

GNU Parallel is required for the parallel execution scripts.

**Linux (Debian-based)**
```bash
sudo apt install parallel
```

**macOS** (requires Homebrew)
```bash
brew install parallel
```

**Windows (via Conda)**

You might add parallel to one of the environments.

```bash
conda activate r_env
conda install -c conda-forge parallel
```

---

## 3. Installation

### 3.1 Clone the Repository

```bash
git clone https://github.com/duerieti/Sensitivity-Analysis-of-a-TIR-Based-CWP-Detection-Algorithm.git
```

### 3.2 Set Up Environments

#### Option A: Using Conda (recommended)

Two separate environments are provided. Note that the environment files were
exported on Linux and may require minor adjustments on other platforms.

**Sensitivity analysis environment**
```bash
conda env create -f sensitivity_env.yml
```

**Algorithm evaluation environment**
```bash
conda env create -f r_env.yml
```

#### Option B: Using R directly

Install the required packages in R or RStudio with `install.packages()`:

```r
install.packages(c(
  "tidyverse",
  "sf",
  "terra",
  "exactextractr",
  "tmap",
  "maptiles",
  "tmaptools",
  "rcpp",
  "data.table",
  "ggplot2",
  "future",
  "furrr",
  "tictoc",
  "logger",
  "stars",
  "raster",
  "lwgeom",
  "leaflet",
  "leafem",
  "sensitivity",
  "stringdist",
  "abind"
))
```
**R version:** 4.5.3

**System dependencies** (required for geospatial packages, Linux/macOS only):
- GDAL ≥ 3.0
- GEOS ≥ 3.8
- PROJ ≥ 6.0


---

## 4. Repository Structure

```
├── comparison/                            # Compare v1, v2, and v3 algorithm versions
│   ├── compute_jaccard.R                  # Compute Jaccard similarities between algorithms
│   ├── param_combi_1/                     # First parameter combination
│   │   ├── current/                       # v1 algorithm output
│   │   ├── new/                           # v2 algorithm output
│   │   └── v3_12_core/                    # v3 algorithm output
│   ├── param_combi_2/                     # Second parameter combination
│   │   ├── current/                       # v1 algorithm output
│   │   └── new/                           # v2 algorithm output
│   └── param_combi_3/                     # Third parameter combination
│       ├── current/                       # v1 algorithm output
│       └── new/                           # v2 algorithm output
│
├── comparison_for_speed/                  # Compare v1 and v2 for runtime benchmarking
│   ├── param_combi_1/                     # First parameter combination
│   │   ├── current/                       # v1 algorithm output
│   │   └── new/                           # v2 algorithm output
│   ├── param_combi_2/                     # Second parameter combination
│   │   ├── current/                       # v1 algorithm output
│   │   └── new/                           # v2 algorithm output
│   └── param_combi_3/                     # Third parameter combination
│       ├── current/                       # v1 algorithm output
│       └── new/                           # v2 algorithm output
│
├── compute_morris/                        # Scripts and data for Morris sensitivity analysis
│   ├── compute_aggregated_parameters.R    # Compute CWP area and temperature delta per Morris run
│   ├── compute_detection_frequency_raster.R  # Compute detection frequency raster
│   ├── compute_morris_indices.R           # Compute Morris indices from CWP area
│   ├── cwp_annotated.csv                  # Annotated CWP locations with classifications
│   ├── lumped_stats_emme_big.csv          # Aggregated statistics per CWP and Morris run
│   ├── sa_object_big.rds                  # Second Morris sampling grid object
│   └── data/                             # Input data
│       ├── derived_data_products/         # Pre-rounded rasters for faster computation
│       └── original_data/                # TIR data from Tonolla and Antonetti
│           ├── Centerlines_FINAL/         # River centerlines
│           └── thermal_rasters_FINAL/    # TIR orthophoto rasters
│
└── functions/                             # Algorithm implementations and utility functions
    ├── Cold_water_patch_detection_anto.R  # Original v1 implementation (Tonolla & Antonetti)
    ├── CWP_code_commented_claude.R        # Commented version of the original v1 implementation
    ├── detect_cwp_v2_internal_rounding.R  # v2 algorithm with internal TIR raster rounding
    ├── detect_cwp_v2.R                    # v2 algorithm without internal rounding
    ├── detect_cwp_v3.R                    # v3 algorithm (IDW-based reference temperature)
    ├── function_current.R                 # v1 algorithm isolated as a standalone function
    ├── idw.cpp                            # IDW function in C++ for v3 (imported via Rcpp)
    ├── jaccard.R                          # Jaccard similarity computation
    ├── sample_parameters_big.R            # Second Morris parameter grid generation
    └── sample_parameters.R               # First Morris parameter grid generation
```
