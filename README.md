# falklands-gentoos

Code and analyses accompanying the manuscript:

> **Multi-colony tracking of a marine central place forager reveals a site specific yet broadly consistent foraging strategy**

## Overview

This repository contains the R code used to analyse movement, diving behaviour, and habitat use of gentoo penguins breeding at multiple colonies in the Falkland Islands. The analyses investigate colony-specific variation and broader patterns in foraging strategies across colonies.

The repository includes scripts for:

* Preparing environmental raster datasets
* Quantifying seafloor depth associated with penguin locations
* Fitting hurdle models of foraging behaviour
* Fitting linear mixed models for trip, dive, and benthic foraging metrics
* Generating figures and outputs used in the manuscript

## Manuscript

**Article:** *Multi-colony tracking of a marine central place forager reveals a site specific yet broadly consistent foraging strategy*

**Publication link:** *[Placeholder – add DOI or journal URL upon publication]*

## Repository Structure

```text
.
├── figures/                               # Figures and graphical outputs
├── hurdle_models.R                        # Hurdle models of foraging metrics
├── linear_mixed_models_dives.R            # Mixed models for dive characteristics
├── linear_mixed_models_prop_benthic.R     # Mixed models for benthic diving proportion
├── linear_mixed_models_trips.R            # Mixed models for trip characteristics
├── prep_rasters.R                         # Preparation of environmental raster layers
├── quantify_seafloor_depth.R              # Seafloor depth extraction and processing
└── README.md
```

## Requirements

Analyses were conducted in **R**.
Required packages are listed in the respective scripts.

## Data Availability

This repository contains the code used for the analyses. Tracking data are available at *[Placeholder – add seabird tracking URL]*

## Citation

If you use this code, please cite:

**Handley, J.M. et al. (2026)** *Multi-colony tracking of a marine central place forager reveals a site specific yet broadly consistent foraging strategy.*, [add JOURNAL NAME, ISSUE, AND PAGES]

DOI: *[Placeholder]*

## Authors

**Jonathan Handley**

**Joshua Wilson**

For questions regarding the analyses or code, please contact Joshua Wilson directly.

