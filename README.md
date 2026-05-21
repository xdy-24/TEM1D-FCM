# TEM1D-FCM

MATLAB implementation for one-dimensional transient electromagnetic inversion with fuzzy C-means clustering constraints.

This repository provides the source code associated with the manuscript:

**Weak Structural Constraints Enhance the Performance of Multi-Parameter Transient Electromagnetic Inversion in Metal Mines**

submitted to *Computers & Geosciences*.

## Description

TEM1D-FCM is a MATLAB implementation for one-dimensional transient electromagnetic inversion considering induced polarization effects. The inversion uses the Cole-Cole model and incorporates fuzzy C-means clustering constraints to improve the recovery of model parameters.

The example provided in this repository uses the `realmodel.txt` file as input data.

## Requirements

- MATLAB R2024a or later is recommended.
- No additional MATLAB toolboxes are required unless otherwise stated in the code comments.

## Repository contents

```text
TEM1D-FCM/
├── README.md
├── LICENSE
├── run_quick_test.m
├── inv_main.m
├── occam1D_inversion_with_clustering.m
├── calculate_jacobian_with_clustering.m
├── calculate_misfit_with_clustering.m
├── classify_layers.m
├── fcm_cluster.m
├── guided_fcm.m
└── realmodel.txt
```

## Main files

- `run_quick_test.m`: Quick-test script for running the example.
- `inv_main.m`: Main inversion program.
- `occam1D_inversion_with_clustering.m`: Main inversion subroutine with clustering constraints.
- `calculate_jacobian_with_clustering.m`: Jacobian calculation with clustering constraints.
- `calculate_misfit_with_clustering.m`: Misfit calculation with clustering constraints.
- `classify_layers.m`: Layer classification function.
- `fcm_cluster.m`: Fuzzy C-means clustering function.
- `guided_fcm.m`: Guided fuzzy C-means clustering function.
- `realmodel.txt`: Example input data file containing time gates and corresponding data.

## Input data

The input file `realmodel.txt` contains two columns:

1. Time gates
2. Corresponding transient electromagnetic data

The current `realmodel.txt` file provides an example dataset for testing the code.

## Quick test

To run a quick test, download or clone this repository, open MATLAB, set the repository folder as the current working directory, and run:

```matlab
run_quick_test
```

Alternatively, run:

```matlab
addpath(genpath(pwd));
inv_main
```

## Expected output

After running the quick test, the program performs the example inversion and generates output files such as:

- `inversion_log.txt`
- `layer_classification.txt`
- `parameter_log.txt`
- `residuals_log.txt`

The program also generates figures showing the inversion results, including resistivity and chargeability fitting curves, decay voltage fitting curve, iteration curve, clustering results, and time-gate residuals.

## Notes

The inversion uses guided fuzzy C-means clustering. Because random initialization is used in the clustering procedure, the inversion results may show slight variations between different runs.

For a new dataset, users should modify the input data file and the corresponding model and inversion parameters in `inv_main.m` and `occam1D_inversion_with_clustering.m`.

## License

This project is released under the MIT License. See the `LICENSE` file for details.

## Contact

Dongyang Xie  
Email: 2024126079@chd.edu.cn

