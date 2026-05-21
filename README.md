# TEM1D-FCM

MATLAB implementation for one-dimensional transient electromagnetic inversion with fuzzy C-means clustering constraints.

This repository provides the computer code associated with the manuscript:

**Weak Structural Constraints Enhance the Performance of Multi-Parameter Transient Electromagnetic Inversion in Metal Mines**

submitted to *Computers & Geosciences*.

## Description

TEM1D-FCM is a MATLAB implementation for one-dimensional transient electromagnetic inversion considering induced polarization effects. The inversion uses the Cole-Cole model and incorporates fuzzy C-means clustering constraints to improve the recovery of model parameters.

The example provided in this repository uses the `realmodel.txt` file as input data.

## Requirements

- MATLAB R2024a or later is recommended.
- A Windows 64-bit operating system is required for the provided `mexBtFwdArbitraryLoop.mexw64` file.
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
├── mexBtFwdArbitraryLoop.mexw64
└── realmodel.txt
```

## Main files

- `run_quick_test.m`: Quick-test script for running the example.
- `inv_main.m`: Main inversion program.
- `mexBtFwdArbitraryLoop.mexw64`: Compiled MATLAB MEX forward-modeling program used to compute the transient electromagnetic response. This file is intended for Windows 64-bit MATLAB environments.
- `occam1D_inversion_with_clustering.m`: Main inversion subroutine with clustering constraints.
- `calculate_jacobian_with_clustering.m`: Jacobian calculation with clustering constraints.
- `calculate_misfit_with_clustering.m`: Misfit calculation with clustering constraints.
- `classify_layers.m`: Layer classification function.
- `fcm_cluster.m`: Fuzzy C-means clustering function.
- `guided_fcm.m`: Guided fuzzy C-means clustering function.
- `realmodel.txt`: Example input data file containing time gates and corresponding transient electromagnetic data.

## Input data

The input file `realmodel.txt` contains two columns:

1. Time gates
2. Corresponding transient electromagnetic data

The current `realmodel.txt` file provides an example dataset for testing the code.

## Quick test

To run the example, users should download or clone the entire repository rather than downloading only `run_quick_test.m`. The quick-test script depends on the MATLAB subroutines, the example input data file `realmodel.txt`, and the compiled MEX forward-modeling routine `mexBtFwdArbitraryLoop.mexw64`.

Users can download the repository by clicking **Code > Download ZIP** on the GitHub page, or by using:

```bash
git clone https://github.com/xdy-24/TEM1D-FCM.git
```

After downloading the repository, open MATLAB, set the repository folder as the current working directory, and run:

```matlab
run_quick_test
```

Alternatively, users can run:

```matlab
addpath(genpath(pwd));
inv_main
```

## Expected output

After running the quick test, the program performs the example inversion and generates output text files, including:

- `inversion_log.txt`
- `layer_classification.txt`
- `parameter_log.txt`
- `residuals_log.txt`

In addition, the program generates five MATLAB figures showing:

1. Resistivity and chargeability fitting curves
2. Time-gate errors
3. Decay voltage fitting curve
4. Iteration curve
5. Clustering results

## Notes

Because the guided fuzzy C-means clustering procedure involves random initialization, repeated runs may produce slightly different numerical results. Therefore, the outputs may not be exactly identical to the example results, although the overall trends should remain consistent.

For a new dataset, users should modify the input data file and the corresponding model and inversion parameters in `inv_main.m` and `occam1D_inversion_with_clustering.m`.

## Note on the MEX forward-modeling program

The file `mexBtFwdArbitraryLoop.mexw64` is a compiled MATLAB MEX forward-modeling program used by the inversion code to calculate the transient electromagnetic response. The provided MEX file is intended for Windows 64-bit MATLAB environments. The MATLAB scripts in this repository implement the inversion framework and the fuzzy C-means clustering constraints used in the manuscript.

## License

This project is released under the MIT License. See the `LICENSE` file for details.

## Contact

Developer: Dongyang Xie  
Maintainer: Zhipeng Qi  
Contact: Zhipeng Qi, qzhipeng@126.com
