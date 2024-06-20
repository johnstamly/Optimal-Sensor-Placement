# Optimal Sensor Placement


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12180851.svg)](https://doi.org/10.5281/zenodo.12180851)



This repository contains a Python script for optimal dynamic sensor placement using the Fisher Spectral Radius criterion. The code was inspired by the paper of Loutas T.H. and Bourikas A [[1]](#1).

## Citation

If you use this software, please cite it as:

```bibtex
@misc{STAMATELATOS2024ODSP,
  author       = {Giannis Stamatelatos},
  title        = {Optimal Dynamic Sensor Placement},
  year         = {2024},
  month        = {June},
  note         = {Release 1.0.0, 20th June 2024},
  url          = {https://github.com/johnstamly/Optimal-Sensor-Placement},
  doi          = {10.5281/zenodo.12180591}
}
```


## Description

The script performs optimal sensor placement in a Finite Element Model (FEM) based on the Fisher Spectral Radius criterion. It aims to find the minimum number of sensors and their optimal locations to maximize information and guarantee the capture of all varieties of modal information by introducing redundancy and convergence criteria.

## Features

- Computes the Fisher Spectral Radius for each node in the FEM model
- Finds the maximum Fisher node and updates the Fisher Information matrix
- Iteratively selects the optimal sensor locations based on the spectral difference ratio
- Provides the percentage of the final SpRI compared to the initial model SpRI
- Calculates the convergence of the last placed sensor
- Displays the results in a formatted table

## Requirements

- Python 3.x
- NumPy
- h5py
- tqdm
- PrettyTable

## Installation

1. Clone the repository:
```
git clone https://github.com/johnstamly/optimal-sensor-placement.git
```

2. Install the required dependencies:
```
pip install numpy h5py tqdm prettytable
```

## Usage

1. Prepare the eigenvector matrices (e.g. `T1.mat` and `T2.mat`) and place them in the same directory as the script.
   (T1.mat and T2.mat are provided as an example. They contain the 10 first strain eigen modes of x and y axis respectively for each d.o.f.
   You can use also z-axis by changing the "dim" parameter)
   
3. Run the script:
```
 python optimal_sensor_placement.py
```

3. Enter the desired convergence error and initial sensor count when prompted.

4. The script will perform the optimal sensor placement and display the results in a formatted table.

## Results

The script will output the following results:

- Sensor positions (node indices)
- Number of sensors placed
- Percentage of the final SpRI compared to the initial model SpRI
- Convergence of the last placed sensor
- Execution time in seconds
- Number of loops performed
- Threshold percentage (Ts%)
- Final SpRI value

## Authors

- Giannis Stamatelatos
- Theodoros Loutas

## Affiliation

University of Patras, Department of Mechanical Engineering & Aeronautics,
Laboratory of Applied Mechanics and Vibrations

## References
<a id="1">[1]</a> 
Loutas, T. H., & Bourikas, A. (2017). Strain sensors optimal placement for vibration-based structural health monitoring. The effect of damage on the initially optimal configuration. Journal of Sound and Vibration, 410, 217-230. https://doi.org/10.1016/j.jsv.2017.08.033

## License

This project is licensed under the [MIT License](LICENSE).






