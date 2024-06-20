# Optimal Sensor Placement

This repository contains a Python script for optimal sensor placement for dynamic messurments (strain, accelerations) using the Fisher Spectral Radius criterion.

## Description

The script performs optimal sensor placement in a Finite Element Model (FEM) based on the Fisher Spectral Radius criterion. It aims to find the optimal locations for placing sensors to maximize the information gained from the system.

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

git clone https://github.com/johnstamly/optimal-sensor-placement.git


2. Install the required dependencies:

pip install numpy h5py tqdm prettytable


## Usage

1. Prepare the eigenvector matrices `T1.mat` and `T2.mat` and place them in the same directory as the script.

2. Run the script:

python optimal_sensor_placement.py


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

## License

This project is licensed under the [MIT License](LICENSE).
