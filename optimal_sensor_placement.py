import h5py
import numpy as np
import time
import os
import multiprocessing
from tqdm import tqdm 
from prettytable import PrettyTable

table = PrettyTable()
table.field_names = ["Description"]

table.add_row(["OPTIMAL SENSOR PLACEMENT"])
table.add_row([""])
table.add_row(["University of Patras, Department of Mechanical Engineering & Aeronautics"])
table.add_row([""])
table.add_row(["Authors: Giannis Stamatelatos, Theodoros Loutas"])
table.add_row([""])

print(table)

num_threads = multiprocessing.cpu_count()
os.environ['NUMEXPR_MAX_THREADS'] = str(num_threads)
os.environ['NUMEXPR_NUM_THREADS'] = str(num_threads)
print(f"\nRunning....\nusing cores/threads: {num_threads}")

########### INPUTS ##################################################
########### INPUT PARAMETERS ###########
converg =  float(input("\nEnter the desired convergence error (e.g. 0.05): "))   # Convergence threshold for the last selected sensor
Ts_initial = 0.0  # Threshold for similar information Fisher matrices
Ns_maximum = 30   # Maximum Number of sensors to search
Ns_initial = int(input("\nEnter initial sensor count: "))  # Initial Number of sensors to search
Dim= 2            # Dimension (e.g., 2D or 3D) is the number of DOFs per node



################# FUNCTIONS #########################################

def fisher_spectral(E, I_init, Dim):
    """
    Compute the Fisher Spectral Radius for each node in the FEM model.

    Args:
        E (numpy.ndarray): Eigenvector matrix.
        I_init (numpy.ndarray): Initial Fisher Information matrix.
        Dim (int): Dimension of the problem (e.g., 2D or 3D).

    Returns:
        tuple: A tuple containing two numpy arrays:
            - Ikx0: Updated Fisher Information matrix for each node.
            - SpRx: Spectral radius for each node.
    """
    num_nodes = int(np.size(E, 0) / Dim)  # Number of nodes in the FEM model
    Ikx0 = []
    SpRx = []

    for k in range(num_nodes):
        # Extract eigenvectors for the current node
        N = E[Dim*k : Dim*k + Dim, :]

        # Compute the updated Fisher Information matrix
        Ikx = np.dot(N.T, N) 
        Ikx0.append(Ikx)

        # Compute the spectral radius (maximum eigenvalue)
        w = np.linalg.eigvals(Ikx + I_init)
        w = np.abs(w)
        SpRx.append(max(w))

    Ikx0 = np.array(Ikx0)
    SpRx = np.array(SpRx)

    return Ikx0, SpRx


def spectral_dif(E, Ikmax, k, Dim):
    """
    Compute the spectral difference ratio for a given node in the FEM model.

    Args:
        E (numpy.ndarray): Eigenvector matrix.
        Ikmax (numpy.ndarray): Maximum Fisher Information matrix.
        k (int): Index of the node.
        Dim (int): Dimension of the problem (e.g., 2D or 3D).

    Returns:
        float: Spectral difference ratio (Rk) for the given node.
    """
    # Extract eigenvectors for the current node
    N = E[Dim*k : Dim*k + Dim, :]

    # Compute the Fisher Information matrix for the current node
    Ik = np.dot(N.T, N)

    # Compute the spectral difference ratio
    Temp1 = Ikmax - Ik
    Temp2 = Ikmax + Ik
    w1 = np.linalg.eigvals(Temp1)
    w1 = max(np.abs(w1))
    w2 = np.linalg.eigvals(Temp2)
    w2 = max(np.abs(w2))
    Rk = w1 / w2

    return Rk


def sensor_placement(Ts, Ns, Dim):
    """
    Perform sensor placement based on the Fisher Spectral Radius criterion.

    Args:
        Ts (float): Threshold for spectral difference ratio.
        Ns (int): Number of sensors to place.
        Dim (int): Dimension of the problem (e.g., 2D or 3D).

    Returns:
        tuple: A tuple containing the following:
            - Sensor_Node_Position (list): List of node indices where sensors are placed.
            - Prcnt_SpRI (float): Percentage of the final SpRI compared to the initial model SpRI.
            - div (float): Convergence of the last sensor placed.
            - SpRI (float): Final SpRI after sensor placement.
    """
    # Loading Eigenvector Matrices
    T1 = h5py.File("T1.mat")
    T2 = h5py.File("T2.mat")  # T1 and T2 are the corresponding x and y eigenvectors

    data = {}  # Creating dictionary to save data values
    for k, v in T1.items():
        data[k] = np.array(v)
    for k, v in T2.items():
        data[k] = np.array(v)

    data["T1"] = data['T1'].T  # Transpose matrices to shape: (DOF x eigenvectors)
    data["T2"] = data['T2'].T

    E = []
    for i in range(data['T1'].shape[0]):  # Combine the T1 and T2 matrices row by row
        row = data['T1'][i, :]
        E.append(row)
        row = data['T2'][i, :]
        E.append(row)
    E = np.asarray(E)

    # Parameters Initialization
    Scx = E.shape[0]  # Total available DOFs
    I = np.zeros((np.size(E, 1), np.size(E, 1)))  # Initialization of I (10x10) zeros
    SpRI = max(np.abs(np.linalg.eigvals(I)))
    Sensor_Node_Position = []  # List of Sensor position (as nodes)
    Ikx0, SpRx = fisher_spectral(E, I, Dim)  # Calculating the model total SpRI
    SpRI_model = sum(SpRx)

    # Find maximum Fisher node through corresponding maximum Spectral Radius and zeroing the eigenvector matrix for each sensor
    # and zeroing also the eigenvectors row of each node that its Rx is lesser than Ts.
    while Ns > 0 and Scx > 0:
        Ikx0, SpRx = fisher_spectral(E, I, Dim)
        idx = np.argmax(SpRx)
        Ikxmax = Ikx0[idx]
        I = I + Ikxmax
        Sensor_Node_Position.append(idx)
        E[idx:idx+Dim, :] = 0  # Zeroing the eigenvalues of the sensor node
        Scx = Scx - 1  # Subtract the corresponding DOF for break condition
        PIN = []  # List of indexes of similar information matrices below Ts threshold
        for k in range(int(E.shape[0] / Dim)):
            if np.any(E[Dim*k:Dim*k + Dim, :] != 0):
                Rk = spectral_dif(E, Ikxmax, k, Dim)
                if Rk < Ts:
                    for i in range(Dim):
                        PIN.append(Dim*k + i)

        E[PIN, :] = 0  # Zeroing all rows with similar info
        Scx = Scx - len(PIN)  # Subtract the corresponding DOFs for break condition
        Ns = Ns - 1  # Subtract from the pool of remaining sensors for placement

        div = 100 * np.abs((SpRx[idx] - SpRI) / SpRx[idx])  # Convergence of the last sensor placed
        SpRI = max(np.abs(np.linalg.eigvals(I)))

    Prcnt_SpRI = SpRI / SpRI_model * 100  # The % of SpRI sensor selection to the total of the model
    return Sensor_Node_Position, Prcnt_SpRI, div, SpRI

######################## MAIN PROGRAM ################################################################

# Start the timer
start_time = time.time()

# Initialize variables
counter , div = 0 , 101

# Main loop for sensor placement
for Ns in range(Ns_initial,Ns_maximum):
    Ts = Ts_initial # starting Ts threshold from 0 on each loop
    with tqdm(total=100, desc=f"{Ns} Sensors Configuration \nTs (Redundancy Threshold) sweeping ", unit="iteration") as pbar:
        while Ts < 1 and div > converg:
            Sensor_Node_Position,Prcnt_SpRI,div,SpRI = sensor_placement(Ts, Ns, Dim)
            Ts += 0.05
            counter += 1
            pbar.update(5)  # Assuming each step represents 5% progress
        if div <= converg:
            break
  

# Calculate the execution time
execution_time = time.time() - start_time

# Print the results in a formatted table
table = PrettyTable()
table.title = "Minimum and Optimal Sensor Placement Results"
table.field_names = ["Description", "Value"]

table.add_row(["The Sensor Positions are", [x + 1 for x in Sensor_Node_Position]])
table.add_row(["Number of sensors", len(Sensor_Node_Position)])
table.add_row(["Percentage of SpRI / SpRI(model)", f"{Prcnt_SpRI:.2f}%"])
table.add_row(["Convergence of the last placed sensor", f"{div:.6f}"])
table.add_row(["Execution time in seconds", f"{execution_time:.2f}"])
table.add_row(["Number of loops", counter])
table.add_row(["Ts%", f"{round(Ts - 0.05, 4)*100}%"])
table.add_row(["SpRI", SpRI])

print(table)

