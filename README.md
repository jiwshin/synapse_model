# synapse_model
Simulation code supporting Shin et al. (2025), eLife

This repository contains MATLAB scripts for simulating synaptic vesicle dynamics, specifically focusing on Synaptotagmin-7 (Syt7) mediated short-term plasticity (STP).

## File Descriptions

### Main Simulation
*   **`ee_2.m`**: The **primary execution script**.
    *   Runs the simulation of synaptic release across multiple trials and stimulation frequencies (e.g., 5Hz, 10Hz, 20Hz, 40Hz).
    *   Calls `ni7_2.m` for the core numerical integration.
    *   Generates key figures.

### Core Functions
*   **`ni7_2.m`**: The numerical integration kernel.
    *   Simulates the time-dependent dynamics of synaptic vesicles, calcium concentration (global and local Gaussian transients), and Syt7 binding states.
    *   Updates state matrices based on transition probabilities (docking, priming, fusion).
*   **`setparam2.m`**: Parameter initialization module.
    *   Sets global parameters, rate constants, and initial conditions for the simulation based on an input parameter vector.
    *   Configures Syt7 forward/backward rates.
*   **`initSyt.m`**: Helper function to calculate the initial steady-state distribution of Syt7 calcium-bound states.

## Usage

1.  **Setup**: Ensure all `.m` files are in the same directory or MATLAB path.
2.  **Run**: Execute the **`ee_2.m`** script in MATLAB.
    *   The script uses default parameters defined within `ee_2.m` and `setparam2.m` if external data tables (`rcvdata`, `stpdata`) or parameter files are not provided.
3.  **Output**:
    *   **Figure 9**: STP curves comparing simulation results with data.
    *   **`mc_.mat`**: Saves the fitted quantal parameters (`Nfit`, `qfit`) and simulation statistics.
