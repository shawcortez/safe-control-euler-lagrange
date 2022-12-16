# Introduction

This repository contains the Python code for constructing safe-by-design control laws for Euler-Lagrange systems. The code uses properties of Euler-Lagrange systems to formally construct zeroing control barrier functions (ZCBFs) to guarantee safety (i.e. position and velocity constraints), while always respecting input constraints. The ZCBFs are automatically constructed and the associated algorithm takes into account robustness margins and sampling time effects. Additionally, two controllers are included. One is a continuous-time QP-based control law to enforce safety under the assumption that the control is locally Lipschitz. Another, sample-data based QP control law is provided to handle sampling time effects and account for the discrete implementation of digital controllers while still ensuring safety. The control laws attempt to track an existing, nominal control law while always satisfying the ZCBF conditions for safety. The system dynamics used in this code are for a 2DOF planar manipulator, however different system dynamics can be incorporated by creating an appropriate dynamics module. See the two_dof_module.py for an example.

The results from this code and associated analysis/guarantees are found in the manuscript entitled "Safe-by-Design Control for Euler-Lagrange Systems" by Wenceslao Shaw Cortez and Dimos V. Dimarogonas, which can be found here: https://doi.org/10.1016/j.automatica.2022.110620. Please reference this manuscript when using this code.

# Technologies

- The Python version used to develop this code is Python 2.7.17
- Required modules:
  - numpy
  - matplotlib
  - yaml
  - sys
  - math
  - cvxopt (https://cvxopt.org/)
  - scipy.integrate
  
# Setup

To install this project, assuming the appropriate Python/modules are installed, simply git clone this repository.

# Use

## Initial run

To use the code, first try running:

python zcbf_simulation_main.py exp2

The code should output a message stating "Running experiment: ZCBF_control_exp2.yaml" and proceed to compute the ZCBF parameters to be used. Finally, the simulation should be run and a plot of the state and input trajectories should appear. Once the plot is closed, a message should appear asking to save the data. Enter 'no' to avoid saving.

## Overview

The two main scripts are "zcbf_simulation_main.py" and "two_dof_dynamics_bounds_main", and the simulation parameters are included in the yaml files. Each yaml file uses the following naming convention: "ZCBF_control_xxxx.yaml". The script "plot_data.py" develops the plots used in the accompanying manuscript. Simply toggle the comments at the beggining of the file to create the two main plots of the manuscript.

To create new simulations using the same 2DOF dynamics module (see two_dof_module.py), create a new .yaml file and copy the contents of one of the existing yaml files. If creating a new simulation with new model parameters OR larger delta value, new system bounds must be computed. To do this save the changes in the new yaml file and run:

python two_dof_dynamics_bounds_main.py xxx

where xxxx is the associated name in the new yaml file: "ZCBF_control_xxxx.yaml". For example, if the new yaml file is titled "ZCBF_control_exp3_new.yaml", the command to run the simulation is: 

python two_dof_dynamics_bounds_main exp3_new

The code should output all the relevant system bounds: km2, kc, kg, c3, c5, c1. Note that c1 is notoriously long to compute (about 40 min for the 2DOF model). Copy these bounds into the respective locations of the yaml file. Most system bound computations are done numerically, i.e., via grid search, by checking over the states defined by the constraint sets. To adjust the step size used in the grid search, adjust the dq_bounds and dv_bounds terms appropriately in the yaml file.

Once the new yaml file is updated with the appropriate system bounds, the simulation is ready to run. To run the new simulation type: 

python zcbf_simulation_main xxxx

where again xxxx is associated with the new yaml file name.

## Additional Info

The yaml files contain all the necessary information to run the simulations. The most important parameters are the "compute_zcbf" and "sim_control".

"compute_zcbf" is a boolean that dictates if the ZCBF construction algorithm will be run ('True'). In doing so, all the ZCBF parameters in the yaml file will be treated as initial parameters for the algorithm. The associated ZCBF design parameters are: "gamma", "nu", "delta", "eta_bar", "alpha", and "beta". The terms "eps_frac" and "nu_frac" are used to define how the epsilon and nu ZCBF parameters are to be chosen. If "compute_zcbf" is set to 'False', the simulation will run with the ZCBF parameters defined in the yaml file and the algorithm for constructing the ZCBF parameters will not run. Note the "alpha" and "beta" terms represent extended class K_infty functions and currently the only supported functions for alpha are: alpha(h) = h, alpha(h) = atan(h), and for beta: beta(h) = h, beta(h) = h³. New extended class K_infty functions can be added at the beginning of "ZCBF_module.py" in the "ExtendedClassK" class.

The parameter "sim_control" dictates which control implementation to run. When set to "cont_time" the continuous time controller is used, and when set to "discrete_time" the sampled-data controller is used. Note that T_sample defines the sampling period and so is only used for the sampled-data controller.

The control laws run by attempting to track a reference "nominal" control law. This is used to represent a human operating the system or an existing control law. Currently, the nominal control parameter "nom_control" can either be a "PD" control law or "ComputedTorque". For the "ComputedTorque" controller, the associated reference trajectory is listed in the existing yaml files.











