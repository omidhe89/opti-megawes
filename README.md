# opti-megawes
An Optimization Framework for Path Planning of Megawatt-Scale Airborne Wind Energy Systems.

Opti-megawes is a MATLAB-based numerical software package developed to solve constrained nonlinear optimal control problems using direct optimization techniques. 
It employs a direct multiple shooting method combined with an implicit Euler discretization scheme to compute optimal flight trajectories, including both circular 
and lemniscate (figure-eight) paths. The AWES model is based on the model used in **MegAWES [1]**, an open-source MATLAB Simulink toolbox, which was developed to evaluate 
the performance and control of an AWES with a forty-meter wing span aircraft. A minimal coordinate model, in which the aircraft is represented by Euler angles, is 
utilized, and the position propagation of the system is described in spherical coordinates. Additionally, a Penalty-based Interior-Point Homotopy (PIPH) approach, similiar to **AWEBOX [2]** is 
adopted to increase the efficiency and robustness of solving the OCP. The tool offers different aerodynamic modelling options. These options correspond to the manner 
by which the multiple aerodynamic coefficients are computed, including the **Vortex Lattice Method (VLM)[3]**, **Actuator Line Method (ALM) in Large Eddy simulation [4]**, and using 
**body resolved Computational Fluid Dynamics (CFD) [5]**.

## Prerequisites
Before you begin, ensure you have the following installed and configured:

* **MATLAB**: This project uses MATLAB scripts to automate the LaTeX compilation process and interact with CasADi.

* **CasADi**: A symbolic framework for numeric optimization. You will need to install the MATLAB interface for CasADi.
    * **Installation**: Download the CasADi `.zip` file from [https://web.casadi.org/get](https://web.casadi.org/get). Extract its contents to a desired folder on your system.
    * **MATLAB Setup**: In MATLAB, navigate to the folder where you extracted the CasADi files. Then, run the following commands in the MATLAB Command Window:
        ```matlab
        addpath(genpath('<path/to/the/folder/in/which/CasADi/files/are/extracted/to>'));
        import casadi.*;
        ```
        **Important**: Replace `<path/to/the/folder/in/which/CasADi/files/are/extracted/to>` with the actual path to your extracted CasADi directory.
    * **Verification**: To check CasADi installation, try executing the following lines in the MATLAB Command Window:
        ```matlab
        x = MX.sym('x');
        disp(jacobian(sin(x),x));
        ```
        If any error is thrown, this implies an irregular installation of CasADi, and the user is advised to follow the above instructions again closely. Also, please keep in mind that CasADi might have to be imported to MATLAB every time MATLAB is restarted, depending on the user's system settings.

* **LaTeX Distribution**: You'll need a full LaTeX distribution (like TeX Live, MiKTeX for Windows, or MacTeX for macOS) installed on your system.
    * **Important**: Make sure the LaTeX compiler (e.g., `pdflatex`, `xelatex`) is added to your system's `PATH` environment variable. This allows MATLAB to find and execute the LaTeX commands.

## Author
Omid Heydarnia
omid.heydarnia@ugent.be

## References
1.  Eijkelhof, D. and R. Schmehl, Six-Degrees-Of-Freedom Simulation Model for Future Multi-Megawatt Airborne Wind Energy Systems. Renewable Energy, 2022. 196.
2.  De Schutter, J., et al., AWEbox: An Optimal Control Framework for Single- and Multi-Aircraft Airborne Wind Energy Systems. Energies, 2023. 16(4): p. 1900.
3.  Drela, M.a.H.Y. {AVL} 3.40b. 2024; Available from: http://web.mit.edu/drela/Public/web/avl/.
4.  Crismer, J.-B., et al. Large-Eddy Simulation of airborne wind energy systems wakes. in Journal of Physics: Conference Series. 2023. IOP Publishing.
5.  Pynaert, N., et al., Wing Deformation of an Airborne Wind Energy System in Crosswind Flight Using High-Fidelity Fluid–Structure Interaction. Energies, 2023. 16(2): p. 602.