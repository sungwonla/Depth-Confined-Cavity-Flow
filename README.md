# Depth-Confined-Cavity-Flow
MATLAB code for simulating the depth-confined 3D lid-driven cavity flow problem in the Brinkman limit (Ault Research Group)

This code gives a numerical solution for the classic lid-driven cavity-flow problem, but under the special condition that the fluid flow is highly confined in the depth direction. The derivation of the governing equations used in the simulation are given in the file "Brinkman_Derivation.pdf", and result in quasi-2D depth-averaged equations for streamfunction, from which the depth-averaged velocities can be found. The simulation works well for the interior points of the mesh, but needs refinement near the boundaries. 
