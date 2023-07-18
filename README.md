# EHT_Project

## Abstract 
The paper presents a modified Extended Hückel method that provides accurate values of heats of formation and structural parameters for hydrocarbons. The method is computationally efficient and can predict reaction enthalpies for various hydrocarbon transformations. The results of the method are compared to experimental data for a wide range of molecules and radicals, showing good agreement with a mean absolute error of 1.90 kcal/mol.

## Method
The total energy of the system is calculated as the sum of the electronic energy and short-range repulsion energy. The electronic energy is obtained by solving eigenvalue equations, while the repulsion energy is approximated using two-center potentials. The method employs an orthogonal atomic orbital basis for computational efficiency.

## Parameters



The parameters of the effective Hamiltonian and the repulsion potential are derived by fitting them to experimental heats of formation and structural parameters for a training set of molecules and radicals. The fitted parameters are then used to calculate the heats of formation for various hydrocarbons.

## Molecule Data Structure 
### Each Atom has an orbital that is in 3D 
EX: Hydrogen has on s-orbital but each paramters are different for S-orbital in x dimension, S-orbital in y dimension, and S-orbital in z dimension: 
The Data Structure for the Atom Information is compirsed of each file containing Atomic Orbital of that Atom like Hydrogen S-orbital in in the x,y,z order from top to bottom. 

Here is an example of each paramter: 
elem_num,x0, y0, z0, alpha,d, l0,l1,l2,gamma,beta,n_s_orb,n_p_orb,U_s_orb,U_p_orb,heat_of_formation,n_s_orb_tot,n_p_orb



1 0.0 0.0 0.0 3.42525091 0.15432897 0 0 0 -7.176 -9 1 0 -13.605 0 52.1 1 0 
1 0.0 0.0 0.0 0.62391373 0.53532814 0 0 0 -7.176 -9 1 0 -13.605 0 52.1 1 0
1 0.0 0.0 0.0 0.16885540 0.44463454 0 0 0 -7.176 -9 1 0 -13.605 0 52.1 1 0


## Compile Instructions
Fill in Paramters as follows: 
```
bash eht.sh 
```

## Utilties Function 
Factorial and Combinations 
