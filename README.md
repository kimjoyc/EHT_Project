# EHT_Project

## Abstract 
The paper presents a modified Extended Hückel method that provides accurate values of heats of formation and structural parameters for hydrocarbons. The method is computationally efficient and can predict reaction enthalpies for various hydrocarbon transformations. The results of the method are compared to experimental data for a wide range of molecules and radicals, showing good agreement with a mean absolute error of 1.90 kcal/mol.

## Method
The total energy of the system is calculated as the sum of the electronic energy and short-range repulsion energy. The electronic energy is obtained by solving eigenvalue equations, while the repulsion energy is approximated using two-center potentials. The method employs an orthogonal atomic orbital basis for computational efficiency.

![IMG_5773 (2)](https://github.com/kimjoyc/EHT_Project/assets/88675769/785474e3-6c67-4f02-9308-bc77a326721d)

## Parameters

The parameters of the effective Hamiltonian and the repulsion potential are derived by fitting them to experimental heats of formation and structural parameters for a training set of molecules and radicals. The fitted parameters are then used to calculate the heats of formation for various hydrocarbons.

![IMG_5770](https://github.com/kimjoyc/EHT_Project/assets/88675769/edcf7e9b-56ba-4fbb-a8d4-ee7d1cececc0)

![IMG_5775 (2)](https://github.com/kimjoyc/EHT_Project/assets/88675769/9fe6d70d-efc3-44e9-bc8d-c8b43c6f1965)

![IMG_5774 (2)](https://github.com/kimjoyc/EHT_Project/assets/88675769/d9785a96-611a-4bb5-841b-2d0d3cacc0b0)

## Molecule Data Structure 
### Each Atom has an orbital that is in 3D 
EX: Hydrogen has on s-orbital but each paramters are different for S-orbital in x dimension, S-orbital in y dimension, and S-orbital in z dimension: 
The Data Structure for the Atom Information is compirsed of each file containing Atomic Orbital of that Atom like Hydrogen S-orbital in in the x,y,z order from top to bottom. 

Here is an example of each paramter: 
elem_num,x0, y0, z0, alpha,d, l0,l1,l2,gamma,beta,n_s_orb,n_p_orb,U_s_orb,U_p_orb,heat_of_formation,n_s_orb_tot,n_p_orb_tot

S-orbital in the x direction: 
elem_num= 1

Cartesian Coordinate info: 
x0=0.0
y0=0.0 
z0=0.0 

Exponents and contraction coefficients for STO-3G for the standard light elements(in atomic units):

exponent:
alpha=3.42525091 

contraction coeffiecient for s-orbital:
d=0.15432897

angular momentum for s-orbital:
l0=0
l1=0
l2=0

Semi-empirical parameters that define the CNDO/2 model:

for H:
gamma=-7.176

for H:
beta=-9

Number of s-orbital
n_s_orb=1

Number of p-orbital
n_p_orb=0

Atomic Parameters Chart in eV:
U_s_orb=-13.605
U_p_orb=0

for Hydrogen: 
heat_of_formation=52.1 kcal/mol 

total s-orbital
n_s_orb_tot=1

total p-orbital
n_p_orb_tot=0

S-orbital x direction:
Sx_info:
1 0.0 0.0 0.0 3.42525091 0.15432897 0 0 0 -7.176 -9 1 0 -13.605 0 52.1 1 0 


S-orbital y direction:
Sy_info:
1 0.0 0.0 0.0 0.62391373 0.53532814 0 0 0 -7.176 -9 1 0 -13.605 0 52.1 1 0


S-orbital z direction:
Sz_info:
1 0.0 0.0 0.0 0.16885540 0.44463454 0 0 0 -7.176 -9 1 0 -13.605 0 52.1 1 0


## Compile Instructions
Fill in Paramters as follows: 
```
bash eht.sh 
```

## Utilties Function 
Factorial and Combinations 

## Conclusion / Results

### Methodology Challenges
The main challenge in implementing the code was related to setting up the H-matrix correctly. Since there are multiple orbitals belonging to carbon atoms, it was necessary to filter the parameters based on the type of atomic bonds (C-H, C-C, H-H) and the interacting orbitals. The approach involved using atomic numbers and angular momentum to separate the different possibilities of C-H bond pairs (ss or sp-pairs) and C-C pairs (pp-pairs). However, as the molecules became larger, it became apparent that the H-matrix was not properly set up to solve for the correct eigenvalues, leading to inaccuracies in the results.

### Results and Limitations
The code implementation only yielded a result for diatomic hydrogen gas, with a calculated heat of formation of 0.104 kcal/mol, slightly deviating from the paper's value of 0.0 kcal/mol. Additionally, using different coordinate systems (atomic units instead of NIST angstrom coordinates) resulted in a slight discrepancy in the calculation for hydrogen gas.

Unfortunately, due to limited information and missing Cartesian coordinate data for larger molecules in the dataset, it was challenging to obtain accurate heat of formation values for the 47 molecules listed in Table 2 of the paper. The paper's source of this information is unclear. The limitations in accessing the necessary data hindered further progress in the project.

### Future Improvements
To overcome the coding challenges, future improvements can be made by reevaluating the data structure and ensuring the availability of accurate and complete data. Implementing internuclear distances instead of relying solely on Cartesian coordinates could lead to more accurate results. Additionally, a potential modification could involve using a mapping approach instead of filtering conditions for atomic orbital pairs to improve the accuracy and efficiency of the code.

Considering the limitations and incomplete data available, the project was concluded at this stage. It should be noted that obtaining the total energy of organic molecules was challenging, with only heat of combustion data (e.g., ethane) being readily available online.
