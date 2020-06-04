Final model for predicting electrophilicity and nucleophilicity


For electrophilicity:
Your file should be in this format,
------------------------------------------------------
Me2N+=CH-CH-Ph  C/[N+](C)=C\C=C\C1=CC=CC=C1  3.733005743 0.43441  0.27663  0.1648   0.14684  0.1378   0.13379  0.01274  0.01206 0.00497  0.00225  0.00189  0.00182  0.00178  0.00075  0.0007  -0.00429  -0.00449 -0.00514 -0.0057  -0.00786 -0.00834 -0.01295 -0.06172 -0.06483 -0.07758 -0.08034
ani+Br+2QM BrC(C(C(Br)=C/1)=O)=CC1=C\C2=CC=C(OC)C=C2  3.406799513  4.8184e-01  1.4127e-01  1.3386e-01  1.3369e-01  1.1928e-01  1.1719e-01 1.0839e-01  7.1530e-02  2.6640e-02  1.9600e-03  1.8900e-03  1.8300e-03 1.7800e-03 -1.0000e-04 -3.8000e-04 -1.5900e-03 -1.7400e-03 -2.3900e-03  -3.5100e-03 -3.5700e-03 -3.9500e-03 -4.5400e-03 -1.4500e-02 -2.3710e-02  -6.1190e-02 -6.2170e-02 -7.1040e-02 -8.6780e-02
...
------------------------------------------------------
  (Name)      (SMILES)   (Adiabatic E.A.)  (Parr function, sorted in descending order)

Name -- Name of the molecule
SMILES -- SMILES expression of molecular structure
Adiabatic E.A. -- Electron affinity calculated from fully relaxed (n+1) electron structure
Parr function -- Spin density of fully relaxed (n+1) electron structure, sorted in descending order. 
                 ~30 Parr function values can be provided. write down to the 30th value if the number of atom in molecule is larger than 30.


For nucleophilicity:
Your input file should be in this format,
------------------------------------------------------
4-aminopyridine  NC1=CC=NC=C1  3.04959102  0.45322  0.30205  0.264    0.26396  0.00329  0.00329  0.00054 -0.00799 -0.00799 -0.01196 -0.01196 -0.12522 -0.12522 
glycineamide   NCC(N)=O   2.847399712  0.89849  0.064    0.06399  0.03633 -0.00104 -0.00112 -0.00352 -0.00428 -0.00825 -0.02198 -0.02261
...
------------------------------------------------------
  (Name)      (SMILES)   (Adiabatic I.E.)  (Parr function, sorted in descending order)

Name -- Name of the molecule
SMILES -- SMILES expression of molecular structure
Adiabatic I.E. -- Ionization energy calculated from fully relaxed (n-1) electron structure
Parr function -- Spin density of fully relaxed (n-1) electron structure, sorted in descending order. 
                 ~30 Parr function values can be provided. write down to the 30th value if the number of atom in molecule is larger than 30.
