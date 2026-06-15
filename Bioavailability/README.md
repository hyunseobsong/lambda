Here we provide a script for calculating bioavailability of Organic Matter from FTICR-MS-derived molecular formula using an extended form of the lambda model with incorporation of terms capturing traits identified by Ahamed et al. (2023).
This script can be called via terminal/command line with python, using the following parameters:

python.exe calc_bioav.py --input "C:\path\to\OM.csv" --output "C:\path\to\output_bioav.csv" --molform-col MolForm --pH  7.0 --Vh  1.0  --OC  1.0

By specifying the input and output paths, the molecular formula column name can be specified (--molform-col), as well as the pH (--pH), harvest volume (--Vh), and concentration of OM (--OC).
Concentration of Organic Matter can be a singular value, making the assumption of uniform concentrations, or if concentrations are known, may be provided as a column of compound-specific concentrations labeled "OM_Concentration" in the input csv.

Currently this script is limited to Aerobic respiration, and has not been extended to alternate electron acceptors.
