# Substrate-Explicit Thermodynamic Modeling

## Release notes
* Substrate-Explicit Thermodynamic Modeling (SXTM; also known as lambda modeling) enables the formulation of stoichiometric and dynamic biogeochemical models based on the molecular formulas of organic compounds.
* The original model was developed with a focus on O₂ as the terminal electron acceptor (Song et al., Frontiers in Microbiology, 2020) and has since been extended to accommodate diverse electron acceptors, including but not limited to O₂ (Veeramani et al., 2025; Zheng et al., 2024).
The original lambda model (O₂-based) was implemented in R, while the extended model (supporting O2 or alternative electron acceptors) has been developed in both Python and R, referred to as LambdaPy and LambdaR, respectively.
* This repository contains the most up-to-date versions of the code for calculating key thermodynamic functions and parameters—such as Gibbs free energy changes and the lambda parameter—as well as the stoichiometries of catabolic, anabolic, and metabolic reactions.
* Which code should I use? For analyses involving aerobic respiration (O₂ as the electron acceptor), please use the code in 
<b>Original code for aerobic respiration</b>:; for analyses involving O2 or alternative electron acceptors, please use the code in <b>Extension to diverse electron acceptors</b>:.

## Citation
Song H-S, Stegen JC, Graham EB, Lee J-Y, Garayburu-Caruso VA, Nelson WC, Chen X, Moulton JD and Scheibe TD (2020) Representing Organic Matter Thermodynamics in Biogeochemical Reactions via Substrate-Explicit Modeling. <i>Front. Microbiol.</i> 11:531756. doi: 10.3389/fmicb.2020.531756

Zheng J, Scheibe TD, Boye K, and Song HS (2024) Thermodynamic control on the decomposition of organic matter across different electron acceptors, Soil Biology and Biochemistry, 193, 109364.

Veeramani M, Kharel S, McCullough HC, Chen X, Zheng J, Stegen JC, Scheibe TD, and Song HS (2025) LambdaPy and LambdaR: Thermodynamics-Based Biogeochemical Reaction Modeling Packages for Integrating High-Resolution Mass Spectrometry Data, Modeling the Microbiome, Edited by Raman K and Sambamoorthy G, Lab Protocol Series: Methods in Molecular Biology, Springer Nature (IN press) 

## Contacts
Hyun-Seob Song (hsong5@unl.edu)

