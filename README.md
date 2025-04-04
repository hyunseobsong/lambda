# Thermodynamic modeling based on the lambda theory

## Release notes
* Substrate-Explicit Thermodynamic Modeling (SXTM; also known as lambda modeling) enables the formulation of stoichiometric and dynamic biogeochemical models based on the molecular formulas of organic compounds.
* The original model was developed with a focus on O₂ as the terminal electron acceptor (Song et al., Frontiers in Microbiology, 2020) and has since been extended to accommodate diverse electron acceptors, including but not limited to O₂ (Veeramani et al., 2025; Zheng et al., 2024).
The original lambda model (O₂-based) was implemented in R, while the extended model (supporting multiple electron acceptors) has been developed in both Python and R, referred to as lambdaPy and lambdaR, respectively.
* This repository contains the most up-to-date versions of the code for calculating key thermodynamic functions and parameters—such as Gibbs free energy changes and the lambda parameter—as well as the stoichiometries of catabolic, anabolic, and metabolic reactions.
* Which code should I use? For analyses involving aerobic respiration (O₂ as the electron acceptor), please use the code in [this folder]; for analyses involving O2 or alternative electron acceptors, please use the code in [this folder].
Click the Code icon (in green) to download the <b>latest version</b> designed for general purpose. 
* <b>[Code for the Frontiers paper](https://github.com/hyunseobsong/lambda/releases/tag/v0.1)</b>: A repository for the publication entitled ``Representing Organic Matter Thermodynamics in Biogeochemical Reactions via Substrate-Explicit Modeling`` (https://doi.org/10.3389/fmicb.2020.531756). In this paper, we proposed a new concept of biogeochemical modeling (termed substrate-explicit modeling) that enables parameterizing organic matter (OM)-specific oxidative degradation pathways and reaction rates based on the thermodynamic properties of OM pools.

## Links
* <b>Substrate-explicit Modeling Tutorial</b> ([KBase narrative](https://kbase.us/n/65526/69/)): It demonstates how to apply a thermodynamic theory (known as λ theory) to convert molecular formulas of compounds (derived from FTICR-MS peaks using [Formularity](https://omics.pnl.gov/software/formularity)) into stoichiometric and kinetic forms of biogeochemical reactions and how to use the resulting kinetic equations to simulate dynamic conversion of compounds in batch and continuous stirred tank reactors.
* KBase app to build stoichiometric reaction models from chemical formulas ([KBase link](https://narrative.kbase.us/#appcatalog/app/ThermoStoichWizard/run_ThermoStoichWizard/), [GitHub](https://github.com/coldfire79/ThermoStoichWizard))
* KBase app to perform λ analysis ([KBase link](https://narrative.kbase.us/#catalog/apps/ThermoStoichWizard/run_lambda_analysis/), [GitHub](https://github.com/coldfire79/ThermoStoichWizard))
* KBase app to simulate biogeochemical models in a batch reactor ([KBase link](https://narrative.kbase.us/#appcatalog/app/BatchBiogeochemicalReactionModel/run_BatchBiogeochemicalReactionModel/), [GitHub](https://github.com/coldfire79/BatchBiogeochemicalReactionModel))
* KBase app to simulate biogeochemical models in a continuous stirred tank reactor (CSTR) ([KBase link](https://narrative.kbase.us/#catalog/apps/BatchBiogeochemicalReactionModel/run_cstr/), [GitHub](https://github.com/coldfire79/BatchBiogeochemicalReactionModel))

## System requirements
R > 3.5.0


## Citation
Song H-S, Stegen JC, Graham EB, Lee J-Y, Garayburu-Caruso VA, Nelson WC, Chen X, Moulton JD and Scheibe TD (2020) Representing Organic Matter Thermodynamics in Biogeochemical Reactions via Substrate-Explicit Modeling. <i>Front. Microbiol.</i> 11:531756. doi: 10.3389/fmicb.2020.531756

## Contacts
Hyun-Seob Song (hsong5@unl.edu); Joon-Yong Lee (joonyong.lee@pnnl.gov)
