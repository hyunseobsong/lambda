#LambdaR script

#Assign your input parameters here ==================================================================

#Specify the location of your input csv, containing columns for each element CHNOPS and rows for each compound
input_data <- read.csv('C:/Users/hmccullough2/Downloads/demo_input.csv')

#Specify your electron acceptor value, from the list below the input parameters section
#Default is oxygen, with value 0
electronAcceptor <- 0

# Assign your pH, numeric values 0-14 
pH <- 7

#Specify your output destination and file name

output_file <- "C:/Users/hmccullough2/Desktop/lambdaR_out.csv"

#=====================================================================================================

#Electron acceptor choices ===========================================================================

# 0)  Oxygen (Default set to O2 as e acceptor) 0 (Default). O2 + 4H+ + 4e- ----\> 2H2O
#   
#   1)  Nitrogen Compounds, Nitrates and Nitrites
# 
#   -   1.1. NO3- + 10H+ + 8e- ---\> NH4+ + 3H2O 
#   
#   -   1.2. NO3- + 2H+ + 2e- ---\> NO2- + H2O
#   
#   -   1.3. NO3- + 6H+ + 5e- ---\> (1/2)N2 + 3H2O
#   
#   -   1.4 NO2- + 4H+ + 3e- ---\> (1/2)N2 + 2H2O
#   
#   -   1.5. NO2- + 8H+ + 6e- ---\> NH4+ + 2H2O
#   
#   -   1.6. N2 + 8H+ + 6e- ---\> 2NH4+
#   
#   2)  Sulphur compounds, Sulphates and Sulphites
# 
#   -   2.1. SO4\^2- + (9/2)H+ + 8e- ---\> (1/2)H2S + (1/2)HS- + 4H2O
#   
#   -   2.2 SO4\^2- + 2H+ + 2e- ---\> SO3\^2- + H2O
#   
#   -   2.3. SO4\^2- + 5H+ + 4e- ---\> (1/2)S2O3\^2- + (5/2)H2O
#   
#   -   2.4. SO4\^2- + 8H+ + 6e- ---\> S + 4H2O
#   
#   -   2.5. SO4\^2- + 9H+ + 8e- --\> HS- + 4H2O
#   
#   -   2.6. SO3\^2- + (15/2)H+ + 6e- ---\> (1/2)H2S + (1/2)HS- + 3H2O
#   
#   3)  Iron compounds, ferrous and ferric
# 
#   -   3.1. Fe(OH)3 + 3H+ + e- --\> Fe2+ + 3H2O
#   
#   -   3.2. FeOOH + 3H+ + e- --\> Fe2+ + 2H2O
#   
#   -   3.3. Fe3O4 + 8H+ + 2e- --\> 3Fe2+ + 4H2O
#   
#   -   3.4. Fe3+ + e- ---\> Fe2+
#   
#   4)  Bicarbonate and Hydrogen ion
# 
#   -   4.1. HCO3- + 9H+ + 8e- --\> CH4 + 3H2O
#   
#   -   4.2. H+ + e- ---\> (1/2)H2
#   
#   5)  Acetate
# 
#   -   5 CH3COO- + 9H+ + 8e- --\> 2CH4 + 2H2O
#   
#   6)  Manganese
# 
#   -   6 MnO2 + 4H+ + 2e- --\> Mn2+ + 2H2O


#=====================================================================================================



#Installation and Loading of necessary libraries =====================================================
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  # If not installed, install 'tidyverse'
  install.packages("tidyverse")
}

library(tidyverse)
#=====================================================================================================
CHEMICAL_ELEMENTS <- c("C", "H", "N", "O", "P", "S")

get_compositions <- function(df) {
  
  chemical_compositions <- NULL
  formulas <- NULL
  
  if ("C" %in% colnames(df)) {
    tdf <- df[df[["C"]] > 0, ]
    if ("C13" %in% colnames(df)) {
      tdf <- tdf[tdf[["C13"]] == 0, ]
    }
    chemical_compositions <- as.matrix(tdf[CHEMICAL_ELEMENTS])
    formulas <- tdf[["MolForm"]]
  } else {
    stop("Either columns for compositions (e.g., C, H, N, ...) or `MolForm` column is required.")
  }
  
  if ("Z" %in% colnames(df)) {
    chemical_compositions <- cbind(chemical_compositions, df[["Z"]])
  } else {
    chemical_compositions <- cbind(chemical_compositions, rep(0, nrow(df)))
  }
  
  list(chemical_compositions = chemical_compositions, formulas = formulas)
}

getThermoStoich <- function(chemForm, ui, pH_Input) {
  a <- chemForm[1]
  b <- chemForm[2]
  c <- chemForm[3]
  d <- chemForm[4]
  e <- chemForm[5]
  f <- chemForm[6]
  z <- chemForm[7]
  
  # Step 1a) stoichD: stoichiometries for an electron donor=====================
  ySource <- -1
  yH2o <- -(3*a + 4*e - d)
  yHco3 <- a
  yNh4 <- c
  yHpo4 <- e
  yHs <- f
  yH <- 5*a + b - 4*c - 2*d + 7*e - f
  yE <- -z + 4*a + b - 3*c - 2*d + 5*e - 2*f
  stoichD <- c(ySource, yH2o, yHco3, yNh4, yHpo4, yHs, yH, yE, rep(0, 20))
  stoichA <- numeric(28)
  stoichA <- rep(0, 28)
  
  # Step 1b) stoichA: stoichiometries for an electron acceptor =================
  if (ui == 0) {
    stoichA[c(9, 7, 8, 2)] <- c(-1, -4, -4, 2)
  } else if (ui == 1.1) {
    stoichA[c(10, 7, 8, 4, 2)] <- c(-1, -10, -8, 1, 3)
  } else if (ui == 1.5) {
    stoichA[c(11, 7, 8, 4, 2)] <- c(-1, -8, -6, 1, 2)
  } else if (ui == 1.6) {
    stoichA[c(12, 7, 8, 4, 2)] <- c(-1, -8, -6, 2, 0)
  } else if (ui == 3.4) {
    stoichA[c(13, 7, 8, 14, 2)] <- c(-1, 0, -1, 1, 0)
  } else if (ui == 4.2) {
    stoichA[c(7, 8, 15)] <- c(-1, -1, 0.5)
  } else if (ui == 1.2) {
    stoichA[c(10, 7, 8, 11, 2)] <- c(-1, -2, -2, 1, 1)
  } else if (ui == 1.3) {
    stoichA[c(10, 7, 8, 12, 2)] <- c(-1, -6, -5, 0.5, 3)
  } else if (ui == 1.4) {
    stoichA[c(11, 7, 8, 12, 2)] <- c(-1, -4, -3, 0.5, 2)
  } else if (ui == 2.1) {
    stoichA[c(16, 7, 8, 17, 6, 2)] <- c(-1, -4.5, -8, 0.5, 0.5, 4)
  } else if (ui == 2.6) {
    stoichA[c(18, 7, 8, 17, 6, 2)] <- c(-1, -7.5, -6, 0.5, 0.5, 3)
  } else if (ui == 2.2) {
    stoichA[c(16, 7, 8, 18, 2)] <- c(-1, -2, -2, 1, 1)
  } else if (ui == 2.4) {
    stoichA[c(16, 7, 8, 19, 2)] <- c(-1, -8, -6, 1, 4)
  } else if (ui == 2.3) {
    stoichA[c(16, 7, 8, 20, 2)] <- c(-1, -5, -4, 0.5, 2.5)
  } else if (ui == 3.1) {
    stoichA[c(21, 7, 8, 14, 2)] <- c(-1, -3, -1, 1, 3)
  } else if (ui == 3.2) {
    stoichA[c(22, 7, 8, 14, 2)] <- c(-1, -3, -1, 1, 2)
  } else if (ui == 3.3) {
    stoichA[c(23, 7, 8, 14, 2)] <- c(-1, -8, -2, 3, 4)
  } else if (ui == 2.5) {
    stoichA[c(16, 7, 8, 6, 2)] <- c(-1, -9, -8, 1, 4)
  } else if (ui == 4.1) {
    stoichA[c(3, 7, 8, 24, 2)] <- c(-1, -9, -8, 1, 3)
  } else if (ui == 5) {
    stoichA[c(25, 7, 8, 24, 2)] <- c(-1, -9, -8, 2, 2)
  } else if (ui == 6) {
    stoichA[c(26, 7, 8, 27, 2)] <- c(-1, -4, -2, 1, 2)
  } else {
    stop("Invalid value of ui. Please provide a valid value.")
  }
  
  # Step 1c) stoichCat: stoichiometries for catabolic reaciton ================= 
  yEd <- stoichD[8]
  yEa <- stoichA[8]
  stoichCat <- sapply(1:length(stoichD), function(i) stoichD[i]-(yEd/yEa)*stoichA[i])
  
  
  # Step 2a) stoichAnStar: stoichiometries for anabolic reaciton (N source = NH4+) =======
  chemFormBiom <- c(1, 1.8, 0.2, 0.5, 0, 0, 0)
  aB <- chemFormBiom[1]
  bB <- chemFormBiom[2]
  cB <- chemFormBiom[3]
  dB <- chemFormBiom[4]
  eB <- chemFormBiom[5]
  fB <- chemFormBiom[6]
  zB <- chemFormBiom[7]
  ySource <- -1
  yH2o <- -(3*aB+4*eB-dB)
  yHco3 <- aB
  yNh4 <- cB
  yHpo4 <- eB
  yHs <- fB
  yH <- 5*aB+bB-4*cB-2*dB+7*eB-fB
  yE <- -zB + 4*aB + bB - 3*cB - 2*dB + 5*eB - 2*fB
  stoichAnStarB <- c(ySource, yH2o, yHco3, yNh4, yHpo4, yHs, yH, yE, rep(0, 20))
  stoichAnStarB[9:28] <- rep(0, 20)
  stoichAnStarB <- -stoichAnStarB
  stoichAnStarB[28] <- stoichAnStarB[1]
  stoichAnStarB[1] <- 0
  
  # Step 2b) "overall" anabolic reaction =======================================
  stoichAnStar <- sapply(1:length(stoichAnStarB), function(i) stoichAnStarB[i] + (1 / a) * stoichD[i])
  yEana <- stoichAnStar[8]
  
  if (yEana > 0) {
    stoichAn <- sapply(1:length(stoichAnStar), function(i) stoichAnStar[i] - (yEana / yEa) * stoichA[i])
  } else if (yEana < 0) {
    stoichAn <- sapply(1:length(stoichAnStar), function(i) stoichAnStar[i] - (yEana / yEd) * stoichD[i])
  } else {
    stoichAn <- stoichAnStar
  }
  
  # Step 3: get lambda =========================================================
  # - estimate delGd0 using LaRowe and Van Cappellen (2011)
  ne <- -z + 4*a +b-3*c-2*d+5*e-2*f  # number of electrons transferred in D 
  nosc <- -ne/a + 4  # nominal oxidation state of carbon 
  delGcox0 <- 60.3-28.5*nosc  # kJ/C-mol
  delGd0 <- delGcox0*a*abs(stoichD[1])  # kJ/rxn
  
  # estimate delGf0 for electron donor
  delGf0_D_zero <- 0
  delGf0_zero <- c(delGf0_D_zero, -237.2, -586.9, -79.37, -1089.1, 12.05, 0, 0, 16.5, -111.3, -32.2, 18.19, -4.6, -78.87, 0, -744.63, -33.4, -486.6, 0, 522.5, -690, -489.8, -1012.6, -34.06, -392, -465.14, -228, -67)
  delGcox0_zero <- sum(delGf0_zero*stoichD)
  delGf0_D_est <- (delGd0-delGcox0_zero)/stoichD[1]  
  delGf0 <- delGf0_zero
  delGf0[1] <- delGf0_D_est
  
  # standard delG at pH=0
  delGcat0 <- sum(delGf0*stoichCat)
  delGan0 <- sum(delGf0*stoichAn)
  
  # standard delG at pH=7
  R <- 0.008314
  T <- 298.15
  iProton <- 7
  pH_value <- 10^(-pH_Input)
  delGd <- delGd0+(R*T*stoichD[iProton])*(log10(pH_value))
  delGcox <- delGd/a
  delGcat <- delGcat0+(R*T*stoichCat[iProton])*(log(pH_value))
  delGan <- delGan0+(R*T*stoichAn[iProton])*(log(pH_value))
  
  # The Thermodynamic Electron Equivalents Model (TEEM) ========================
  eta <- 0.43
  delGsyn <- 200
  
  if (is.nan(delGan0) && is.nan(delGan)) {
    lambda0 <- NaN
    lambda_ <- NaN
    stoichMet <- rep(NaN, length(stoichCat))
    delGdis0 <- NaN
    delGdis <- NaN
  } else {
    if (delGan < 0) {
      m <- 1
    } else {
      m <- -1
    }
    if (delGan0 < 0) {
      m0 <- 1
    } else {
      m0 <- -1
    }
    lambda0 <- ((delGan0*(eta^m0))+delGsyn)/(-delGcat0*eta)
    lambda_ <- (delGan*(eta^m)+delGsyn)/(-delGcat*eta)
    if (lambda_ > 0) {
      stoichMet <- numeric(length(stoichCat))
      for (i in seq_along(stoichCat)) {
        stoichMet[i] <- lambda_*stoichCat[i] + stoichAn[i]
      }
    } else {
      stoichMet <- stoichAn
    }
    delGdis0 <- sum(lambda0*-delGcat0) - delGan0
    delGdis <- sum(lambda_*-delGcat) - delGan
  }
  
  #Calculating CUE, NUE, TER ===================================================
  CUE <- stoichMet[length(stoichMet)]*1 / (abs(stoichMet[1])*a)
  NUE <- stoichMet[length(stoichMet)]*0.2 / (abs(stoichMet[1])*c + abs(stoichMet[4])*1)
  TER <- (NUE/CUE)*(1/0.2)
  
  #Assembling output data format ===============================================
  prefixes <- c("stoichD_", "stoichA_", "stoichCat_", "stoichAn_", "stoichMet_")
  suffixes <- c("donor", "h2o", "hco3", "nh4", "hpo4", "hs", "h", "e", "O2", "NO3", "NO2", "N2", "Fe3+", "Fe2+", "H2", "SO4", "H2S", "SO3", "S", "S2O3", "Fe(OH)3", "FeOOH", "Fe3O4", "CH4", "CH3COO-", "MnO2", "Mn2+", "biom")
  lengths <- c(length(stoichD), length(stoichA), length(stoichCat), length(stoichAn), length(stoichMet))
  names_list <- c("CUE", "NUE", "TER", "delGcox0", "delGd0", "delGcat0", "delGan0", "delGdis0", "lambda0", "delGcox", "delGd", "delGcat", "delGan", "delGdis", "lambda")
  for (i in 1:length(prefixes)) {
    names_list <- c(names_list, outer(paste0(prefixes[i], ""), suffixes[1:lengths[i]], paste0))
  }
  
  all_values <- c(CUE, NUE, TER, delGcox0, delGd0, delGcat0, delGan0, delGdis0, lambda0, delGcox, delGd, delGcat, delGan, delGdis, lambda_)
  all_values <- c(all_values, stoichD, stoichA, stoichCat, stoichAn, stoichMet)
  
  lambda <- data.frame(matrix(all_values, nrow = 1, byrow = TRUE), stringsAsFactors = FALSE)
  colnames(lambda) <- names_list
  return(lambda)
}

get_lambda <- function(compositions, ui, pH_Input) {
  lambda1 <- data.frame()
  for (i in 1:nrow(compositions)) {
    lambda <- getThermoStoich(unname(compositions[i, ]), ui, pH_Input)
    lambda1 <- rbind(lambda1, lambda)
  }
  return(lambda1)
}

CHEMICAL_ELEMENTS <- c("C", "H", "N", "O", "P", "S")
composition <- get_compositions(input_data)
lambda_out <- get_lambda(composition$chemical_compositions, electronAcceptor, pH)
rownames(lambda_out) <- composition$formulas
write.csv(lambda_out, file = output_file)