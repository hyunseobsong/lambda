" 
Compute the lambda from chemical compositions
"
library(tidyverse)

CHEMICAL_ELEMENTS = c("C","H","N","O","P","S")

# functions ------------------------------------------------------

# extract chemical compositions from the table
get_compositions <- function(df) {
  chemical_compositions <- NULL
  formulas <- NULL
  if ("C" %in% colnames(df)) {
    tdf <- df %>% 
      filter(C > 0)
    if ("C13" %in% colnames(df)) {
      tdf <- tdf %>% 
        filter(C13 == 0)
    }
    chemical_compositions <- as.matrix(tdf[CHEMICAL_ELEMENTS])
    formulas <- tdf$MolForm
  } else if ("MolForm" %in% colnames(df)) {
    tdf <- df %>% 
      drop_na(MolForm) %>%
      filter(MolForm != "")
    parse_output <- parse_formulas(tdf$MolForm)
    formulas <- tdf$MolForm[parse_output$is_valid]
    chemical_compositions <- parse_output$composition[parse_output$is_valid,]
    warning("`MolForm` column is parsed to get the chemical compositions")
  } else {
    error("Either columns for compositions (e.g., C, H, N, ...) or `MolForm` column is required.")
  }
  
  if ("Z" %in% colnames(df)) {
    chemical_compositions <- cbind(chemical_compositions, "Z"=df$Z)
  } else {
    chemical_compositions <- cbind(chemical_compositions, "Z"=0)
  }
  
  list(
    "chemical_compositions" = chemical_compositions,
    "formulas" = formulas
  )
}

# parse formulas into compositions
parse_formulas <- function(formulas) {
  rst <- array(0, dim=c(length(formulas), length(CHEMICAL_ELEMENTS)))
  is_valid <- array(TRUE, dim=c(length(formulas)))
  for (k in 1:length(formulas)){
    formula <- formulas[k]
    ge <- gregexpr("[A-Z]\\d*", formula, perl=TRUE)
    s_index <- ge[[1]]
    s_len <- attr(s_index, "match.length")
    for (i in 1:length(s_len)){
      token <- substr(formula, s_index[i], s_index[i] + s_len[i] - 1)
      element <- substr(token, 1, 1)
      if (grepl(element, "CHNOSP")) {
        idx = which(CHEMICAL_ELEMENTS %in% element)
        if (rst[k, idx] > 0) {   # same element again? (e.g., C13)
          if (token != "C13") {  
            warning(paste0(formula,": wrong format"))
          }
          is_valid[k] = FALSE
          next
        }
        if (s_len[i] == 1) {
          rst[k, idx] = 1
        } else {
          num_element <- try(strtoi(substr(formula, s_index[i] + 1, s_index[i] + s_len[i] - 1)))
          if (class(num_element)=="integer"){
            rst[k, idx] = num_element
          } else {
            print(paste("[ERROR] an unknown chemical element found:", token, "in", formula))
          }
        }
      } else {
        print(paste("[ERROR] an unknown chemical element found:", element, "in", formula))
      }
    }
  }
  colnames(rst) <- CHEMICAL_ELEMENTS
  list("composition"=rst, "is_valid"=is_valid)
}

# compute thermodynamic properties and lambda values
getThermoStoich <- function(chemForm) {
  a <- chemForm[1]
  b <- chemForm[2]
  c <- chemForm[3]
  d <- chemForm[4]
  e <- chemForm[5]
  f <- chemForm[6]
  z <- chemForm[7]
  
  # Step 1a) stoichD: stoichiometries for an electron donor
  ySource <- -1
  yH2o <- -(3*a+4*e-d)
  yHco3 <- a
  yNh4 <- c
  yHpo4 <- e
  yHs <- f
  yH <- 5*a+b-4*c-2*d+7*e-f
  yE <- -z+4*a+b-3*c-2*d+5*e-2*f
  stoichD <- c(ySource,yH2o,yHco3,yNh4,yHpo4,yHs,yH,yE)
  stoichD[c(9,10)] <- 0 # add additional components: e-acceptor and biomass
  
  # Step 1b) stoichA: stoichiometries for an electron acceptor (i.e., oxygen)
  stoichA <- rep(0, 10)
  stoichA[9] <- -1  # oxygen
  stoichA[7] <- -4  #  h+
  stoichA[8] <- -4  #  e-
  stoichA[2] <- 2  #  h2o
  
  # Step 1c) stoichCat: stoichiometries for catabolic reaciton 
  yEd <- stoichD[8]
  yEa <- stoichA[8]
  stoichCat <- stoichD-(yEd/yEa)*stoichA
  
  # Step 2a) stoichAnStar: stoichiometries for anabolic reaciton 
  #          (N source = NH4+)
  chemFormBiom <- c(1, 1.8, 0.2, 0.5, 0, 0, 0)  # C H_1.8 N_0.2 O_0.5
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
  yE <- -zB+4*aB+bB-3*cB-2*dB+5*eB-2*fB
  stoichAnStarB <- c(ySource,yH2o,yHco3,yNh4,yHpo4,yHs,yH,yE)
  stoichAnStarB[c(9,10)] <- 0  # add additional components: e-acceptor and biomass
  stoichAnStarB <- -stoichAnStarB
  stoichAnStarB[10] <- stoichAnStarB[1]
  stoichAnStarB[1] <- 0
  
  # Step 2b) "overall" anabolic reaction
  stoichAnStar <- stoichAnStarB+(1/a)*stoichD
  yEana <- stoichAnStar[8]
  if (yEana > 0)
    stoichAn <- stoichAnStar-yEana/yEa*stoichA
  else if (yEana < 0)
    stoichAn <- stoichAnStar-yEana/yEd*stoichD
  else
    stoichAn <- stoichAnStar
  
  # Step 3: get lambda
  
  # - estimate delGd0 using LaRowe and Van Cappellen (2011)
  ne <- -z+4*a+b-3*c-2*d+5*e-2*f  # number of electrons transferred in D 
  nosc <- -ne/a+4  # nominal oxidataion state of carbon 
  delGcox0 <- 60.3-28.5*nosc  # kJ/C-mol
  delGd0 <- delGcox0*a*abs(stoichD[1])  # kJ/rxn
  
  # - estimate delGf0 for electron donor
  delGf0_D_zero <- 0
  # delGf0_zero <- c(delGf0_D_zero, -237.2, -586.8, -79.4, -1096.1, 12.1, 0, 0, 16.4, -67)
  delGf0_zero <- c(delGf0_D_zero, -237.2, -586.9, -79.5, -1089.1, 12.0, 0, 0, 16.5, -67)
  delGcox0_zero <- drop(delGf0_zero %*% stoichD)
  delGf0_D_est <- (delGd0-delGcox0_zero)/stoichD[1]
  # - finally, delGf0
  delGf0 <- delGf0_zero
  delGf0[1] <- delGf0_D_est
  
  # - standard delG at pH=0
  delGcat0 <- drop(delGf0 %*% stoichCat)
  delGan0 <- drop(delGf0 %*% stoichAn)
  
  # - stadard delG at pH=7
  R <- 0.008314  # kJ/(K.mol)
  T <- 298.15  # K
  iProton <- 7  # [eD,h2o,hco3-,nh4+,hpo4^2-,hs-,h+,e-,eA,biom]
  delGd <- delGd0+R*T*stoichD[iProton]*log(1e-7)
  delGcox <- delGd / a
  delGcat <- delGcat0+R*T*stoichCat[iProton]*log(1e-7)
  delGan <- delGan0+R*T*stoichAn[iProton]*log(1e-7)
  
  # The Thermodynamic Electron Equivalents Model (TEEM)
  # --------
  eta <- 0.43
  delGsyn <- 200  # kJ/(mol.X)
  if (is.nan(delGan0) & is.nan(delGan)) {
    lambda0 <- NaN
    lambda <- NaN
    stoichMet <- array(NaN, dim=length(stoichCat))
    delGdis0 <- NaN
    delGdis <- NaN
  } else {
    if (delGan < 0)
      m <- 1
    else
      m <- -1
    
    if (delGan0 < 0)
      m0 <- 1
    else
      m0 <- -1
    
    lambda0 <- (delGan0*eta^m0+delGsyn)/(-delGcat0*eta)
    lambda <- (delGan*eta^m+delGsyn)/(-delGcat*eta)
    
    if (lambda > 0)
      stoichMet <- lambda*stoichCat+stoichAn
    else
      stoichMet <- stoichAn
    
    delGdis0 <- drop(lambda0%*%(-delGcat0)) - delGan0
    delGdis <- drop(lambda%*%(-delGcat)) - delGan
  }
  
  c(delGcox0,delGd0,delGcat0,delGan0,delGdis0,lambda0,
    delGcox,delGd,delGcat,delGan,delGdis,lambda,
    stoichD,stoichA,stoichCat,stoichAn,stoichMet)
}

# compute in batch
get_lambda <- function(formula_matrix) {
  nrows = nrow(formula_matrix)
  lambda_rst <- array(0, dim=c(nrows, 62))
  for(i in 1:nrows) {
    lambda_rst[i,] <- getThermoStoich(formula_matrix[i,])
  }
  lambda_rst
}


# user parameters ------------------------------------------------------

outfile <- "demo_input_out.txt"
fticr_data <- read_csv("demo_input.csv")


# main run -------------------------------------------------------------

info <- get_compositions(fticr_data)
out <- get_lambda(info$chemical_compositions)

# build data frame
df <- as.data.frame(out)

# build col names
names <- rep("", 62)

names[1:12] <- c("delGcox0","delGd0","delGcat0","delGan0","delGdis0","lambda0",
                 "delGcox","delGd","delGcat","delGan","delGdis","lambda")

stoich_colnames <- c("donor","h2o","hco3","nh4","hpo4","hs","h","e","acceptor","biom")
stoich_types <- c("stoichD","stoichA","stoichCat","stoichAn","stoichMet")

for (i in 1:length(stoich_types)) {
  names[((i-1)*10+13):(i*10+12)] <- array(sapply(stoich_types[i], paste, stoich_colnames, sep="_"))
}
colnames(df) <- names
df['MolForm'] <- info$formulas

write.table(df, file = outfile, row.names=FALSE, sep = "\t")
