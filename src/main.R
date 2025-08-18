### Script for benefit-risk analysis of maternal RSV vaccination in South Africa
### created by Ayaka Monoi 


# load the require packages
if (!require(pacman)){ #load packages
  install.packages("pacman")
}

library(pacman)
pacman::p_load(char = c("tidyverse", "here","epitools", 
                        "BayesianTools", "knitr", "ggplot2",  "MASS", 
                        "dplyr", "conflicted", "ggtext", "MCMCpack", "readxl",
                        "fitR","lattice","gridExtra","ggpubr"))

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

dir.create("output", showWarnings = FALSE)

# estimate waning protection during infancy from maternal RSV vaccination (main, including Fig1)
source("src/ve_fit.R")


# estimate RSV disease burden in South Africa (main)
source("src/RSV_burden.R")

# estimate RSV-associated infant deaths averted through vaccination (main)
source("src/benef_calc.R")

# estimate gestational age (GA)-specific neonatal mortality in South Africa (main)
source("src/nmr.R")

# look at GA-specific births in South African component of RSVpreF phase 3 trial (main, including Fig2B)
source("src/ptb_posterior.R")

# estimate excess neonatal deaths among infants born to vaccinated mothers compared to unvaccinated mothers in South Africa (main, including Fig2C)
source("src/risk_calc.R")

# look at GA-specific births born to mothers vaccinated at 27-36 weeks in South African component of RSVpreF phase 3 trial (main, including Fig3)
source("src/ptb_posterior_vacc.R")

# combine risk and benefit to calculate risk-benefit ratio (main)
source("src/benefit_risk.R")

# reproduce figures (main+supplement, including Fig2, Fig3, Fig4, FigS5, FigS6)
source("src/reprod_fig.R")

# comparison of waning efficacy models (supplement)
source("src/ve_fit_compare.R")

# include neonatal mortality of infants born at 37+ weeks from pooled LMIC data (Supplement)
source("src/nmr_sens.R")

# Sensitivity analysis for GA dating (Supplement)
source("src/nmr_ga_dating.R")
