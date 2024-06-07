# This code allow users to set all relevant parameters in one location. 
# It will be sourced first during any computer run and parameters will be updated from these values

# Set code that will be used to name the results file and the model summary file
# It must be within the double quotes (" ") 
# and must be one that is recognised by R (no blanks, for example)
run_ID <- "_beta_tests"

# Set number of bootstap simulations to run with the disturbance scenario
nboot <- 20

# set species code. Here are the valid values: 
# BND = bottlenose dolphin
# GS = grey seal
# HP = harbour porpoise
# HS = harbour seal
# MW = minke whale
spec <- 'HP'
# minit.file <- 'MakeInitScenarioHP_DEB_beta.R'

# Set proportions of females in the population 
propfemale <- 0.5

# Set size of population that is being studied
pmean <-5000

# Specify number of females to simulate
sim_number <- 500

# decide whether daily energy intake will vary among individuals and days
stochastic_R <- TRUE
if (stochastic_R) {load("iPCoD_DEB_HP_parABC_StochasticR.RData")} else {load("iPCoD_DEB_HP_parABC_noStochasticR.RData")}

# include density dependence?
density_dependence <- FALSE
if (density_dependence) {# read sample values for slope of density dependence relationship
  vh_samples <- readRDS("iPCoD_DEB_DD_parABC_SlopeRelationship_Stochastic.rds")
  dd_slope_values <- vh_samples$preyslope
  dd_slope_values <- dd_slope_values[-which(is.na(dd_slope_values))]
}

# should disturbance duration be fixed?
fixed_disturbance_duration <- FALSE
# if it is fixed, the default is to use the iPCoD value of 6hr
# set mean duration of disturbance
if (fixed_disturbance_duration) disturbance_duration <- 6 else disturbance_duration <- 2
# set default shape parameter for Erlang distribution for use if fixed_disturbance_duration is FALSE
Erlang.shape <- 3

### SPECIFY DETAILS OF PILING OPERATIONS
# Set number of years on which piling will occur. 
# Set this to zero if there is no piling or you only want to model the effect of collisions
pile_years <- 3

# Input proportion of animals in the vulnerable components of the population
# The default is that the entire population is vulnerable
# if sum(vulnmean) < 1.0, then the remainder of the population (1-sum(vulnmean)) will not be affected by any operation
vulnmean <- c (0.3,0.6)
nvulnmean <- length(vulnmean)
rescale <- sum(vulnmean)


  # Specify name of csv file that contains information on days on which piling will occur 
  # the name should be placed between the quotation marks: " "
  piling.file <- "MultPilingOpsMultYears.csv"
  
# Set number of piling Operations to be modelled
  pilesx1 <- 3

# "vulnpile" is a matrix indicating which columns of piling.file are to be combined 
# to predict the effects of piling on each vulnerable component of the population
vulnpile <- matrix(0, nrow = nvulnmean, ncol = pilesx1)

# Indicate which Operations will affect each vulnerable component of the population
# A separate line is required for each vulnerable component of the population.
# Add vulnpile[2,] <- c(,) if there are two vulnerable components, 
# vulnpile[3,] <- c(,) if there are three vulnerable components, and so on.
# In this example, animals in the first vulnerable component are vulnerable to the effects of 
# operations 1 & 3, and animals in the second component are vulnerable to the effects of operation 2
# NB EACH PILING OPERATION SHOULD ONLY AFFECT ONE VULNERABLE LOCAL POPULATION 
# AN ERROR MESSAGE WILL BE PRODUCED IF THISIS NOT THE CASE
vulnpile[1, ] <- c(0,1,0)
vulnpile[2, ] <- c(1,0,1)
if(length(which(colSums(vulnpile)>1))){stop('inappropriate allocation of piling operations to vulnerable local populations')}

# Input number of animals predicted to experience disturbance during 1 day of piling for each Operation
# This version of iPCoD includes an option to specify that the number of animals that are disturbed
# and experience PTS varies through the year because of changes in density by specifying 
# a value of 4 for the  variable "seasons". If this is the case, separate values can be specified
# for Winter, Spring, Summer and Autumn.  
# The default is that these numbers are the same throughout the year
seasons <- 1
numDt <- numPt <- matrix(c(0),nrow=seasons,ncol=pilesx1)
if (seasons == 1){
  numDt[1,] <- c(0,0,0)   # disturbance
  numPt[1,] <- c(0,0,0)   # PTS
  } else {
    # input numbers for winter (December, January and February)
      numDt[1,] <- c(50,100,150)
      numPt[1,] <- c(1,2,3)
    # input numbers for spring (March, April and May)
      numDt[2,] <- c(100,200,300)
      numPt[2,] <- c(2,4,6)
    # input numbers for summer (June, July and August)
      numDt[3,] <- c(100,200,300)
      numPt[3,] <- c(2,4,6)
    # input numbers for autumn (September, October and November)
      numDt[4,] <- c(50,100,150)
      numPt[4,] <- c(1,2,3)
    }

# set number of initial years that will be discarded 
burnin_end <- 5
# specify number of years to simulate after burnin_end
years <- 25
years <- burnin_end + years - 1


# these files contain the initial values for the DEB simulations
load("initial_reserves&calf_ages.RData")

