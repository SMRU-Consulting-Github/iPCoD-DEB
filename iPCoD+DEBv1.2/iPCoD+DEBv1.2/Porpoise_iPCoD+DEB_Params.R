

##################################################
###### SET KEY PARAMETERS
# set maximum age an individual can attain
max_age <- 40
# calculate maximum age in days 
maximum_age <- max_age*365
# set maximum duration of simulation
maximum_days <- years*365
# "years" is number of years for simulation from config 


# set assumed reserve level of calves at weaning 
rho_j <- 0.3243

# steepness of assimilation response
#eta <- 15
eta <- 20
# relative cost of maintaining reserves
Theta_F <- 0.2
# shape parameter for starvation-induced reduction in milk delivery
# xi_m <- -2
xi_m <- -3
# set length of gestation
Tp <- 305
#Tp <- 335

# set foetal mortality.  This value is based on Murphy et al 2015 
foetal_mortality <- 0.197

# Kleiber constant for estimating field maintenance costs
K <- 0.294
# length of lactation. Assumed to be 8 months
Tl <- 250

##################################################
# create vector of dates coded as Julian days, starting on 1 June, no leap years
birthday <- 152
age_day1 <- 
day1 <- birthday
# create vector of dates coded as Julian days, starting on 1 June + Tl (i.e. end of lactation period) no leap years
# day1 <- birthday + Tl - 365

julian_days <- c(day1:365)
remaining <- length(julian_days) + 1
julian_days[remaining:365] <- c(1:(day1-1))
# replicate the Julian days vectors
max_years <- ceiling(maximum_days/365) + 1
julian_days <- rep(julian_days,max_years)

# variations in target body condition over year

#rho_max <- 0.45
rho_max <- 0.35
#rho_min <- 0.3
rho_min <- 0.24
rho_annual <- rep(rho_max, 365)
# steady decline in rho from 1 March (Julian 60) to 1 July (Julian 182)
rho_annual[60:182] <- rho_max - (rho_max-rho_min)*c(1:123)/(182-59)
# steady rise in rho from 2 July (Julian 183) to 1 December (Julian 335)
#rho_annual[183:335] <- rho_min + (rho_max-rho_min)*c(1:153)/(335-182)
####### 9 March modification
# rho_min remains low until 15 September (Julian 258) then rises to 1 December (Julian 335)
rho_annual[183:258] <- rho_min
rho_annual[259:335] <- rho_min + (rho_max-rho_min)*c(1:77)/(335-258)

# annual cycle in catabolic and anabolic reserve conversion efficiencies in MJ/kg
epsi_minus_annual <- rep(20,365)
epsi_minus_annual[60:182] <- 35
epsi_plus_annual <- rep(55,365)
epsi_plus_annual[60:182] <- 28

# rescale rho and epsi's so that they match females time scale and calculate rho_s
rho_t <- rho_annual[julian_days[1:365]]
epsi_minus <- epsi_minus_annual[julian_days[1:365]]
epsi_plus <- epsi_plus_annual[julian_days[1:365]]

rho_t <- rep(rho_t,max_years)
epsi_minus <- rep(epsi_minus, max_years)
epsi_plus <- rep(epsi_plus, max_years)
# calculate starvation threshold
#rho_s <- rho_t - 0.1
rho_s <- rho_t - 0.11
#########################################################

#######################################################
######### GROWTH PARAMETERS FOR FEMALE ################
# length at infinity 
Linf <- 160
# length at birth 
Lb <- 70
# growth rate 
k <- 0.0015
# mass length scaling 
omega1 <- 5.9*10^(-5)
# mass length scaling exponent
omega2 <- 2.67

# energy efficiency per unit structural mass 
Sigma_G <- 30 

############## calculate length at age
ages <- c(1:(maximum_age+1))
La <- Linf - (Linf-Lb)*exp(-k*ages)
# calculate core weight at age 
Sa <- omega1*La^omega2
######### calculate growth costs for all age classes
CGa <- Sigma_G*(Sa[2] - Sa[1])
CGa[2:maximum_age] <- sapply(2:maximum_age, function(X) Sigma_G*(Sa[X+1] - Sa[X]))

#########################################################
############## LACTATION PARAMETERS #####################
# shape parameter for decrease in milk supply with calf age
xi_c <- 0.5

### calculate change in % of calf's demand delivered by female over course of lactation
# extended beyond Tl to allow modelling of calf's growth post weaning if required
milk_days <- c(1:(Tl+365))
temp_prop <- (1 - (milk_days-Tn)/(Tl-Tn))/(1 - xi_c*(milk_days-Tn)/(Tl-Tn))
milk_prop <- ifelse(temp_prop>1,1,temp_prop)
milk_prop[(Tl+1):(Tl+365)] <- 0

# milk equivalent of R
phi_L <- 3.55
# efficiency of conversion of mother's reserves to milk
Sigma_L <- 0.86

# calculate feeding efficiency on each day over all possible ages
assim_effic <- sapply(1:maximum_age, function(X) X^upsilon/(Tr^upsilon + X^upsilon))



## CALCULATE AGE-SPECIFIC SURVIVAL RATES ####################################################
# calculate survival rates for foetus
daily_foetal_surv <- (1-foetal_mortality)^(1/Tp)
cum_Surv_foetus <- sapply(1:maximum_age,function(i) daily_foetal_surv^i)
cum_Surv_foetus[maximum_age] <- 0

# age-specific cumulative survival loosely based on Sinclair et al. (2020) 
# but adult survival modified to allow multiple age classes (up to max_age) after age 10
# see initial_age_structure.R, which generates "initial_demographics.RData" for more details
# that file also contains the stable age structure (age_structure) associated with the demographic rates
# specified in the file and the population growth rate

load("initial_demographics.RData")
# reset Calf_surv so that population rate of increase without any starvation mortality is positive
#Calf_surv <- 0.79
#J_surv <- 0.85
#Ad_surv <- 0.93
Ad_surv <- DEB_adult_survival

index <-(1:maximum_age)
daily_surv <- 0
daily_surv[1:365] <- Calf_surv^(1/365)
daily_surv[366:(age2*365)] <- J_surv^(1/365)
daily_surv[(age2*365+1):maximum_age] <- Ad_surv^(1/365)

cum_Surv <- 1
cum_Surv <- cumprod(daily_surv)
cum_Surv[maximum_age] <- 0

#### Virtual experts' predictions of effect of PTS from posteriors
# 1. survival as a result of PTS at 2-10 kHz Mature Females
MFsurv_2_10_probs <- matureSurvival_2_10kHz$density
MFsurv_2_10_X <- matureSurvival_2_10kHz$X/100
PTS.Ad_surv <- 1 - sample(MFsurv_2_10_X, 1, MFsurv_2_10_probs, replace = T)

# 3. Juveniles survival as a result of PTS at 2-10 kHz
Jsurv_2_10_probs <- juvenileSurvival_2_10kHz$density
Jsurv_2_10_X <- juvenileSurvival_2_10kHz$X/100
PTS.J_surv <- 1 - sample(Jsurv_2_10_X, 1, Jsurv_2_10_probs, replace = T)
 
# Dependents survival as a result of PTS at 2-10 kHz
Dsurv_2_10_probs <- dependentSurvival_2_10kHz$density
Dsurv_2_10_X <- dependentSurvival_2_10kHz$X/100
PTS.Calf_surv <- 1 - sample(Dsurv_2_10_X, 1, Dsurv_2_10_probs, replace = T)

PTS.daily_surv <- 0
PTS.daily_surv[1:365] <- (Calf_surv*PTS.Calf_surv)^(1/365)
PTS.daily_surv[366:(age2*365)] <- (J_surv*PTS.J_surv)^(1/365)
PTS.daily_surv[(age2*365+1):maximum_age] <- (Ad_surv*PTS.Ad_surv)^(1/365)

PTS.cum_Surv <- 1
PTS.cum_Surv <- cumprod(PTS.daily_surv)
PTS.cum_Surv[maximum_age] <- 0

##################### PREGNANCY PARAMETERS #############
# age at sexual maturity, this calculation ensures that first birth occurs on the age2th birthday
min_age <- age2 - 1
min_age <- min_age*365


# if rho is above threshold for pregnancy on the day of the year when pregnancy is possible, implantation can occur
preg_possible <- rep(c(0),maximum_days+1)
implants <- seq((365-Tp+1),maximum_days+1,365)
preg_possible[implants] <- 1

# probability that implantation will occur
#fert_success <- 0.95
fert_success <- 1

#### Virtual experts' predictions of effect of PTS from posteriors
# 2. fecundity as a result of PTS at 2-10 kHz
MFrepro_2_10_probs <- matureReproduction_2_10kHz$density
MFrepro_2_10_X <- matureReproduction_2_10kHz$X/100
#PTS.fert_success <- 1 - sample(MFrepro_2_10_X, 1, MFrepro_2_10_probs, replace = T)
PTS.fert_success <- 0.9

###### calculate length and weight at age for foetus (assumes linear growth in length)
gest_t <- c(1:(Tp+1))
# NB L_foetus[Tp+1] is calculated solely for the purpose of determining growth costs
# L_foetus[Tp] is equivalent to La[0], i.e. length at birth in the cB growth curve
L_foetus <- gest_t*Lb/Tp
S_foetus <- omega1*L_foetus^omega2


first_days <- seq(1,(maximum_age-1), 365)

f.foetusandcalf_life <- function(i){
  # determine foetal and calf life expectancy
  # with these values, ~20% of calves die from causes other than starvation before weaning
  check1 <- runif(1)
  check2 <- runif(1,cum_Surv[maximum_age-1],1)
  life_expect_f <- min(which(cum_Surv_foetus<=check1))
#  if (life_expect_f <= Tp) {life_expect_c <- 0} else
  if (life_expect_f < Tp) {life_expect_c <- 0} else
    {life_expect_c <- min(which(cum_Surv<=check2))
     life_expect_f <- Tp}
  # calf growth and energy expenditure is modelled up to age 365 days
  # return values
  list (life_expect_f = life_expect_f, life_expect_c = life_expect_c)
  }

f.starve_death <- function(rhoC,rho_s){
  # determine whether a simulated animal dies (return 0) as a result of starvation on a particular day
  mort <- mu_s*(rho_s/rhoC - 1)
  rbinom(1,1,exp(-mort))
  }

##### parameters for beta distribution if modelling stochasticity in resources
mu <- 0.25
SD <- 0.055
b <- mu*(1-mu)^2/SD^2 + mu - 1
a <- mu*b/(1-mu)

