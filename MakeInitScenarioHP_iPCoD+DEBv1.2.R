# Load in the appropriate values from the results of the expert elicitation on PTS
load ("outputList.robj")
matureSurvival_2_10kHz <- outputList$HP_MF_surv_6dB
# adjust entry for matureSurvival_2_10kHz$density[1] and matureReproduction 
# so that they aren't infinity!
matureSurvival_2_10kHz$density[1] <- 250
matureReproduction_2_10kHz <- outputList$HP_fert_6dB
juvenileSurvival_2_10kHz <- outputList$HP_depjuv_6dB
dependentSurvival_2_10kHz <- outputList$HP_depjuv_6dB

## SET DEMOGRAPHIC RATES TO THOSE RECOMMENDED BY SINCLAIR ET AL FOR iPCoD
# Set juvenile survival
J_surv <- 0.85
# Set age at which an average female gives birth to her first calf
age2 <- 5
# Set adult survival
Ad_surv <- 0.931

# Set calf/pup survival
Calf_surv <- 0.8455
# and calf age at independence
age1 <- 1
Fertility <- 0.34*propfemale

age_classes <- 40

# determine initial stable age structure for population from Leslie matrix
L                  <- array(0, dim = c(age_classes, age_classes))
L[1, age2]          <- Fertility * J_surv
L[1, (age2 + 1):age_classes] <- Fertility * Ad_surv
index3             <- array(c(2:(age1+1), 1:age1), dim = c( age1 ,2) )
L[index3]            <- Calf_surv
mature             <- ifelse(age2 == 9, 9, age2 + 1)
index1             <- array(c((mature + 1):age_classes, (mature):(age_classes-1)), dim = c((age_classes - mature), 2))
L[index1]          <- Ad_surv
L[age_classes,age_classes] <- ifelse(age_classes==10,Ad_surv,0)
index2             <- array(c((age1+2):(age2 + 1), (age1+1):age2), dim = c((age2 - age1), 2))
L[index2]          <- J_surv

ev              <- eigen(L)
population_growth_rate <- ev$val[1]
# age_structure holds the proportion of animals in each age class (excluding calves) for a stable age structure  
age_structure   <- as.numeric(ev$vec[2:age_classes, 1] / (sum(ev$vec[2:age_classes, 1])))

## Generate age structure for DEB that results in the same mortality rate for 9+ females as the original iPCoD model
DEB_adult_survival <- Ad_surv
#DEB_age_structure <- sapply(0:30, function(j) DEB_adult_survival^j)
#DEB_age_structure <- age_structure[9]*DEB_age_structure/sum(DEB_age_structure)
#age_structure <- age_structure[-9]
#age_structure <- c(age_structure,DEB_age_structure)

save(age_structure, Calf_surv, J_surv, Ad_surv, DEB_adult_survival, Fertility, age2, population_growth_rate, file = "initial_demographics.RData")

# newvulnmean includes proportion of animals in undisturbed remainder of population, if there is one!
  newvulnmean <- vulnmean
  if(sum(vulnmean)== 1){newvulnmean <- newvulnmean} else {newvulnmean[nvulnmean+1] <- 1 - sum(newvulnmean)}

# if pile_years is zero, there is no need to process "pile" 
if (pile_years > 0) {
  
# read csv file with schedule of piling activities
pile <- read.csv(file = piling.file, header = TRUE)
  
# remove day labels from first column of csv file
library(stringr)
yvec   <- str_match(colnames(pile), 'Operation')
npiles <- length(yvec[!is.na(yvec)])
pilesx <- which(str_match(colnames(pile), 'Operation') %in% 'Operation') # remove any unwanted columns, i.e. row names that come over from excel
pile  <- pile[, pilesx, drop = FALSE]
pilesx <- which(str_match(colnames(pile), 'Operation') %in% 'Operation') # repeat to get the vector dimensions correct.

# check that number of piling operations specified in config matches number in piling file
if(pilesx1 != ncol(pile)){stop('Number of Piling Operations do not match')}

# check that number of days in piling file is an exact multiple of 365
if(nrow(pile)/365 != pile_years){stop('Number of days in piling file is not an exact multiple of 365')}

# add extra days to beginning of piling schedule so that animals are exposed to disturbance from 
# 1 June in year before piling starts to 31 May in first year of piling
for (i_extra in 1:214){pile <- rbind(rep(c(0),npiles),pile)}
# add extra days to end of piling schedule so that animals are exposed to disturbance from 
# 1 June in last year of piling to 31 May in year after piling ends
for (i_extra in 1:151){pile <- rbind(pile,rep(c(0),npiles))}

# template for  piling operations should ensure that number of rows is an exact multiple of 365
pile_years <- nrow(pile) / 365 
# yearvec indicates year number for each day in pile
yearvec <- rep(1:(pile_years), each = 365)

pile <- data.frame(pile, pbdays = rowSums(pile), pvec = rowSums(pile))
# pile$pbdays indicates number of piling events on a particular day, pile$pvec is a flag indicating whether or not piling has occurred
pile$pvec[pile$pvec > 1] <- 1 

# NDt will contains number of animals predicted to be disturbed during 1 day of activity for each piling operation
# NPt will contains number of animals predicted to be experience PTS during 1 day of activity for each piling operation
# values are repeated for each day in "pile"
if (seasons == 1){
  NDt <- matrix(rep(numDt[1,1:npiles],nrow(pile)), nrow(pile), npiles, byrow = TRUE)
  NPt <- matrix(rep(numPt[1,1:npiles],c(365)),nrow(pile), npiles, byrow = TRUE)
  } else {
  day_seq <- c(92,91,90,92)
  cum_day_seq <- c(1,93,184,274)
  N_rows <- c(3,4,1,2)
  daily_NDt <- daily_NPt <- matrix(c(0),365,npiles)
  for (i_piles in 1:npiles){
    daily_NDt[,i_piles] <- unlist(sapply(1:4, function(X) rep((numDt[N_rows[X],i_piles]),day_seq[N_rows[X]])))
    daily_NPt[,i_piles] <- unlist(sapply(1:4, function(X) rep((numPt[N_rows[X],i_piles]),day_seq[N_rows[X]])))
  }
  ## duplicate daily_ matrices so that they cover the entire time period when piling can occur
  NDt <- daily_NDt
  NPt <- daily_NPt
  if (pile_years > 1) {
    for (i_annual in 2:pile_years){
      NDt <- rbind(NDt,daily_NDt)
      NPt <- rbind(NPt,daily_NPt)}
  }
  # end of calculations if seasons > 1  
}

#creates a column for each vulnerable sub-population indicating the days on which there is piling that will affect it
for(i in 1:nrow(vulnpile)){
  pvec <- as.matrix(pile[, pilesx]) %*% as.matrix(vulnpile[i, ])
  pvec[pvec > 1] <- 1
  colnames(pvec) <- paste('vuln', i, 'pvec', sep = '')
  pile <- cbind(pile, pvec)              
}


  # p_disturb_base and p_pts_base contain baseline probabilities for disturbance and PTS for each vulnerable sub-pop
  # on each day. They are multiplied by the value of sdNDt for each bootstrap to get the values for each set of simulations
  p_disturb_base  <- matrix(0, nrow = pile_years * 365, ncol = nvulnmean)
  colnames(p_disturb_base) <- paste("sub-pop", 1:nvulnmean, sep = '')
  p_pts_base <- p_disturb_base
  
  for (j in 1:(nvulnmean)) { 
    for(k in 1:nrow(p_disturb_base))  {
      p_disturb_base[k, j] <- ifelse( sum((pile[k, pilesx] * vulnpile[j, pilesx] * NDt[k, pilesx]) / (vulnmean[j] * pmean)) > 1.0, 1.0, 
                                      sum((pile[k, pilesx] * vulnpile[j, pilesx] * NDt[k, pilesx]) / (vulnmean[j] * pmean))) 
      
      p_pts_base[k, j]     <- sum(pile[k, pilesx] * vulnpile[j, pilesx] * NPt[k, pilesx]) / sum(pile[k, pilesx] * vulnpile[j, pilesx] * NDt[k, pilesx])
    }
  }
  
  p_pts_base[is.na(p_pts_base)] <- 0  
  
}

  
  ### convert p_disturb_base and p_pts_base into an array with as many columns as years
  ## for each vulnerable population
  p_disturb_by_year <- p_pts_by_year <- array(c(0),dim = c(365,pile_years,nvulnmean))
  p_disturb_by_year[,,1:nvulnmean]<- p_disturb_base
  p_pts_by_year[,,1:nvulnmean] <- p_pts_base
  
    
  # matrix to store values that will be retained from all simulations
  # age.out stores full age structure for disturbed and undisturbed populations at start of breeding season age.out[,1&2,]
  
  dat.out <- array(NA, dim = c((years+1), 5, nboot))

  # Ndist is the total population size for the disturbed population each year, NNotDist is the equivalent for the matching undisturbed population
  colnames(dat.out) <- c('Impacted', 'Unimpacted','unimpacted-impacted','average pa% change Impacted','average pa% change Unimpacted')
  
  # population_size and population_size_d contain detailed information about what happens in the simulated population
  population_size <- array(c(0),dim = c((years+1), 12, nboot))
  pop.names <- c("calves born","juv f","mature f","1+ f","preg rate","birth rate","calf surv","mean days disturb","new PTS", "PTS prop","starve deaths","Rsim")
  colnames(population_size) <- pop.names
  population_size_d <- population_size
  


