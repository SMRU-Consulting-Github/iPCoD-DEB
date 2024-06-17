

# SET UP RELEVANT MATRICES 
# F_sim is total reserves, F_preg is threshold level of reserves to become pregnant, 
Sa_sim <- F_sim <- rho_sim <- Ir_sim <- Cm_sim <- CGest_sim <- starv_prop <- matrix(c(0),nrow = (maximum_days+1),ncol = sim_number)
age_days <- F_preg <- matrix(c(0),nrow = (maximum_days+1),ncol = sim_number)
# summary statistics for each simulated female
life_expectancy <- rep(c(0),sim_number)
# R_sim.mean contains the mean simulated value for Rmean in each year
R_sim.mean <- 0

# set up matrices for calf related variables,
# each column contains data on all calves produced by a female with the equivalent columm number 
# e.g. F_C[,5] contains the reserve weights of every calf produced by the female represented by Sa_sim[,5]
F_C <- rho_C <- Ir_C <- calf_demand <- Im_C <- Cm_C <- calf_age <- matrix(c(0), nrow = (maximum_days+1), ncol = sim_number)

# matrices to store values for one-off events
# conceive_record also indicates day on which each female actually gives birth
conceive_record <- calf_deaths <- calf_days_record <- calf_life_expectancy <- matrix(c(0),nrow = (years+1), ncol = sim_number)

########## caclulate pregnancy costs used in determining threshold for pregnancy
# calculate maintenance cost of foetal tissue for a 10 year old female
CM_foetus <- sapply(1:Tp, function(X) K*Sigma_M*((Sa[(10*365+X)]+S_foetus[X])^0.75 - (Sa[(10*365+X)])^0.75))
# calculate growth cost of foetus
CG_foetus <- 0
CG_foetus[1] <- Sigma_G*S_foetus[2]
CG_foetus[2:Tp] <- sapply(2:Tp, function(X) Sigma_G*(S_foetus[X+1] - S_foetus[X]))
# calculate calf's total reserves at birth = minimum value of rho_s
#F_start <- rho_s[365-age_day1]*S_foetus[Tp+1]/(1 - rho_s[365-age_day1])
F_start <- (rho_min - 0.11)*S_foetus[Tp+1]/(1 - (rho_min - 0.11))

# calculate daily and total cost of gestation
Cost_gestation <- sum(CG_foetus[1:Tp]) + sum(CM_foetus[1:Tp]) + max(epsi_plus)*F_start
# calculate extra cost of adding F_start and spread this over entire gestation period
extra_fat_cost <- Cost_gestation/(Cost_gestation - max(epsi_plus)*F_start)
CGest <- (CG_foetus + CM_foetus)*extra_fat_cost

# calculate options for threshold reserves for start of pregnancy
#F_neo <- 0.6*sum(CG_foetus[1:Tp])/max(epsi_minus)
#F_neo <- sum(CG_foetus[1:Tp])/max(epsi_minus)
F_neo <- (sum(CG_foetus[1:Tp])+max(epsi_plus)*F_start)/max(epsi_minus)
#F_neo <- Cost_gestation/max(epsi_minus)

#### start with stable age structure
initial_age_structure <- rmultinom(1,sim_number,age_structure)

##### THIS CODE DETERMINES THE INITIAL AGES, BODY CONDITION AND REPRODUCTIVE STATUS OF FEMALE TO BE SIMULATED
#### based on values from a simulation run of 1000 animals with one_age_class set to TRUE
##### so all simulated animals survive to maximum_age (14600 days)
##### initial_reserves&calf_ages.RData contains reserve levels calf_life_expectancy and calf_age at start of each year of life for those 1000 females
#save(initial_F, initial_calf_life_expectancy, initial_calf_age, file = "initial_reserves&calf_ages.RData")
#load("initial_reserves&calf_ages.RData")

begin <- 1
for (i in 1:length(initial_age_structure)){
  if (initial_age_structure[i]>0)
  {end <- begin + initial_age_structure[i] - 1
# age_days [1,] will contain age (in days) of each female at the start of the simulation
  age_days[1,begin:end] <- rep(i,initial_age_structure[i])*365 + 1
  sampled_females <- sample(c(1:ncol(initial_F)),initial_age_structure[i],replace=FALSE)
# the amount of reserves associated with each of these females will go into F_sim[1,]
  F_sim[1,begin:end] <- initial_F[i,sampled_females]
# calf_age[1,] == 1 indicates which females will give birth to a calf on day 1 of the simulation
  calf_age[1,begin:end] <- initial_calf_age[i,sampled_females]
  calf_life_expectancy[1,begin:end] <- initial_calf_life_expectancy[i,sampled_females]
  begin <- end + 1}
}

#### determine life expectancy of each simulated female and fill appropriate rows of age_days with her age in days on those days
i_day <- 1
i_years <- 1

for (i_life in 1:sim_number){
  check <- runif(1,cum_Surv[maximum_age-1],cum_Surv[age_days[1,i_life]])
life_expectancy[i_life] <- min(which(cum_Surv<=check))
ages_in_days <- c(age_days[1,i_life]:life_expectancy[i_life])
# check if simulated female is predicted to be still alive at the end of the simulation
last_day <- i_day + life_expectancy[i_life] - age_days[1,i_life]
# if so, only simulate the days remaining in the simulation
if (last_day>nrow(age_days)) last_day <- nrow(age_days)
age_days[(i_day:last_day),i_life] <- ages_in_days[1:length(i_day:last_day)]
# fill in the structural mass of individual i_life on those days
Sa_sim[(i_day:last_day),i_life] <- Sa[ages_in_days[1:length(i_day:last_day)]]

# now fill appropriate rows of calf_age and Sa_C
if (calf_age[1,i_life]==1){
  if (calf_life_expectancy[i_years, i_life] < 365) {calf_days <- calf_life_expectancy[i_years, i_life]} else {calf_days <- 365}
  calf_ages_in_days <- c(1:calf_days)
  last_day <- i_day + calf_days - 1
  # check if calf is predicted to be still alive at the end of the simulation
  if (last_day>nrow(calf_age)) last_day <- nrow(calf_age)
  calf_age[(i_day:last_day),i_life] <- calf_ages_in_days[1:length(i_day:last_day)]
#  Sa_C[(i_day:last_day),i_life] <- Sa[calf_ages_in_days[1:length(i_day:last_day)]]
  }
}

##### SIMULATE OVER NUMBER OF YEARS REQUIRED
begin_year <- seq(1,maximum_days,365)
end_year <- seq(365,maximum_days,365)
implant_days <- which(preg_possible > 0)

# LOOP OVER YEARS
for (i_years in 1:years){
  
  ### create burned in matrices for use by PopDyn_DEB_disturbed_beta.R if burn in is complete
  if (i_years==burnin_end)  {
    # new adult matrices
    Sa_sim_d <- Sa_sim
    F_sim_d <- F_sim
    rho_sim_d <- rho_sim
    Ir_sim_d <- Ir_sim
    Cm_sim_d <- Cm_sim
    CGest_sim_d <- CGest_sim
    starv_prop_d <- starv_prop
    age_days_d <-age_days
    F_preg_d <- F_preg
    # new calf matrices
 #   Sa_C_d <- Sa_C
    F_C_d <- F_C
    rho_C_d <- rho_C
    Ir_C_d <- Ir_C
    calf_demand_d <- calf_demand
    Im_C_d <- Im_C
    Cm_C_d<- Cm_C
    calf_age_d <- calf_age
    calf_life_expectancy_d <- calf_life_expectancy
    calf_deaths_d <- calf_deaths

    
    #### copy vectors in their entirity 
    conceive_record_d <- conceive_record
    life_expectancy_d <- life_expectancy
    
    #### copy population size details
    population_size_d[,,boot_number] <- population_size[,,boot_number]
  }
  
  ### determine population size on the first day of simulated year
  alives <- which(age_days[begin_year[i_years],]>0)
  matures <- length(which(age_days[begin_year[i_years],]>min_age))
  population_size[i_years,4,boot_number] <- length(alives)
  calves_born <- length(which(calf_age[begin_year[i_years],]==1))
  population_size[i_years,1,boot_number] <- calves_born
  population_size[i_years,6,boot_number] <- calves_born/matures
  ########## CHECK THAT THERE ARE STILL SOME FEMALES ALIVE ######
  if (length(alives)>0) {
  Rmean_dd <- Rmean*(dd_slope*length(alives)/sim_number + dd_intercept)
    #### set up and fill the R_sim matrix (further modification will be necessary to model disturbance)
  R_sim <- matrix(Rmean_dd, nrow = 365, ncol = length(alives))
  # if stochastic_R = TRUE fill annual R_sim matrix with random values
  if (stochastic_R) R_sim <- sapply(1:length(alives), function(j) Rmean_dd*rbeta(365,a,b)/mu)
  R_sim.mean[i_years] <- mean(R_sim)
  population_size[i_years,12,boot_number] <- R_sim.mean[i_years]
  
  day_index <- 0
  # LOOP OVER DAYS IN YEAR i_years
  
  for (i_day in begin_year[i_years]:end_year[i_years]){
    day_index <- day_index+1
    
    # LOOP OVER FEMALES THAT ARE ALIVE AT START OF i_years
    for (i_sim in 1:length(alives)){
      
      #### CHECK THAT FEMALE IS STILL ALIVER
      if (age_days[i_day,alives[i_sim]] > 0){
        
        #### begin calculations for females that are still alive
        
      # check if female can become pregnant 
        if (age_days[i_day,alives[i_sim]] > min_age & preg_possible[i_day] > 0){
      F_preg[i_day,alives[i_sim]] <- F_neo + rho_s[i_day]*Sa[age_days[i_day,alives[i_sim]]]/(1-rho_s[i_day])
      if (F_sim[i_day,alives[i_sim]] >= F_preg[i_day,alives[i_sim]]){
        
        ##### START OF LOOP FOR FEMALES THAT ARE ABLE TO BECOME PREGNANT ON i_day ############ 
        # reset pregnancy status to 0 for rest of fertile period
        # preg_state[(i_day+1):(i_day+implant_period),i_sim] <- 0
        # check if ovum implants successfully
        if (rbinom(1,1,fert_success)==1){
          
          ###### start of loop for females that implant successfully
          tmp  <- f.foetusandcalf_life(i_sim) 
          conceive <- i_day
          foetal_life <- tmp$life_expect_f
          # store foetal life for simulated female
          conceive_record[i_years,alives[i_sim]] <- foetal_life
          calf_life_expectancy[(i_years+1), alives[i_sim]] <- tmp$life_expect_c
          
          # "birth_day" is day on which female gives birth OR day on which foetus dies
          birth_day <- conceive + foetal_life
          ####### check that female survives until birth_day
          if (age_days[birth_day,alives[i_sim]]==0) calf_life_expectancy[(i_years+1), alives[i_sim]] <- 0
          ####### end of check on female survival
            

          ######## check if there is enough time for female to give birth before the simulation ends
          if (birth_day > maximum_days){
            CGest_sim[conceive:maximum_days, alives[i_sim]] <- CGest[1:length(conceive:maximum_days)]} else 
              {
            CGest_sim[conceive:(birth_day-1),alives[i_sim]] <- CGest[1:foetal_life]
          
            ######## check if foetus or female dies before the end of the pregnancy (in which case calf life expectancy is 0)
            if (calf_life_expectancy[(i_years+1), alives[i_sim]]>0){
            ######### determine how many days of the calf's life will be simulated  
            if (calf_life_expectancy[(i_years+1), alives[i_sim]] < 365) {calf_days <- calf_life_expectancy[(i_years+1), alives[i_sim]]} else {calf_days <- 365}
              if ((birth_day+calf_days) > maximum_days+1) calf_days <- maximum_days - birth_day + 1
              # insert calf ages, calf core mass and calf growth costs 
    #          calf_days_record[i_years+1,alives[i_sim]] <- calf_days
              calf_age[birth_day:(birth_day+calf_days-1),alives[i_sim]] <- c(1:calf_days)
   #           Sa_C[birth_day:(birth_day+calf_days-1),alives[i_sim]] <- Sa[1:calf_days]
              #CGa_C[birth_day:(birth_day+calf_days-1),alives[i_sim]] <- CGa[1:calf_days]
              ######### end of calculations for foetuses that are actually born
              ############################################################### 
            }
          }
          
          ###### end of calculations for females that implant successfully
          }
        
        ##### end of all calculations relating to females who can become pregnant on i_day
      }
      }
      
      
      # calculate energy intake of mother
      rho_sim[i_day,alives[i_sim]] <- (F_sim[i_day,alives[i_sim]])/(F_sim[i_day,alives[i_sim]] + Sa_sim[i_day,alives[i_sim]])
      if (rho_sim[i_day,alives[i_sim]] <= 0) rho_sim[i_day,alives[i_sim]] <- 0.000001
      Ir_sim[i_day,alives[i_sim]] <- assim_effic[age_days[i_day,alives[i_sim]]]*R_sim[day_index,i_sim]*Sa_sim[i_day,alives[i_sim]]^(2/3)/(1 + exp(-eta*(rho_t[i_day]/rho_sim[i_day,alives[i_sim]] - 1)))
      Cm_sim[i_day,alives[i_sim]] <- K*Sigma_M*(Sa_sim[i_day,alives[i_sim]]+Theta_F*F_sim[i_day,alives[i_sim]])^0.75
      starv_prop[i_day,alives[i_sim]] <- (1 - xi_m)*(rho_sim[i_day,alives[i_sim]] - rho_s[i_day])/((rho_t[i_day] - rho_s[i_day]) - xi_m*(rho_sim[i_day,alives[i_sim]] - rho_s[i_day]))
      if (rho_sim[i_day,alives[i_sim]] < rho_s[i_day]) starv_prop[i_day,alives[i_sim]] <- 0
      
      #####STARVATION-RELATED MORTALITY OF 1+ FEMALES REMOVED IN VERSION 1.3 FOR COMPATABILITY WITH ABC OUTPUTS
      # check if reserves fall below threshold for starvation-related mortality
#     if (rho_sim[i_day,alives[i_sim]] < rho_s[i_day]) {
#        starv_prop[i_day,alives[i_sim]] <- 0
#        is_female_dead <- f.starve_death(rho_sim[i_day,alives[i_sim]],rho_s[i_day])
#        if (is_female_dead == 0) {
#          population_size[i_years,11,boot_number] <- population_size[i_years,11,boot_number] + 1
#          age_days[(i_day+1):maximum_days,alives[i_sim]] <- 0
#          life_expectancy[alives[i_sim]] <- age_days[i_day,alives[i_sim]]
#        }
#     }
      
      #### end of first set of calculations for females that are still alive
      }
      
      
      ###### calculate appropriate values for current calf (if any)
      if (calf_age[i_day,alives[i_sim]] > 0){
        
        #### start of lactation calculations
        Sa_C_sim <- Sa[calf_age[i_day,alives[i_sim]]]
        # set calf's energy reserves to initial level if calf_age = 1
        if (calf_age[i_day,alives[i_sim]]==1) F_C[i_day,alives[i_sim]] <- F_start
#        rho_C[i_day,alives[i_sim]] <- F_C[i_day,alives[i_sim]]/(F_C[i_day,alives[i_sim]] + Sa_C[i_day,alives[i_sim]])
        rho_C[i_day,alives[i_sim]] <- F_C[i_day,alives[i_sim]]/(F_C[i_day,alives[i_sim]] + Sa_C_sim)
        if (rho_C[i_day,alives[i_sim]] <= 0) rho_C[i_day,alives[i_sim]] <- 0.000001
        # check for starvation-related death of calf
        if (rho_C[i_day,alives[i_sim]] < rho_s[i_day]) {
          is_calf_dead <- f.starve_death(rho_C[i_day,alives[i_sim]],rho_s[i_day])
          if (is_calf_dead == 0) 
          {calf_deaths[i_years,alives[i_sim]] <- calf_age[i_day,alives[i_sim]]
          calf_age[(i_day+1):maximum_days,alives[i_sim]] <- 0
          calf_life_expectancy[i_years,alives[i_sim]] <- calf_age[i_day,alives[i_sim]]
          }
        }
        # determine demand for milk
        calf_demand[i_day,alives[i_sim]] <- 1/(1 + exp(-eta*(rho_t[i_day]/rho_C[i_day,alives[i_sim]] - 1)))
        calf_demand[i_day,alives[i_sim]] <- ifelse(calf_demand[i_day,alives[i_sim]]>1,1,calf_demand[i_day,alives[i_sim]])
        # determine milk intake
#        Im_C[i_day,alives[i_sim]] <- phi_L*calf_demand[i_day,alives[i_sim]]*milk_prop[calf_age[i_day,alives[i_sim]]]*starv_prop[i_day,alives[i_sim]]*Sa_C[i_day,alives[i_sim]]^(2/3)
        Im_C[i_day,alives[i_sim]] <- phi_L*calf_demand[i_day,alives[i_sim]]*milk_prop[calf_age[i_day,alives[i_sim]]]*starv_prop[i_day,alives[i_sim]]*Sa_C_sim^(2/3)
        
        
        # determine calf's food intake
 #       Ir_C[i_day,i_sim] <- R_sim[day_index,i_sim]*assim_effic[calf_age[i_day,alives[i_sim]]]*Sa_C[i_day,alives[i_sim]]^(2/3)/(1 + exp(-eta*(rho_t[i_day]/rho_C[i_day,alives[i_sim]] - 1)))
        Ir_C[i_day,alives[i_sim]] <- R_sim[day_index,i_sim]*assim_effic[calf_age[i_day,alives[i_sim]]]*Sa_C_sim^(2/3)/(1 + exp(-eta*(rho_t[i_day]/rho_C[i_day,alives[i_sim]] - 1)))
        # determine cost of maintenance. 
#        Cm_C[i_day,alives[i_sim]] <- K*Sigma_M*(Sa_C[i_day,alives[i_sim]] + F_C[i_day,alives[i_sim]]*Theta_F)^0.75
        Cm_C[i_day,alives[i_sim]] <- K*Sigma_M*(Sa_C_sim + F_C[i_day,alives[i_sim]]*Theta_F)^0.75
        # determine calf's surplus
        surplus_C <- Im_C[i_day,alives[i_sim]] + Ir_C[i_day,alives[i_sim]] - Cm_C[i_day,alives[i_sim]] - CGa[calf_age[i_day,alives[i_sim]]]
        efficC <- ifelse(surplus_C > 0, epsi_plus[i_day], epsi_minus[i_day])
        # modify calf's energy reserves
        F_C[(i_day+1),alives[i_sim]] <- F_C[i_day,alives[i_sim]] + surplus_C/efficC
        
        #### end of lactation calculation
        }
      
      # determine surplus for females that are still alive and modify their reserves accordingly
      if (age_days[i_day,alives[i_sim]]>0){
        
      #### start of surplus calculations for females
      
      CLact_sim <- Im_C[i_day,alives[i_sim]]/Sigma_L
      CGa_sim <- CGa[age_days[i_day,alives[i_sim]]]
      surplus_sim  <- Ir_sim[i_day,alives[i_sim]] - Cm_sim[i_day,alives[i_sim]] -  CLact_sim - CGest_sim[i_day,alives[i_sim]] - CGa_sim
      effic <- ifelse(surplus_sim > 0, epsi_plus[i_day], epsi_minus[i_day])
      F_sim[(i_day+1),alives[i_sim]] <- F_sim[i_day,alives[i_sim]] + surplus_sim/effic
      if (i_day > maximum_days) F_sim[(i_day+1),alives[i_sim]] <- 0
      if (F_sim[(i_day+1),alives[i_sim]] < 0) F_sim[(i_day+1),alives[i_sim]] <- 0
      
            #### end surplus calculations for females
            }
          
          ### end of loop over females
          }

      ## end of loop over days of year
      }

  #### determine which calves will still be alive at the beginning of the next year 
  #calves_alive <- which(calf_life_expectancy[i_years,]>365)
  calves_alive <- which(calf_age[i_day,]>0)
  population_size[i_years,3,boot_number] <- length(which(age_days[implant_days[i_years],]>min_age))
  population_size[i_years,2,boot_number] <- population_size[i_years,4,boot_number] - population_size[i_years,3,boot_number]
  ## randomly assign a gender to each surviving calf
  calf_gender <- rbinom(length(calves_alive),1,propfemale)
  female_calves_alive <- calves_alive[which(calf_gender==1)]
  population_size[i_years,7,boot_number] <- length(calves_alive)/population_size[i_years,1,boot_number]
  conceptions <- length(which(conceive_record[i_years,]>0))
  population_size[i_years,5,boot_number] <- conceptions/population_size[i_years,3,boot_number]
  
  #### check that there are some female calves still alive
  if(length(female_calves_alive)>0){
  #### determine indices of columns that need to be added to accommodate these calves
  extra_columns <- ncol(F_sim) + c(1:length(female_calves_alive))
  ### matrix to augment "big" matrices with 1000's of rows
  augment <- matrix(c(0), nrow=(maximum_days+1), ncol = length(female_calves_alive))
  #### augment one_plus matrices
  Sa_sim <- cbind(Sa_sim, augment)
  F_sim <- cbind(F_sim, augment)
  rho_sim <- cbind(rho_sim, augment)
  Ir_sim <- cbind(Ir_sim, augment)
  Cm_sim <- cbind(Cm_sim, augment)
  CGest_sim <- cbind(CGest_sim, augment)
  starv_prop <- cbind(starv_prop, augment)
  age_days <-cbind(age_days, augment)
  F_preg <- cbind(F_preg, augment)
  # augment calf matrices
#  Sa_C <- cbind(Sa_C, augment)
  F_C <- cbind(F_C, augment)
  rho_C <- cbind(rho_C, augment)
  Ir_C <- cbind(Ir_C, augment)
  calf_demand <- cbind(calf_demand, augment)
  Im_C <- cbind(Im_C, augment)
  Cm_C <- cbind(Cm_C, augment)
  calf_age <- cbind(calf_age, augment)
  
  ##### Transfer relevant calf vectors to age one-plus matrices
  F_sim[(i_day+1),extra_columns] <- F_C[(i_day+1),female_calves_alive]
  F_C[(i_day+1),female_calves_alive] <- 0
  life_expectancy[extra_columns] <- calf_life_expectancy[(i_years), female_calves_alive]
  # transfer calf_ages of each female calf
    for (i_age_days in extra_columns[1]:extra_columns[length(extra_columns)]){
      ages_in_days <- c(366:life_expectancy[i_age_days])
      
      # check if simulated female is predicted to be still alive at the end of the simulation
      # if so, only simulate the days remaining in the simulation
      last_day <- i_day + length(ages_in_days)
      if (last_day>nrow(age_days)) last_day <- nrow(age_days)
      age_days[((i_day+1):last_day),i_age_days] <- ages_in_days[1:length((i_day+1):last_day)]
      
      # fill in the structural mass of individual i_life on those days
      Sa_sim[((i_day+1):last_day),i_age_days] <- Sa[ages_in_days[1:length((i_day+1):last_day)]]
      # end of i_age_days loop
      }

  ### matrix to augment "small" matrices with 10's of rows
  augment2 <- matrix(c(0), nrow=(years+1), ncol = length(female_calves_alive))
  # augment the matrices
  conceive_record <- cbind(conceive_record, augment2)
  calf_life_expectancy <- cbind(calf_life_expectancy, augment2)
  calf_deaths <- cbind(calf_deaths, augment2)
#  calf_days_record <- cbind(calf_days_record,augment2)
  #### end of matrix augmentation
  }
  
  ##### END OF JUMP OUT WHEN THERE ARE NO FEMALES ALIVE #####
  }
  # end of loop over years
  }    
population_size[(i_years+1),4,boot_number] <- length(which(age_days[(maximum_days+1),]>0))


