
###### THIS IS THE CODE FOR EVALUATING EFFECTS OF DISTRUBANCE 
##### IT WILL ONLY WORK AFTER PopDyn_DEB_undisturbed_beta.R HAS BEEN RUN!!!!!

# days_disturbed will contain number of days on which a female is disturbed in each year when there is disturbance
days_disturbed <- matrix(c(0), ncol=ncol(conceive_record_d),nrow=(years+1))
# affected_by_PTS is a flag indicating whether or not a female has been affected 
# affected_by_PTS = value of i_day on which PTS occurred
# fert_success_d documents affect of PTS (if any) on fertility
affected_by_PTS <- rep(c(0),ncol(conceive_record_d))
died_from_starvation <- rep(c(0),ncol(conceive_record_d))
fert_success_d <- rep(c(1),ncol(conceive_record_d))
## assign all females to a vulnerable local population
subpop_number <- rep(c(1),ncol(conceive_record_d))
## rescale vulnmean to account for animals not affected by any operation
new_vulnmean <- vulnmean/rescale
if (nvulnmean>1){
  all_females <- c(1:length(subpop_number))
  numbers <- rmultinom(1,length(subpop_number),new_vulnmean)
  for (i_subs in 1:(nvulnmean-1)) {
    subpop <- sample(all_females,numbers[i_subs])
    subpop_number[subpop] <- i_subs
    all_females <- all_females[-which(all_females%in%subpop)]
    }
    subpop_number[all_females] <- nvulnmean
}

## determine begin and end points of each piling year in pile
begin_dist <- seq(1,nrow(pile),365)
end_dist <- seq(365,nrow(pile),365)

dist_year <- 1
# LOOP OVER YEARS
for (i_years in burnin_end:years){

  ### determine population size on the first day of simulated year
  alives <- which(age_days_d[begin_year[i_years],]>0)
  matures <- length(which(age_days_d[begin_year[i_years],]>min_age))
  population_size_d[i_years,4,boot_number] <- length(alives)
  calves_born_d <- length(which(calf_age_d[begin_year[i_years],]==1))
  population_size_d[i_years,1,boot_number] <- calves_born_d
  population_size_d[i_years,6,boot_number] <- calves_born_d/matures
  calves_affected_by_PTS <- rep(c(0),ncol(conceive_record_d))
  ########## CHECK THAT THERE ARE STILL SOME FEMALES ALIVE ######
  if (length(alives)>0) {
    Rmean_dd <- Rmean*(dd_slope*length(alives)/sim_number + dd_intercept)
        #### set up and fill the R_sim matrix
  R_sim <- matrix(Rmean_dd, nrow = 365, ncol = length(alives))
   # if stochastic_R = TRUE fill annual R_sim matrix with random values
  if (stochastic_R) R_sim <- sapply(1:length(alives), function(j) Rmean_dd*rbeta(365,a,b)/mu)

  if (dist_year < pile_years)  {
    # check if there is actually any piling in dist_year
    if (sum(pile$pvec[begin_dist[dist_year]:end_dist[dist_year]])>0) {
          # new.PTS will contain first day of PTS for females that acquire PTS in current year
          new.PTS <- new.PTS.calves <- rep(c(0),ncol(conceive_record_d))
          # PTS.by_day documents any PTS events within each year
          PTS.by_day <- disturbed.by_day <- array(c(0), dim = c(365,ncol(conceive_record_d),nvulnmean))
          # loop over vulnerable populations
              for (i_pop in 1:nvulnmean){
              # extract the days in piling year when there is piling that may affect members of the vulnerable population
              disturbed_days <- which(pile[begin_dist[dist_year]:end_dist[dist_year],(5+i_pop)]>0)
              ########## check if the vulnerable population is exposed to piling
                if (length(disturbed_days>0)){
                  disturbees <- which(subpop_number[alives]==i_pop)
                    ########### loop over days when there is piling
                    for (X in 1:length(disturbed_days)) {
                    # determine who is disturbed on each disturbed day
                    effects <- rbinom(length(disturbees),1,(1- p_disturb_by_year[disturbed_days[X],dist_year,i_pop]))
                    disturbed.by_day[disturbed_days[X],alives[disturbees[which(effects==0)]],i_pop] <- 1
                    # determine if PTS might occur
                    PTS.effects <- rbinom(length(which(effects==0)),1,p_pts_by_year[disturbed_days[X],dist_year,i_pop])
                    PTS.by_day[disturbed_days[X],alives[disturbees[which(effects==0)]],i_pop] <- PTS.effects
          
                    if (fixed_disturbance_duration) {effects[which(effects==0)] <- (24-disturbance_duration)/24} else {
                      ##### code to allow the response duration to vary as predicted by an Erlang distribution
                      effects[which(effects==0)] <- (24-rgamma(length(which(effects==0)),shape= Erlang.shape, rate = Erlang.shape/disturbance_duration))/24
                      effects[which(effects<0)] <- 0
                    ##### end of variable response duration calculation
                    }
 
                           # modify the appropriate values of R_sim
                  R_sim[disturbed_days[X],disturbees] <- R_sim[disturbed_days[X],disturbees] * effects
          
                  ############ end of X loop over days when there is piling 
                  }
              ########## end of loop for vulnerable populations that are exposed to piling in i_years
              } 
        ##### JUMP to here if no piling in a particular year
        
      # determine total number of days of disturbance experienced by each female in local population i_pop
      sum.disturbed <- sapply(c(1:ncol(conceive_record_d)), function(D) sum(disturbed.by_day[,D,i_pop]))
      days_disturbed[i_years,] <- days_disturbed[i_years,] + sum.disturbed
    
      # determine total number of days on which each female in local population i_pop was exposed to sound loud enough to cause PTS
      sum.PTS <- sapply(c(1:ncol(conceive_record_d)), function(all_females) sum(PTS.by_day[,all_females,i_pop]))
      # determine first day on which PTS is experienced
      min.PTS_day <- sapply(which(sum.PTS>0), function(Z) min(which(PTS.by_day[,Z,i_pop]>0)))
      sum.PTS[which(sum.PTS>0)] <- min.PTS_day
          if(length(which(sum.PTS>0)) > 0){
          # assign PTS to calves of all females that were exposed to sound loud enough to cause PTS
          sum.PTS.calves <- sum.PTS
          # remove females that don't have calves
          sum.PTS.calves[which(calf_age_d[begin_year[i_years],]!=1)] <- 0
          # remove calves that had already been weaned before first PTS exposure
          sum.PTS.calves[which(sum.PTS.calves>Tl)] <- 0
          new.PTS.calves[which(sum.PTS.calves>0)] <- sum.PTS.calves[which(sum.PTS.calves>0)]
          # remove females that already have PTS
          sum.PTS[which(affected_by_PTS>0)] <- 0
          # update new.PTS
          new.PTS[which(sum.PTS>0)] <- sum.PTS[which(sum.PTS>0)]
          }
        ##### end of i_pop loop over different vulnerable local populations
        }
    # identify females experiencing PTS level exposure for the first time
    PTS.females <- which(new.PTS>0)
    # update life expectancy and fertility of these females
    if (length(PTS.females>0)){
      for (i_life.PTS in 1:length(PTS.females)){
        # current last simulated day for female 
        end.life <- max(which(age_days_d[,PTS.females[i_life.PTS]]>0))
        # revised life expectancy
        new.life <- min(which(PTS.cum_Surv<=cum_Surv[life_expectancy_d[PTS.females[i_life.PTS]]]))
        # update life expectancy
        life_expectancy_d[PTS.females[i_life.PTS]] <- new.life
        # zero days that will no longer be lived
          if(new.life < age_days_d[end.life,PTS.females[i_life.PTS]]) {
          days_to_remove <- age_days_d[end.life,PTS.females[i_life.PTS]] - new.life
          age_days_d[(end.life - days_to_remove + 1):end.life,PTS.females[i_life.PTS]] <- 0
          }
        }
      }
    fert_success_d[PTS.females] <- PTS.fert_success
    
    PTS.calves <- which(new.PTS.calves>0)
    # update life expectancy of calves and calf_days
    if (length(PTS.calves>0)){
      for (c_life in 1:length(PTS.calves)){
        # check if calf was still alive when it is predicted to get PTS 
          if (new.PTS.calves[PTS.calves[c_life]] > calf_life_expectancy_d[i_years,PTS.calves[c_life]]){
          # new life expectancy
          new.life <- min(which(PTS.cum_Surv<=cum_Surv[calf_life_expectancy_d[i_years,PTS.calves[c_life]]]))
          calf_life_expectancy_d[i_years,PTS.calves[c_life]] <- new.life}
        # check if calf_age_d needs to be reduced because calf is predicted to die in first year
        if(new.life < 365) calf_age_d[begin_year[i_years]:(begin_year[i_years]+new.life-1),PTS.calves[c_life]] <- c(1:new.life)
        #### end of loop over calves with PTS
        }
      }
    affected_by_PTS[PTS.females] <- new.PTS[PTS.females]
    calves_affected_by_PTS[PTS.calves] <- new.PTS.calves[PTS.calves]
    if (length(days_disturbed[i_years,which(days_disturbed[i_years,]>0)])){
    population_size_d[i_years,8,boot_number] <- mean(days_disturbed[i_years,which(days_disturbed[i_years,]>0)])}
    population_size_d[i_years,9,boot_number] <- length(PTS.females)
        
    ###### end of calculations in years when there is some piling
    }

        dist_year <- dist_year + 1
  ###### end of calculations during years when there might be piling
  }

  R_sim.mean[i_years] <- population_size_d[(i_years),12,boot_number] <- mean(R_sim)
  
  day_index <- 0
  #### LOOP OVER DAYS IN YEAR i_years
  for (i_day in begin_year[i_years]:end_year[i_years]){
    day_index <- day_index+1
    
    ###### LOOP OVER FEMALES THAT ARE ALIVE AT START OF i_years
    for (i_sim in 1:length(alives)){
      
      #### CHECK THAT FEMALE IS STILL ALIVER
      if (age_days_d[i_day,alives[i_sim]] > 0){
        
        ##### START OF LOOP FOR FEMALES THAT ARE ABLE TO BECOME PREGNANT ON i_day ############ 
        
        ###### check if female is old enough to become pregnant and implantation is possible
        if (age_days_d[i_day,alives[i_sim]] > min_age & preg_possible[i_day] > 0){
          ####### check if female has enough reserves to become pregnant
          F_preg_d[i_day,alives[i_sim]] <- F_neo + rho_s[i_day]*Sa[age_days_d[i_day,alives[i_sim]]]/(1-rho_s[i_day])
          if (F_sim_d[i_day,alives[i_sim]] >= F_preg_d[i_day,alives[i_sim]]){
        
            ######## check if ovum implants successfully
            if (rbinom(1,1,fert_success_d[alives[i_sim]])==1){
          
              ###### start of calculations for females that implant successfully
              tmp  <- f.foetusandcalf_life(i_sim) 
              conceive <- i_day
              foetal_life <- tmp$life_expect_f
              # store foetal life for simulated female
              conceive_record_d[i_years,alives[i_sim]] <- foetal_life
              calf_life_expectancy_d[(i_years+1), alives[i_sim]] <- tmp$life_expect_c
          
              # "birth_day" is day on which female gives birth OR day on which foetus dies
              birth_day <- conceive + foetal_life
              ####### check that female survives until birth_day
              if (age_days_d[birth_day,alives[i_sim]]==0) calf_life_expectancy_d[(i_years+1), alives[i_sim]] <- 0
              ####### end of check on female survival
            
              ######## check if there is enough time for female to give birth before the simulation ends
              if (birth_day > maximum_days){
                CGest_sim_d[conceive:maximum_days, alives[i_sim]] <- CGest[1:length(conceive:maximum_days)]} else 
                  {
                  CGest_sim_d[conceive:(birth_day-1),alives[i_sim]] <- CGest[1:foetal_life]
          
                  ######## check if foetus or female dies before the end of the pregnancy (in which case calf life expectancy is 0)
                  if (calf_life_expectancy_d[(i_years+1), alives[i_sim]]>0)
                    {
                    ######### determine how many days of the calf's life will be simulated  
                    if (calf_life_expectancy_d[(i_years+1), alives[i_sim]] < 365) {calf_days <- calf_life_expectancy_d[(i_years+1), alives[i_sim]]} else {calf_days <- 365}
                    if ((birth_day+calf_days) > maximum_days+1) calf_days <- maximum_days - birth_day + 1
                    # insert calf ages, calf core mass and calf growth costs  
                    calf_age_d[birth_day:(birth_day+calf_days-1),alives[i_sim]] <- c(1:calf_days)
                    ######### end of calculations for foetuses that are actually born
                    ############################################################### 
                    }
                  }
          
          ###### end of calculations for females that implant successfully
          }
        
        ####### end of calculations relating to females who have enough reserves to become pregnant on i_day
        }
      ###### end of calculation for all females who could become pregnant
      }
      
      
      # calculate energy intake of mother
      rho_sim_d[i_day,alives[i_sim]] <- (F_sim_d[i_day,alives[i_sim]])/(F_sim_d[i_day,alives[i_sim]] + Sa_sim_d[i_day,alives[i_sim]])
      if (rho_sim_d[i_day,alives[i_sim]] <= 0) rho_sim_d[i_day,alives[i_sim]] <- 0.000001
      Ir_sim_d[i_day,alives[i_sim]] <- assim_effic[age_days_d[i_day,alives[i_sim]]]*R_sim[day_index,i_sim]*Sa_sim_d[i_day,alives[i_sim]]^(2/3)/(1 + exp(-eta*(rho_t[i_day]/rho_sim_d[i_day,alives[i_sim]] - 1)))
      Cm_sim_d[i_day,alives[i_sim]] <- K*Sigma_M*(Sa_sim_d[i_day,alives[i_sim]]+Theta_F*F_sim_d[i_day,alives[i_sim]])^0.75
      starv_prop_d[i_day,alives[i_sim]] <- (1 - xi_m)*(rho_sim_d[i_day,alives[i_sim]] - rho_s[i_day])/((rho_t[i_day] - rho_s[i_day]) - xi_m*(rho_sim_d[i_day,alives[i_sim]] - rho_s[i_day]))
      if (rho_sim_d[i_day,alives[i_sim]] < rho_s[i_day]) starv_prop_d[i_day,alives[i_sim]] <- 0
      
      ####### STARVATION-RELATED MORTALITY OF 1+ FEMALES ONLY OCCURS IN YEARS WHEN THERE IS PILING FOR COMPATABILITY WITH ABC OUTPUTS
              # check if reserves fall below threshold for starvation-related mortality
          if (rho_sim_d[i_day,alives[i_sim]] < rho_s[i_day] & (dist_year-1) < pile_years)  {
          is_female_dead <- f.starve_death(rho_sim_d[i_day,alives[i_sim]],rho_s[i_day])
          if (is_female_dead == 0) {
            died_from_starvation[alives[i_sim]] <- i_day
            population_size_d[i_years,11,boot_number] <- population_size_d[i_years,11,boot_number] + 1
            age_days_d[(i_day+1):maximum_days,alives[i_sim]] <- 0
            life_expectancy_d[alives[i_sim]] <- age_days_d[i_day,alives[i_sim]]
            }
          }
      ###### end of first set of calculations for females that are still alive
        }
      
        ###### calculate appropriate values for current calf (if any) and costs of lactation
        if (calf_age_d[i_day,alives[i_sim]] > 0){
          Sa_C_sim <- Sa[calf_age_d[i_day,alives[i_sim]]]
        
          # set calf's energy reserves to initial level if calf_age_d = 1
          if (calf_age_d[i_day,alives[i_sim]]==1) F_C_d[i_day,alives[i_sim]] <- F_start
          rho_C_d[i_day,alives[i_sim]] <- F_C_d[i_day,alives[i_sim]]/(F_C_d[i_day,alives[i_sim]] + Sa_C_sim)
          if (rho_C_d[i_day,alives[i_sim]] <= 0) rho_C_d[i_day,alives[i_sim]] <- 0.000001
        
          ####### check for starvation-related death of calf
          if (rho_C_d[i_day,alives[i_sim]] < rho_s[i_day]) {
            is_calf_dead <- f.starve_death(rho_C_d[i_day,alives[i_sim]],rho_s[i_day])
            if (is_calf_dead == 0) 
              {calf_deaths_d[i_years,alives[i_sim]] <- calf_age_d[i_day,alives[i_sim]]
              calf_age_d[(i_day+1):maximum_days,alives[i_sim]] <- 0
              calf_life_expectancy_d[i_years,alives[i_sim]] <- calf_age_d[i_day,alives[i_sim]]
              }
              ####### end of starvation loop
            }
          
        # determine demand for milk
        calf_demand_d[i_day,alives[i_sim]] <- 1/(1 + exp(-eta*(rho_t[i_day]/rho_C_d[i_day,alives[i_sim]] - 1)))
        calf_demand_d[i_day,alives[i_sim]] <- ifelse(calf_demand_d[i_day,alives[i_sim]]>1,1,calf_demand_d[i_day,alives[i_sim]])
        # determine milk intake
        Im_C_d[i_day,alives[i_sim]] <- phi_L*calf_demand_d[i_day,alives[i_sim]]*milk_prop[calf_age_d[i_day,alives[i_sim]]]*starv_prop_d[i_day,alives[i_sim]]*Sa_C_sim^(2/3)
        
        # determine calf's food intake
        Ir_C_d[i_day,alives[i_sim]] <- R_sim[day_index,i_sim]*assim_effic[calf_age_d[i_day,alives[i_sim]]]*Sa_C_sim^(2/3)/(1 + exp(-eta*(rho_t[i_day]/rho_C_d[i_day,alives[i_sim]] - 1)))
        # determine cost of maintenance. 
        Cm_C_d[i_day,alives[i_sim]] <- K*Sigma_M*(Sa_C_sim + F_C_d[i_day,alives[i_sim]]*Theta_F)^0.75
        # determine calf's surplus
        surplus_C <- Im_C_d[i_day,alives[i_sim]] + Ir_C_d[i_day,alives[i_sim]] - Cm_C_d[i_day,alives[i_sim]] - CGa[calf_age_d[i_day,alives[i_sim]]]
        efficC <- ifelse(surplus_C > 0, epsi_plus[i_day], epsi_minus[i_day])
        # modify calf's energy reserves
        F_C_d[(i_day+1),alives[i_sim]] <- F_C_d[i_day,alives[i_sim]] + surplus_C/efficC
        
        ###### end of calf-related calculations
        }
      
      ###### start of surplus calculations for females that are still alive and modify their reserves accordingly
      if (age_days_d[i_day,alives[i_sim]]>0){
        CLact_sim <- Im_C_d[i_day,alives[i_sim]]/Sigma_L
        CGa_sim <- CGa[age_days_d[i_day,alives[i_sim]]]
        surplus_sim  <- Ir_sim_d[i_day,alives[i_sim]] - Cm_sim_d[i_day,alives[i_sim]] -  CLact_sim - CGest_sim_d[i_day,alives[i_sim]] - CGa_sim
        effic <- ifelse(surplus_sim > 0, epsi_plus[i_day], epsi_minus[i_day])
        F_sim_d[(i_day+1),alives[i_sim]] <- F_sim_d[i_day,alives[i_sim]] + surplus_sim/effic
        if (i_day > maximum_days) F_sim_d[(i_day+1),alives[i_sim]] <- 0
        if (F_sim_d[(i_day+1),alives[i_sim]] < 0) F_sim_d[(i_day+1),alives[i_sim]] <- 0
         ###### end surplus calculations for females
        }
  
          ###### end of JUMP OUT if female dies during year
          }
        ##### end of loop over i_sim = all females in alives
        }
      #### end of loop over i_days of year
      }

  #### determine which calves will still be alive at the beginning of the next year 
  population_size_d[i_years,3,boot_number] <- length(which(age_days_d[implant_days[i_years],]>min_age))
  population_size_d[i_years,2,boot_number] <- population_size_d[i_years,4,boot_number] - population_size_d[i_years,3,boot_number]
  calves_alive <- which(calf_age_d[i_day,]>0)
  population_size_d[i_years,7,boot_number] <- length(calves_alive)/population_size_d[i_years,1,boot_number]
  ## randomly assign a gender to each surviving calf
  calf_gender <- rbinom(length(calves_alive),1,propfemale)
  female_calves_alive <- calves_alive[which(calf_gender==1)]
  conceptions <- length(which(conceive_record_d[i_years,]>0))
  population_size_d[i_years,5,boot_number] <- conceptions/population_size_d[i_years,3,boot_number]
  population_size_d[i_years,10,boot_number] <- length(which(affected_by_PTS[alives]>0))/population_size_d[i_years,4,boot_number] 
  
  #### if there are some female calves still alive augment the relevant matrices
  if(length(female_calves_alive)>0){
  #### determine indices of columns that need to be added to accommodate these calves
  extra_columns <- ncol(F_sim_d) + c(1:length(female_calves_alive))
  ### matrix to augment "big" matrices with 1000's of rows
  augment <- matrix(c(0), nrow=(maximum_days+1), ncol = length(female_calves_alive))
  #### augment one_plus matrices
  Sa_sim_d <- cbind(Sa_sim_d, augment)
  F_sim_d <- cbind(F_sim_d, augment)
  rho_sim_d <- cbind(rho_sim_d, augment)
  Ir_sim_d <- cbind(Ir_sim_d, augment)
  Cm_sim_d <- cbind(Cm_sim_d, augment)
  CGest_sim_d <- cbind(CGest_sim_d, augment)
  starv_prop_d <- cbind(starv_prop_d, augment)
  age_days_d <-cbind(age_days_d, augment)
  F_preg_d <- cbind(F_preg_d, augment)
  # augment calf matrices
#  Sa_C_d <- cbind(Sa_C_d, augment)
  F_C_d <- cbind(F_C_d, augment)
  rho_C_d <- cbind(rho_C_d, augment)
  Ir_C_d <- cbind(Ir_C_d, augment)
  calf_demand_d <- cbind(calf_demand_d, augment)
  Im_C_d <- cbind(Im_C_d, augment)
  Cm_C_d <- cbind(Cm_C_d, augment)
  calf_age_d <- cbind(calf_age_d, augment)
  
  ##### Transfer relevant calf vectors to age one-plus matrices and vectors
  F_sim_d[(i_day+1),extra_columns] <- F_C_d[(i_day+1),female_calves_alive]
  F_C_d[(i_day+1),female_calves_alive] <- 0
  life_expectancy_d[extra_columns] <- calf_life_expectancy_d[i_years, female_calves_alive]
  # assign surviving female calves to their mother's vulnerable local population
  subpop_number[extra_columns] <- subpop_number[female_calves_alive]
  # transfer PTS record of the mothers of the surviving female calves so they can't get PTS again
  affected_by_PTS[extra_columns] <- calves_affected_by_PTS[female_calves_alive]
  calves_affected_by_PTS[extra_columns] <- 0
  died_from_starvation[extra_columns] <- 0
  fert_success_d[extra_columns] <- 1
  if (length(which(affected_by_PTS[extra_columns]>0))) {
    fert_success_d[extra_columns[which(affected_by_PTS[extra_columns]>0)]] <- PTS.fert_success}
  
  # fill in ages_in_days for transferred calves
    for (i_age_days in extra_columns[1]:extra_columns[length(extra_columns)]){
      ages_in_days <- c(366:life_expectancy_d[i_age_days])
      
      # check if simulated female is predicted to be still alive at the end of the simulation
      # if so, only simulate the days remaining in the simulation
      last_day <- i_day + length(ages_in_days)
      if (last_day>nrow(age_days_d)) last_day <- nrow(age_days_d)
      age_days_d[((i_day+1):last_day),i_age_days] <- ages_in_days[1:length((i_day+1):last_day)]
      
      # fill in the structural mass of individual i_life on those days
      Sa_sim_d[((i_day+1):last_day),i_age_days] <- Sa[ages_in_days[1:length((i_day+1):last_day)]]
          #### end of i_age_days loop
          }

  ### matrix to augment "small" matrices with 10's of rows
  augment2 <- matrix(c(0), nrow=(years+1), ncol = length(female_calves_alive))
  # augment the matrices
  conceive_record_d <- cbind(conceive_record_d, augment2)
  days_disturbed <- cbind(days_disturbed, augment2)
  calf_life_expectancy_d <- cbind(calf_life_expectancy_d, augment2)
  calf_deaths_d <- cbind(calf_deaths_d, augment2)
        ### end of matrix augmentation
        }
  
  
   # end of loop over i_years
  }
population_size_d[(i_years+1),4,boot_number] <- length(which(age_days_d[(maximum_days+1),]>0))




