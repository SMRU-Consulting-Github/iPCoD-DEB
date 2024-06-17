for (boot_number in 1:nboot){
  
  # select ABC parameter values
  if (stochastic_R){
    ABC.sample <- sample(c(1:nrow(resFiltSR)),1)
    # shape parameter for calf foraging efficiency
    upsilon <- resFiltSR$upsilon[ABC.sample]
    # age (days) at which resource foraging becomes 50% efficient
    Tr <- resFiltSR$Tr[ABC.sample]
    # parameter defining strength of starvation related mortality
    mu_s <- resFiltSR$mu_s[ABC.sample]
    # day of lactation when female contribution begins to decrease
    Tn <- resFiltSR$Tn[ABC.sample]
    # mean resource density
    Rmean <- resFiltSR$Rmean[ABC.sample]
    # Field metabolic maintenance scalar
    Sigma_M <- resFiltSR$Sigma_M[ABC.sample]
  } else {
    ABC.sample <- sample(c(1:nrow(resFiltNoSR)),1)
    # shape parameter for calf foraging efficiency
    upsilon <- resFiltNoSR$upsilon[ABC.sample]
    # age (days) at which resource foraging becomes 50% efficient
    Tr <- resFiltNoSR$Tr[ABC.sample]
    # parameter defining strength of starvation related mortality
    mu_s <- resFiltNoSR$mu_s[ABC.sample]
    # day of lactation when female contribution begins to decrease
    Tn <- resFiltNoSR$Tn[ABC.sample]
    # mean resource density
    Rmean <- resFiltNoSR$Rmean[ABC.sample]
    # field metabolic maintenance scalar
    Sigma_M <- resFiltNoSR$Sigma_M[ABC.sample]
  }
  
  dd_slope <- 0
  # sample density dependent effects
  if (density_dependence) dd_slope <- sample(dd_slope_values,1)
  dd_intercept <- 1 - dd_slope
  
  source('Porpoise_iPCoD+DEB_Params.R')
  
  ######## RUN ACTUAL SIMULATION #################
  source('PopDyn_iPCoD+DEB_undisturbed.R')
  source('PopDyn_iPCOD+DEB_disturbed.R')
  ### store simulation results
  # dat.out[,1,] is the total size of the disturbed population immediately before there are any births (i.e. minimum population size)
  # dat.out[,2,] is the size of the undisturbed population, dat.out[,3,] is the difference between the two populations
  # dat.out[,3,] is the difference between the two
  t_index <- c(burnin_end:(years+1))
  dat.out[t_index,2, boot_number] <- round(population_size[t_index,4,boot_number]*pmean/sim_number)
  # if a proportion of the disturbed population is not affected by disturbance (i.e. 1 - (sum(vulnmean) < 1))
  # add a suitable portion of the undisturbed population to populations_size_d to make up the numbers
  dat.out[t_index,1, boot_number] <- round(population_size_d[t_index,4,boot_number]*rescale*pmean/sim_number) + round((1 - rescale)*dat.out[t_index,2, boot_number])
  dat.out[t_index,3,boot_number] <- dat.out[t_index,2,boot_number] - dat.out[t_index,1,boot_number]
  dat.out[t_index,4,boot_number] <- 100*((dat.out[t_index,1,boot_number]/(pmean))^(1/(t_index-1))-1)
  dat.out[t_index,5,boot_number] <- 100*((dat.out[t_index,2,boot_number]/(pmean))^(1/(t_index-1))-1)  
  
  #### end of loop over bootstraps
}
