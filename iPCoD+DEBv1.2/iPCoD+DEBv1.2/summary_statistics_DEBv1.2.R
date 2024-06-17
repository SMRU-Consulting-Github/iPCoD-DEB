# ---------------------------------------------------------
# BEFORE RUNNING THIS R SCRIPT DO THE FOLLOWING STEP
# 1) run "pcodControl.R"

# ---------------------------------------------------------
years <- length(dat.out[,1,1])

Undist_mean <- Dist_mean <- Diff_mean <- backtonormal <- max_diff <- 0
# x <- y <- 0

load("HP_ver5test_Output.rdata")
for (i in 1:(years-1)){Dist_mean[i] <- mean(dat.out[i,1,])}
for (i in 1:(years-1)){Undist_mean[i] <- mean(dat.out[i,2,])}
for (i in 1:(years-1)){Diff_mean[i] <- mean(dat.out[i,3,])}

# determine maximum difference between disturbed and undisturbed population for each simulation
# for (i in 1:nboot){max_diff[i] <- -min(dat.out[1:years,3,i])}
(max_diff <- apply(dat.out[1:years,3,],1,min))

(backtonormal <- which((apply(dat.out[1:years,3,],1,min))>-.01*mean(dat.out[1,2,])))

# this R file has the functions for creating figures and table of outputs from iPCoD. 
source("SUMMARY_FUNCTIONS_RJ.R")

# ---------------------------------------------------------
# Are there negative populations in the dat.out object/PCoD projections? If so, how many?
sum(dat.out[,2,]<0,na.rm=T) # number <0 in baseline populations
sum(dat.out[,1,]<0,na.rm=T) # number <0 in impacted populations
 

# ---------------------------------------------------------

### LINES ###
# create line plots of un-impacted and impacted population trajectories over the 25 year simulation

dev.off()
	par(mfrow=c(1,3),oma=c(3,5.5,1.5,1),mar=c(.5,.5,.5,.5))
	LINES(dat.out,simulated.years=2016:(2016+years-2),Dist_mean,Undist_mean, nlines=100) 

# save the LINES plots as a png file in the working folder
	dev.copy(png, 'LINES.png')
	dev.off()
	
	
	
# ---------------------------------------------------------
	
### POLYGONS ###
# create a graph of the mean un-impacted and impacted population trajectories and their 95% confidence intervals.
	
dev.off()
	par(mfrow=c(1,3),oma=c(3,5.5,1.5,1),mar=c(.5,.5,.5,.5))
	POLYGONS(dat.out,simulated.years=2016:(2016+years-2),Dist_mean=Dist_mean,Undist_mean=Undist_mean)

# save the POLYGON plots as a png file in the working folder
	dev.copy(png, 'POLYGONS.png')
	dev.off()
	
	
	
# ------------------------------------------------------------------
(possible.time.points=1:length(dat.out[!is.na(dat.out[,1,i]),1,i]))
# ------------------------------------------------------------------

	
	
# ---------------------------------------------------------
	
### CUSTOM.TABLE ###
# create a custom csv file with un-impacted and impacted population sizes in the user specified years.
# The following code gives a table for the mean un-impacted and impacted population size 
 # in year 1, 6, 12, 18 and 24, 95% CI and the median population size
# note that time.point=1 is the start year (no possible impact yet)
# time.point=2 means the start of 2nd year, after 1 year of impact
# the user can change the values in "time.points", "percent.CI" and "get.median"

customtable <- CUSTOM.TABLE(dat.out, time.points=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), percent.CI=95, get.median=F)

# output table as "csv" file
  write.csv(round(customtable),"custom_table.csv")



# ---------------------------------------------------------

### HISTOGRAMS ###
# create histograms of population size for the unimpacted and impacted populations
# note that time.point=1 is the start year (no possible impact yet)
# time.point=2 means the start of 2nd year, after 1 year of impact

dev.off()
par(mfrow=c(1,3),oma=c(3,5.5,1.5,3),mar=c(.5,.5,.5,.5))
HISTOGRAMS(dat.out, time.point=26, axis_breaks=seq(0,3000,100))

# save the HISTOGRAM plots as a png file in the working folder
  dev.copy(png, 'HISTOGRAMS.png')
  dev.off()

	
# ---------------------------------------------------------
	
#### POPULATION SIZE RATIOS ###
# create ratios of impacted:un-impacted population sizes.
# note that time.point=1 is the start year (no possible impact yet)
# time.point=2 means the start of 2nd year, after 1 year of impact

# Calculate the ratio of impacted:un-impacted population size in years 2 (after 1 yr impact), 7 (after 6 yrs impact), 13, 19 and 25
  ratio1<-dat.out[2,1,]/dat.out[2,2,]
  ratio2<-dat.out[7,1,]/dat.out[7,2,]
  ratio3<-dat.out[13,1,]/dat.out[13,2,]
  ratio4<-dat.out[19,1,]/dat.out[19,2,]
  ratio5<-dat.out[25,1,]/dat.out[25,2,]

  ratio1table<-summary(ratio1)
  ratio2table<-summary(ratio2)
  ratio3table<-summary(ratio3)
  ratio4table<-summary(ratio4)
  ratio5table<-summary(ratio5)

# bind the results across all specified years into 1 table object called ratiotable
  ratiotable<-rbind(ratio1table, ratio2table, ratio3table, ratio4table, ratio5table)
# view the table of results across all years specified
  ratiotable
# output the table of results across all specified years as a csv
  write.csv((ratiotable),"ratio_pop_size_table.csv")
	
# Create histogram of population size ratios
  dev.off()
	par(mfrow=c(2,3),oma=c(4.3,6.4,3.5,.5),mar=c(2.5,2.5,.5,.5))
	RATIOS(dat.out, time.point=2) 
	RATIOS(dat.out, time.point=7)
	RATIOS(dat.out, time.point=13)
	RATIOS(dat.out, time.point=19)
	RATIOS(dat.out, time.point=25)
	
# save the POPULATION SIZE RATIOS plots as a png file in the working folder
	dev.copy(png, 'RATIOS_POP_SIZE.png')
	dev.off()
	

	
	
# ---------------------------------------------------------
	
#### GROWTH RATE RATIOS ###
# the ratio of impacted:un-impacted annual growth rate
# note that time.point=1 is the start year (no possible impact yet)
# time.point=2 means the start of 2nd year, after 1 year of impact
	
# ratios for each time.point between 2 and 26 are calculated in the SUMMARY_FUNCTIONS file
	
# summary table
	GRratiotable
	
# output the table of results across all specified years as a csv
	write.csv((GRratiotable),"ratio_growth_rate_table.csv")
	
# plot histograms of the growth rate ratios for time.points 2, 6, 11, 15, 21 & 26
	dev.off()
	par(mfrow=c(2,3))
	hist(GRtime.point2, breaks=10)
	hist(GRtime.point6, breaks=10)
	hist(GRtime.point11, breaks=10)
	hist(GRtime.point15, breaks=10)
	hist(GRtime.point21, breaks=10)
	hist(GRtime.point26, breaks=10)

	# save the GROWTH RATE RATIOS plots as a png file in the working folder
	dev.copy(png, 'RATIOS_GROWTH_RATE.png')
	dev.off() 
	
	
	
# ---------------------------------------------------------

### CENTILES
# removed from current version of code due to errors
	
# ---------------------------------------------------------