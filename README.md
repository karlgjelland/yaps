# YAPS - (Yet Another Positioning Solver) 
Based on the original [YAPS](https://www.nature.com/articles/s41598-017-14278-z.pdf) presented in Baktoft, Gjelland, Økland & Thygesen (2017): Positioning of aquatic animals based on time-of-arrival and random walk models using YAPS (Yet Another Positioning Solver). [DOI:10.1038/s41598-017-14278-z](https://www.nature.com/articles/s41598-017-14278-z.pdf)  

<h3>This branch introduces the option to do <b>on-the-fly synchronization</b> of the hydrophone array </h3>

YAPS now includes an (so far very experimental!) inline sync-function that enables on-the-fly synchronization of the hydrophone array. The hope is, that this will enable more users to give YAPS a try. For the time being, the sync option will use the hydro with most detections as “time-keeper” and adjust time drift on other hydros using a second order polynomial function. Input to YAPS is still a TOA-matrix with one column per ping (including missed pings to ensure burst interval consistency) and one row per hydro. 

So far the on-the-fly synchronization have only been tested on simulated data! However, it should work on real data as well, assuming they adhere to the TOA-requirements.  Also note, that for large projects and especially for long tracks, it is still recommended to do a separate synchronization of the array prior to running YAPS. Simulations indicate that the inline-sync works ok, but it will in all cases be better to do a separate array synchronization before running data through YAPS. A few things to consider when using inline-sync:
  * The model will be a bit more unstable and jiggly, so proper (=lucky?) initial values gets more important. Quite often this can be overcome, by rerunning getInp() using sdInits = 1 as this will assign new initial parameter values.
  * Overall performance will decrease as more data are used to estimate uninteresting parameters instead of fish tracks. This is especially true in situations where number of hydros detecting each ping is relatively low (i.e. mean number < 3). Even more so if using transmitters with random burst interval compared to stable burst interval.
  * If time drift on individual hydros is larger than 0.5*transmitter burst interval, detections in TOA-matrix might shift one or more phase. This will naturally affect performance negatively – the extent will depend on burst interval and fish behavior. 


To use on own data, compile a toa-matrix based on hydrophone data and replace the hydros dataframe with actual hydrophone positions. 

The yaps package requires [devtools](https://cran.r-project.org/web/packages/devtools/index.html) and [TMB](https://github.com/kaskr/adcomp).  
Please see for [instructions](https://github.com/kaskr/adcomp/wiki/Download) on TMB installation. Remember to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) as specified in the TMB documentation.

Installation:
```
install.packages('devtools')
devtools::install_github("baktoft/yaps", ref="inline_sync")
```

Usage example:

```
rm(list=ls())	
library(yaps)
set.seed(42)

# Simulate true track of animal movement of n seconds
trueTrack <- simTrueTrack(model='crw', n = 15000, deltaTime=1, shape=1, scale=0.5, addDielPattern=TRUE)

# Simulate telemetry observations from true track.
# Format and parameters depend on type of transmitter burst interval (BI) - stable (sbi) or random (rbi).
pingType <- 'rbi'

if(pingType == 'sbi') { # stable BI
	sbi_mean <- 15; sbi_sd <- 1e-4;
	teleTrack <- simTelemetryTrack(trueTrack, ss='rw', pingType=pingType, sbi_mean=sbi_mean, sbi_sd=sbi_sd)
} else if(pingType == 'rbi'){ # random BI
	pingType <- 'rbi'; rbi_min <- 20; rbi_max <- 40;
	teleTrack <- simTelemetryTrack(trueTrack, ss='rw', pingType=pingType, rbi_min=rbi_min, rbi_max=rbi_max)
}

# Simulate hydrophone array
hydros <- simHydros(auto=TRUE, trueTrack=trueTrack)

# If TRUE, drift of hydrophone internal clocks will be introduced in the simulations. Make sure that this is also corrected for in estimations by setting inline_sync = TRUE
clockDrift <- TRUE
inline_sync <- clockDrift


# Get TOA-matrix
toa_list <- simToa(teleTrack, hydros, pingType, sigmaToa=1e-4, pNA=0.25, pMP=0.01, clockDrift=clockDrift)
toa <- toa_list$toa

# Specify whether to use ss_data from measured water temperature (ss_data_what <- 'data') or to estimate ss in the model (ss_data_what <- 'est')
ss_data_what <- 'est'
if(ss_data_what == 'data') {ss_data <- teleTrack$ss} else {ss_data <- 0}

if(pingType == 'sbi'){
	inp <- getInp(hydros, toa, E_dist="Mixture", n_ss=20, pingType=pingType, sdInits=0, ss_data_what=ss_data_what, ss_data=ss_data, inline_sync=inline_sync)
} else if(pingType == 'rbi'){
	inp <- getInp(hydros, toa, E_dist="Mixture", n_ss=20, pingType=pingType, sdInits=0, rbi_min=rbi_min, rbi_max=rbi_max, ss_data_what=ss_data_what, ss_data=ss_data, inline_sync=inline_sync)
} 
str(inp)

pl <- c()
maxIter <- ifelse(pingType=="sbi", 500, 5000)
outTmb <- runTmb(inp, maxIter=maxIter, getPlsd=TRUE, getRep=TRUE)
str(outTmb)

# Estimates in pl
pl <- outTmb$pl

# Error estimates in plsd
plsd <- outTmb$plsd

# plot the resulting estimated track
plot(y~x, data=trueTrack, type="l", xlim=range(hydros$hx), ylim=range(hydros$hy), asp=1)
lines(y~x, data=teleTrack)
points(hy~hx, data=hydros, col="green", pch=20, cex=3)
lines(pl$Y~pl$X, col="red")

# check the on-the-fly synchronization - based on second order polynomial function: TOA_correct = TOA + offset + slope*toa + slope2*toa*toa
plot((toa_list$sync_offset-toa_list$sync_offset[inp$datTmb$timeKeeper+1] ~ pl$sync_offset)); abline(a=0,b=1)
plot((toa_list$sync_slope-toa_list$sync_slope[inp$datTmb$timeKeeper+1] ~ pl$sync_slope)); abline(a=0,b=1)
plot((toa_list$sync_slope2-toa_list$sync_slope2[inp$datTmb$timeKeeper+1] ~ pl$sync_slope2)); abline(a=0,b=1)
```
