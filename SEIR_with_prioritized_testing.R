
# Note: "rep" is currently set to 1000 simulations for each parameter set; 
# this is very slow (~24 hours). To run the code faster, set "rep" to a lower number.


rm(list = ls())
library("purrr")
library("dplyr") #SMA
library(pROC) # for calculating AUC 

filename = "SEIR_model.pdf" 

##### upload CPR validation data #####
{
  data <- read.csv("MainModel.csv") #to see results for Model A, B, or C, use ModelA.csv, etc. 
  data0 <- data[data$true == 0,]
  data1 <- data[data$true == 1,]
}

##### specify parameters #####
{
# simulation parameters
Tmax <- 2*365 # time limit of simulations 
step.size <- 1 # time steps (days)
rep <- 1000 # number of stochastic simulations to run for prioritized/random testing 
n.seed.events <- 15 # number of infected people at start of simulation
initial.state <- c(S=3200000-n.seed.events, E=0, I=n.seed.events, R=0, Ts=0, Te=0, Ti=0, Tr=0) 
# NOTE: cannot start off with any nonzero Ts, Te, Ti, or Tr, because they won't be in the "waiting.people" list

# model is frequency dependent, so we modify beta based on the total population size
beta.divisor <- as.numeric(initial.state["S"]+initial.state["E"]+initial.state["I"]+initial.state["R"])

##### testing parameters #####
w <- c(0.0013, 0.0013, 0.072, 0.00084) # proportion of S, E, I, R that want testing
Ntests <- 1000 # number of tests available
theta <- 2 # test turnaround time (days)
eta <- 0.2 # reduction in transmissibility for infected awaiting test results


# define colors
{
  RedCustom <- "#6D1F00"
  RedCustom.transp <- "#6D1F0040"
  BlueCustom <- "#003366"
  BlueCustom.transp <- "#00336626"
  GreyCustom <- "#858585"
  GreyCustom.transp <- "#85858540"
  YellowCustom <- "#B09500"
  YellowCustom.transp <- "#B0950040"
  LightRedCustom <- "#ff895a"
  LightRedCustom.transp <- "#ff895a40"
  DarkGreen = "#3A5431"
}

}

################################################################################

runSEIR <- function( prior.test=TRUE, final.only=FALSE) {
  
  # SEIR parameters
  {
  sigma <- 1/5.2      # sigma = 1/5.2 per day E -> I
  gamma.days <- runif(1, min=4, max=7)
  gamma <- 1/gamma.days        ## gamma = 1/6 per day I -> R
  beta <- gamma*R0
  
  #create SEIR parameter vector.
  param <- c(beta=beta/beta.divisor, sigma=sigma, gamma=gamma)
  
  t <- 0
  y <- c(initial.state)
  seir.output <- matrix(ncol=13, nrow=1)
  colnames(seir.output) <- c("time", "S", "E", "I","R","Ts", "Te", "Ti", "Tr","n0","n1","AUC","R0")
  seir.output[1,] <- c(t,y,NA,NA,NA,NA)
  
  waiting.people <- matrix(ncol=7,nrow=0)
  colnames(waiting.people) <- c(colnames(data0),"states")
  waiting.people[,"states"] <- factor(c("S","E","I","R"))
  }
  
  while (t < Tmax){
    t <- t+step.size 
    #print(t)
    
    ### STEP 1 - NEW PEOPLE GET TESTED #####################################################
    ##### COVID TESTING SIMULATION #####
    Test.eligible.S <- rbinom(1,y["S"],w[1])
    Test.eligible.E <- rbinom(1,y["E"],w[2])
    Test.eligible.I <- rbinom(1,y["I"],w[3])
    Test.eligible.R <- rbinom(1,y["R"],w[4])
    n0 <- Test.eligible.S + Test.eligible.E + Test.eligible.R # nubmer of '0's wanting testing
    n1 <- Test.eligible.I # number of '1's wanting testing
    
    # set up patient "pool" of people who all want testing
    patients0 <- data0[sample(nrow(data0),n0,replace=TRUE),]
    patients1 <- data1[sample(nrow(data1),n1,replace=TRUE),]
    
    states <- as.factor(rep(c("S","E","R","I"), c(Test.eligible.S,Test.eligible.E,Test.eligible.R,Test.eligible.I))) # vector of associated states
    levels(states)[c("S","E","I","R")] <- c("S","E","I","R")
    # need "0 states" listed first, followed by "1 states", since that's how they are combined below
    
    samp <- rbind(patients0, patients1)
    samp <- cbind(samp,states)
    if(n1 != 0 & n0 != 0){
      AUC <- as.numeric(suppressMessages(auc(true ~ pred_RF, data=samp)))
    }else{
      AUC <- NA
    }

    if (prior.test){ # if doing prioritized testing
      samp <- samp[order(-samp$pred_RF),]
      tested <- samp[1:min(nrow(samp),Ntests),]
    }else{ # if testing randomly
      tested <- sample_n(samp, min(nrow(samp),Ntests),replace=FALSE)
    }
    
    # number of people in each state that get tested:
    dS.test <- sum(tested$states == "S")
    dE.test <- sum(tested$states == "E")
    dI.test <- sum(tested$states == "I")
    dR.test <- sum(tested$states == "R")
    
    # move people into the "waiting for test results" bins
    {
    y["S"] <- y["S"]-dS.test
    y["E"] <- y["E"]-dE.test
    y["I"] <- y["I"]-dI.test
    y["R"] <- y["R"]-dR.test
    y["Ts"] <- y["Ts"]+dS.test
    y["Te"] <- y["Te"]+dE.test
    y["Ti"] <- y["Ti"]+dI.test
    y["Tr"] <- y["Tr"]+dR.test
    }
    
    # add record of tested people to list of people awaiting test results (need to track their states explicitly)
    waiting.people <- rbind(waiting.people,tested)
    
    ### STEP 2 - DISEASE DYNAMICS ##################################################33
    
    # dynamics for SEIR states
    p.expose <- 1-exp(-step.size*(param["beta"]*y["I"] + param["beta"]*eta*y["Ti"]))
    p.infect <- 1-exp(-step.size*param["sigma"])
    p.recover <- 1-exp(-step.size*param["gamma"])
    
    exposed.cases <- rbinom(1, y["S"],p.expose) # S->E
    incident.cases <- rbinom(1, y["E"],p.infect) # E->I
    recovered.cases <- rbinom(1, y["I"], p.recover) # I->R
    
    # dynamics for Ts,Te,Ti,Tr states (same as p.expose above, except with the additional eta)
    p.test.expose <- 1-exp(-step.size*(eta*param["beta"]*y["I"] + eta*param["beta"]*eta*y["Ti"])) 
    p.test.infect <- 1-exp(-step.size*param["sigma"]) 
    p.test.recover <- 1-exp(-step.size*param["gamma"]) 

    exposed.tested.cases <- rbinom(1, y["Ts"],p.test.expose) # Ts->Te
    incident.tested.cases <- rbinom(1, y["Te"],p.test.infect) # Te->Ti
    recovered.tested.cases <- rbinom(1, y["Ti"],p.test.recover) # Ti->Tr
    
    # disease dynamics move people between bins
    {
    y["S"] <- y["S"]-exposed.cases
    y["E"] <- y["E"]+exposed.cases-incident.cases
    y["I"] <- y["I"]+incident.cases-recovered.cases
    y["R"] <- y["R"]+recovered.cases
    y["Ts"] <- y["Ts"]-exposed.tested.cases
    y["Te"] <- y["Te"]+exposed.tested.cases-incident.tested.cases
    y["Ti"] <- y["Ti"]+incident.tested.cases-recovered.tested.cases
    y["Tr"] <- y["Tr"]+recovered.tested.cases
    }
    # change state labels on people in waiting.people list (important for when they get test results back)
    
    if(exposed.tested.cases > 0){
      if(length(which(waiting.people$states == "S"))==1 & exposed.tested.cases==1){
        waiting.people[which(waiting.people$states=="S"),]$states="E"}
      if(length(which(waiting.people$states == "S"))>1){
        waiting.people[sample(which(waiting.people$states == "S"),exposed.tested.cases), ]$states="E"}
    }
    
    if(incident.tested.cases > 0){
      if(length(which(waiting.people$states == "E"))==1 & incident.tested.cases==1){
        waiting.people[which(waiting.people$states=="E"),]$states="I"}
      if(length(which(waiting.people$states == "E"))>1){
        waiting.people[sample(which(waiting.people$states == "E"),incident.tested.cases), ]$states="I"}
    }
    
    if(recovered.tested.cases > 0){
      if(length(which(waiting.people$states == "I"))==1 & recovered.tested.cases==1){
        waiting.people[which(waiting.people$states=="I"),]$states="R"}
      if(length(which(waiting.people$states == "I"))>1){
        waiting.people[sample(which(waiting.people$states == "I"),recovered.tested.cases), ]$states="R"}
    }
      

    ### STEP 3 - PEOPLE GET TEST RESULTS BACK AND RETURN TO SEIR BINS ##################################################
    p.testresults <- 1-exp(-step.size*(1/theta)) # probability of getting test results back on a given day
    
    # number of people in each state Ts,Te,Ti,Tr getting results back
    testresults.S <- rbinom(1, y["Ts"],p.testresults) 
    testresults.E <- rbinom(1, y["Te"],p.testresults) 
    testresults.I <- rbinom(1, y["Ti"],p.testresults) 
    testresults.R <- rbinom(1, y["Tr"],p.testresults) 
    
    if(testresults.S > 0){
      if(length(which(waiting.people$states == "S"))==1 & testresults.S==1){
        waiting.people <- waiting.people[-which(waiting.people$states=="S"),]}
      if(length(which(waiting.people$states == "S"))>1){
        waiting.people <- waiting.people[-sample(which(waiting.people$states == "S"),testresults.S), ]}
    }
    TestS2S <- testresults.S
    
    if(testresults.E > 0){
      if(length(which(waiting.people$states == "E"))==1 & testresults.E==1){
        waiting.people <- waiting.people[-which(waiting.people$states=="E"),]}
      if(length(which(waiting.people$states == "E"))>1){
        waiting.people <- waiting.people[-sample(which(waiting.people$states == "E"),testresults.E), ]}
    }
    TestE2E <- testresults.E
    
    if(testresults.I > 0){
      if(length(which(waiting.people$states == "I"))==1 & testresults.I==1){
        I.results <- waiting.people[which(waiting.people$states=="I"),]
        waiting.people <- waiting.people[-which(waiting.people$states=="I"),]}
      if(length(which(waiting.people$states == "I"))>1){
        I.ind <- sample(which(waiting.people$states == "I"),testresults.I)
        I.results <- waiting.people[I.ind,]
        waiting.people <- waiting.people[-I.ind, ]}
      TestI2I <- nrow(I.results[I.results$true == 0,]) # people now in I who got tested when they weren't in I
      TestI2R <- nrow(I.results[I.results$true == 1,]) # people now in I who got tested when they were in I
    } else {
      TestI2I<-0
      TestI2R<-0
      }
    
    if(testresults.R > 0){
      if(length(which(waiting.people$states == "R"))==1 & testresults.R==1){
        R.results <- waiting.people[which(waiting.people$states=="R"),]
        waiting.people <- waiting.people[-which(waiting.people$states=="R"),]}
      if(length(which(waiting.people$states == "R"))>1){
        R.ind <- sample(which(waiting.people$states == "R"),testresults.R)
        R.results <- waiting.people[R.ind,]
        waiting.people <- waiting.people[-R.ind, ]}
    }
    TestR2R <- testresults.R
    
    # test results move people between bins
    {
    y["S"] <- y["S"]+TestS2S
    y["E"] <- y["E"]+TestE2E
    y["I"] <- y["I"]+TestI2I
    y["R"] <- y["R"]+TestI2R+TestR2R
    y["Ts"] <- y["Ts"]-TestS2S
    y["Te"] <- y["Te"]-TestE2E
    y["Ti"] <- y["Ti"]-TestI2I-TestI2R
    y["Tr"] <- y["Tr"]-TestR2R
    }
    
    #trick to speed up the code
    if (!final.only) {
      seir.output<-rbind(seir.output, c(t,y,n0,n1,AUC,R0))
    }
  }
  
  if(final.only) {
    seir.output[1,]<-c(t,y,n0,n1,auc,R0)
  }
  return(seir.output)
}

###################### RUN MODEL WITH PRIORITIZED TESTING ########
pdf(file=filename,width=7.5, height=11)
par(oma=c(2,2,1,2)) # bottom, left, top, right
m <- matrix(c(1,1,2,3,4,5,6,7,8,9),nrow=5,ncol=2,byrow=TRUE)
layout(mat = m,heights = c(0.125,0.3,0.3,0.3,0.3))

### add legend at top of plot
{
par(mar = c(1,1,1,1))# bottom, left, top, right
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x="topright", inset=c(0,0.025),
       #text.width=c(0.154,0.15,0.1),
       legend=c("S", "E", "I", "R"),
       col=c(BlueCustom, LightRedCustom, RedCustom, GreyCustom),
       bty="n",lwd=2, cex=1.5, horiz = TRUE)

legend(x="bottomright", inset=c(0,0.1),
       #text.width=c(0.154,0.15,0.1),
       legend=c("prioritized testing", "indiscriminate testing"),
       col=c("black","black"),lty= c("solid","21"),
       bty="n",lwd=2, cex=1.5, horiz=TRUE)
}

counter = 1
hosp.dat <- list() # saves hospital data
names = c("A","B","C","D*","E","F","G*","H")

# run through each of the R and wI parameters listed here, as in Fig. 3
R0.list <- c(2.5, 2.25, 2, 1.75, 1.5)
wI.list <- c(0.029, 0.072, 0.144)

for (counter in seq(1,8)){
    
  print(counter)
  if (counter <= 5){
    R0 <- R0.list[counter]
    param.txt <- bquote(paste('R'['e']*' = ', .(R0),', w'['I']*' = ', .(w[3])))
  } else {
    R0 <- 1.75
    w[3] <- wI.list[counter-5]
    #if (counter==8){
    #  Tmax <- 2*365 # final plot needs to run it longer
    #}
    param.txt <- bquote(paste('R'['e']*' = ', .(R0),', w'['I']*' = ', .(w[3])))
  }
  
    
# set up storage vectors/matrices
{
output.table.p.t <- list()
S.vals.p.t <- numeric()
E.vals.p.t <- numeric()
I.vals.p.t <- numeric()
R.vals.p.t <- numeric()
n0.vals.p.t <- numeric()
n1.vals.p.t <- numeric()
AUC.mean.p.t <- numeric()
R0.vals.p.t <- numeric()
}

# run stochastic simulations
for(i in 1:rep){
  #print(i)
  output.table.p.t[[i]] <- as.data.frame(runSEIR(prior.test=TRUE, final.only=FALSE))
  S.vals.p.t <- cbind(S.vals.p.t,output.table.p.t[[i]][,2]+output.table.p.t[[i]][,6])
  E.vals.p.t <- cbind(E.vals.p.t,output.table.p.t[[i]][,3]+output.table.p.t[[i]][,7])
  I.vals.p.t <- cbind(I.vals.p.t,output.table.p.t[[i]][,4]+output.table.p.t[[i]][,8])
  R.vals.p.t <- cbind(R.vals.p.t,output.table.p.t[[i]][,5]+output.table.p.t[[i]][,9])
  n0.vals.p.t <- cbind(n0.vals.p.t,output.table.p.t[[i]][,10])
  n1.vals.p.t <- cbind(n1.vals.p.t,output.table.p.t[[i]][,11])
  AUC.mean.p.t <- cbind(AUC.mean.p.t,mean(na.omit(output.table.p.t[[i]][,12])))
  R0.vals.p.t <- cbind(R0.vals.p.t,output.table.p.t[[i]][,13])
}

# save mean values
{
  S.mean.p.t <- rowMeans(S.vals.p.t)
  E.mean.p.t <- rowMeans(E.vals.p.t)
  I.mean.p.t <- rowMeans(I.vals.p.t)
  R.mean.p.t <- rowMeans(R.vals.p.t)
  n0.means.p.t <- rowMeans(n0.vals.p.t)
  n1.means.p.t <- rowMeans(n1.vals.p.t)
  AUC.means.p.t <- mean(AUC.mean.p.t)
  #print(AUC.means.p.t)
  R0.means.p.t <- rowMeans(R0.vals.p.t)
}
results<-output.table.p.t 


#### calculate lower/upper 95th percentiles ####
{
  Svals = NULL
  Evals = NULL
  Ivals = NULL
  Rvals = NULL
  
  for (i in 1:rep) {
    Svals = rbind(results[[i]][, 2]+results[[i]][, 6], Svals)
    Evals = rbind(results[[i]][, 3]+results[[i]][, 7], Evals)
    Ivals = rbind(results[[i]][, 4]+results[[i]][, 8], Ivals)
    Rvals = rbind(results[[i]][, 5]+results[[i]][, 9], Rvals)
  }
  
  Svals <- apply(Svals, 2, sort)
  Evals <- apply(Evals, 2, sort)
  Ivals <- apply(Ivals, 2, sort)
  Rvals <- apply(Rvals, 2, sort)
  
  S_lower <- Svals[max(ceiling(rep * 0.025), 1), ] / 1000000
  S_upper <- Svals[floor(rep * 0.975), ] / 1000000
  E_lower <- Evals[max(ceiling(rep * 0.025), 1), ] / 1000000
  E_upper <- Evals[floor(rep * 0.975), ] / 1000000
  I_lower <- Ivals[max(ceiling(rep * 0.025), 1), ] / 1000000
  I_upper <- Ivals[floor(rep * 0.975), ] / 1000000
  R_lower <- Rvals[max(ceiling(rep * 0.025), 1), ] / 1000000
  R_upper <- Rvals[floor(rep * 0.975), ] / 1000000
}


### PLOT ############################################

par(mar = c(3,3,3,3)) # margin: bottom, left, top, right

plot(x=results[[i]]$time, y=(results[[i]]$S/1000000), type="n",
     xlab="", ylab="",
     ylim=c(0,sum(initial.state)/1000000),cex.axis=1.25,cex.lab=1.25)

polygon(c(results[[i]]$time,
          rev(results[[i]]$time)),
        c(S_upper,
          rev(S_lower)),
        col=BlueCustom.transp, border = NA)
polygon(c(results[[i]]$time,
          rev(results[[i]]$time)),
        c(E_upper,
          rev(E_lower)),
        col=LightRedCustom.transp, border = NA)
polygon(c(results[[i]]$time,
          rev(results[[i]]$time)),
        c(I_upper,
          rev(I_lower)),
        col=RedCustom.transp, border = NA)
polygon(c(results[[i]]$time,
          rev(results[[i]]$time)),
        c(R_upper,
          rev(R_lower)),
        col=GreyCustom.transp, border = NA)

lines(S.mean.p.t/1000000, col=BlueCustom, lwd=2)
lines(E.mean.p.t/1000000, col=LightRedCustom, lwd=2)
lines(I.mean.p.t/1000000, col=RedCustom, lwd=2)
lines(R.mean.p.t/1000000, col=GreyCustom, lwd=2)


######### NOW WITH RANDOM TESTING ##########################
{
  output.table.r.t <- list()
  S.vals <- numeric()
  E.vals <- numeric()
  I.vals <- numeric()
  R.vals <- numeric()
  n0.vals <- numeric()
  n1.vals <- numeric()
  AUC.means.r.t <- numeric()
  R0.vals <- numeric()
}

for(i in 1:rep){
  output.table.r.t[[i]] <- as.data.frame(runSEIR(prior.test=FALSE, final.only=FALSE))
  S.vals <- cbind(S.vals,output.table.r.t[[i]][,2]+output.table.r.t[[i]][,6])
  E.vals <- cbind(E.vals,output.table.r.t[[i]][,3]+output.table.r.t[[i]][,7])
  I.vals <- cbind(I.vals,output.table.r.t[[i]][,4]+output.table.r.t[[i]][,8])
  R.vals <- cbind(R.vals,output.table.r.t[[i]][,5]+output.table.r.t[[i]][,9])
  n0.vals <- cbind(n0.vals,output.table.r.t[[i]][,10])
  n1.vals <- cbind(n1.vals,output.table.r.t[[i]][,11])
  AUC.means.r.t <- cbind(AUC.means.r.t,mean(na.omit(output.table.r.t[[i]][,12])))
  R0.vals <- cbind(R0.vals,output.table.r.t[[i]][,13])
}

{
  S.mean <- rowMeans(S.vals)
  E.mean <- rowMeans(E.vals)
  I.mean <- rowMeans(I.vals)
  R.mean <- rowMeans(R.vals)
  n0.means <- rowMeans(n0.vals)
  n1.means <- rowMeans(n1.vals)
  AUC.means.r.t <- mean(AUC.means.r.t)
  print(AUC.means.r.t)
  R0.means <- rowMeans(R0.vals)
}

results<-output.table.r.t #SMA

##############################3
# calculate lower and upper 95th percentiles for each curve
{
  Svals = NULL
  Evals = NULL
  Ivals = NULL
  Rvals = NULL
}

for(i in 1:rep){
  Svals = rbind(results[[i]][,2]+results[[i]][,6], Svals)
  Evals = rbind(results[[i]][,3]+results[[i]][,7], Evals)
  Ivals = rbind(results[[i]][,4]+results[[i]][,8], Ivals)
  Rvals = rbind(results[[i]][,5]+results[[i]][,9], Rvals)
}

#### calculate confidence bounds ####
{
Svals <- apply(Svals, 2, sort)
Evals <- apply(Evals,2,sort)
Ivals <- apply(Ivals,2,sort)
Rvals <- apply(Rvals,2,sort)

S_lower <- Svals[max(round(rep*0.025),1),]/1000000
S_upper <- Svals[round(rep*0.975),]/1000000
E_lower <- Evals[max(round(rep*0.025),1),]/1000000
E_upper <- Evals[round(rep*0.975),]/1000000
I_lower <- Ivals[max(round(rep*0.025),1),]/1000000
I_upper <- Ivals[round(rep*0.975),]/1000000
R_lower <- Rvals[max(round(rep*0.025),1),]/1000000
R_upper <- Rvals[round(rep*0.975),]/1000000
}

# add randomized testing lines to the plot
{
polygon(c(results[[i]]$time,
          rev(results[[i]]$time)),
        c(S_upper,
          rev(S_lower)),
        col=BlueCustom.transp, border = NA)
polygon(c(results[[i]]$time,
          rev(results[[i]]$time)),
        c(E_upper,
          rev(E_lower)),
        col=LightRedCustom.transp, border = NA)
polygon(c(results[[i]]$time,
          rev(results[[i]]$time)),
        c(I_upper,
          rev(I_lower)),
        col=RedCustom.transp, border = NA)
polygon(c(results[[i]]$time,
          rev(results[[i]]$time)),
        c(R_upper,
          rev(R_lower)),
        col=GreyCustom.transp, border = NA)

lines(S.mean.p.t/1000000, col=BlueCustom, lwd=2)
lines(E.mean.p.t/1000000, col=LightRedCustom, lwd=2)
lines(I.mean.p.t/1000000, col=RedCustom, lwd=2)
lines(R.mean.p.t/1000000, col=GreyCustom, lwd=2)
lines(S.mean/1000000, col=BlueCustom, lwd=2,lty="21")
lines(E.mean/1000000, col=LightRedCustom, lwd=2,lty="21")
lines(I.mean/1000000, col=RedCustom, lwd=2,lty="21")
lines(R.mean/1000000, col=GreyCustom, lwd=2,lty="21")
}

#### Output statistics from each parameter set #####
{
# calculate difference in I "peak" timing
max.p.t <- which(I.mean.p.t == max(I.mean.p.t))
max.random <- which(I.mean == max(I.mean))

# calculate difference in I "peak" height
height.diff <- max(I.mean)-max(I.mean.p.t)

# and difference in final R values
R.diff <- R.mean[Tmax+1]-R.mean.p.t[Tmax+1]

{
print(paste("Difference in peak timing is ",max.p.t-max.random," days"))
print(paste("Peak of random testing curve is ",height.diff," people higher than for prioritized testing."))
print(paste("That is ",100*round(height.diff/max(I.mean.p.t),2),"% higher"))
print(paste("By end of simulation, ",R.diff," fewer people were infected under the prioritized testing scenario than random."))
print(paste("That means a",round(100*(1-R.mean.p.t[Tmax+1]/R.mean[Tmax+1])),"% reduction in the total cumulative number of infected people."))
}

### Calculate limitations in terms of hospital capacity ###
Hosp.rate <- 10*0.005 # hospitalization rate is 10*IFR (so 10*0.005=0.05 --> 5% )
N.hospitalizations.p.t <- Hosp.rate*I.mean.p.t
D1 <- sum(N.hospitalizations.p.t > 4869) # days exceeding hostpital capacity
D3 <- pmax(N.hospitalizations.p.t-4869,0)

N.hospitalizations <- Hosp.rate*I.mean
D2 <- sum(N.hospitalizations > 4869) 
D4 <- pmax(N.hospitalizations-4869,0)
{
print(paste("The number of patients seeking hospital beds will exceed capacity for",D2-D1,"fewer days, if using prioritized testing"))
print(paste("That is a",100*(D2-D1)/D2,"% reduction in days"))
print(paste("This saves",sum(D4)-sum(D3),"people-days above hospital capacity"))
print(paste("That is a",100*(sum(D4)-sum(D3))/sum(D4),"% reduction in people-days above hospital capacity"))
}
### Now ICU capacity
# assume 14% of hospitalized patients need ICU
N.ICU.p.t <- N.hospitalizations.p.t*0.14
N.ICU <- N.hospitalizations*0.32
D5 <- sum(N.ICU.p.t > 687)
D6 <- sum(N.ICU > 687)
D7 <- pmax(N.ICU.p.t-687,0)
D8 <- pmax(N.ICU-687,0)
{
print(paste("The number of patients needing ICU beds will exceed capacity for",D6-D5,"fewer days, if using prioritized testing"))
print(paste("That is a",100*(D6-D5)/D6,"% reduction in days"))
print(paste("This saves",sum(D8)-sum(D7),"people-days above ICU capacity"))
print(paste("That is a",100*(sum(D8)-sum(D7))/sum(D8),"% reduction in people-days above hospital capacity"))
}
}
mtext(names[counter], side=3, adj=0, line=.5, cex=1.25); 
mtext(param.txt, side=3, adj=1, line=.5, cex=1.25)

# save values of N.hospitalizations, N.ICU for each parameter set
hosp.dat[[counter]] <- rbind(N.hospitalizations.p.t, N.hospitalizations, N.ICU.p.t, N.ICU)

counter = counter+1
}
{
mtext("time (days)",side=1,line=0,outer=TRUE,cex=1.3)
mtext("number of people (millions)",side=2,line=0,outer=TRUE,cex=1.3,las=0)
}
dev.off()




###########################################################3
#### plot hospital/ICU capacity ####
pdf(file="hospital_ICU_demand_with_delayed_results.pdf",width=7.5, height=11*.8)
par(oma=c(2,2,1,2)) # bottom, left, top, right
m <- matrix(c(1,1,2,3,4,5,6,7),nrow=4,ncol=2,byrow=TRUE)
layout(mat = m,heights = c(0.1255,0.3,0.3,0.3))

### add legend at top of plot
{
  par(mar = c(1,1,1,1))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x="topright", inset=c(0,0.025),
         legend=c("daily hospital occupancy", "daily ICU occupancy"), 
         col=c("grey60","DarkGreen"),
         bty="n",lwd=2, cex=1.5, horiz = TRUE)
  
  legend(x="bottomright", inset=c(0.04,0.1),
         legend=c("prioritized testing", "indiscriminate testing"),
         col=c("black","black"),lty= c("solid","21"),
         bty="n",lwd=2, cex=1.5, horiz=TRUE)
}

names = c("A","B","C","D","E")
ICU.capacity <- 687
hosp.capacity <- 4869

counter <- 1
while (counter <= 5){
  hosp.df <- hosp.dat[[counter]]
  N.hospitalizations.p.t <- hosp.df[1,]
  N.hospitalizations <- hosp.df[2,]
  N.ICU.p.t <- hosp.df[3,]
  N.ICU <- hosp.df[4,]

  par(mar = c(3,3,3,3)) # margin: bottom, left, top, right
  times <- seq(1,length(N.hospitalizations))
  plot(times, N.hospitalizations[times], type="n",
     xlab="", ylab="",
     ylim=c(0,21000),cex.axis=1.25,cex.lab=1.25)
  lines(times,N.hospitalizations.p.t[times], col="grey60", lwd=2)
  lines(times,N.hospitalizations[times],col="grey60",lty="21",lwd=2)
  abline(h=hosp.capacity, col="grey60",lwd=2)

  lines(times,N.ICU.p.t[times],col="DarkGreen",lwd=2)
  lines(times,N.ICU[times],col="DarkGreen",lty="21",lwd=2)
  abline(h=ICU.capacity, col="DarkGreen",lwd=2)
  mtext(names[counter], side=3, adj=0, line=.5, cex=1.25); 
  
  
  param.txt <- bquote(paste('R'['e']*' = ', .(R0.list[counter])))
  mtext(param.txt, side=3, adj=1, line=.5, cex=1.25)
  counter <- counter+1
}

mtext("time(days)",side=1,line=0,outer=TRUE,cex=1.3)
mtext("number of people",side=2,line=0,outer=TRUE,cex=1.3,las=0)
dev.off()

