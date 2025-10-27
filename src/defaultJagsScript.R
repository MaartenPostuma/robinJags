#Regular JAGS model
library(rjags)
library(R2jags)
library(ggplot2)

args<-commandArgs(trailingOnly=TRUE)

#load files:
#setwd("/home/rlexmond/Documents/Data")

countData <- read.csv(args[1],sep=",")#"countsSimple.csv", sep=",")

hourlyData <- read.csv(args[2],sep=",")#"hourlyDataWithCovariates2025July.csv", sep=",")


########################################################################
# Make a list with variables (containing numbers only) used in the model
#24 Added daylight as a variable

jagsData <- list(counts = countData$Order,
                 nSamples = nrow(countData),
                 nHours = nrow(hourlyData),
                 site = hourlyData$site,
                 tAnom = hourlyData$normNewTAnomaly,
                 julian = hourlyData$normJulian,
                 element = hourlyData$element,
                 Variable = hourlyData$Variable,
                 startRow = tapply(1:nrow(hourlyData),hourlyData$sample,min),
                 endRow = tapply(1:nrow(hourlyData),hourlyData$sample,max),
                 minutes = hourlyData$minutes)

# Construct the model 
# In dit model wordt het aantal minuten van onvolledige uren gewogen meegenomen.
sink(args[5])
cat("
model{
  # Per sample: how well do observed counts and predicted (accumulated over many hours) counts match
  
  for( i in 1:nSamples){
    counts[i] ~ dnegbin(p[i], r)
    mu[i] <- (sum(predictedHourlyCounts[startRow[i]:endRow[i]]))
    p[i] <- r/(r + mu[i])   
  }
  
  # Count prediction for every hour-sample combination based on the to-be-evaluated model parameter values
  
  for(i in 1:nHours){
    predictedHourlyCounts[i] <- exp(linearModel[i])*minutes[i]/60
    linearModel[i] <- alphaSite[site[i]] + betaT*tAnom[i] + betaJ1*julian[i] + betaJ2*julian[i]*julian[i] + betaVar*Variable[i] + betaE[element[i]]
    vr[i]<- exp( 2* predictedHourlyCounts[i]+ lvar ) * (exp(lvar)-1 )
  }
  
  # Priors
  alpha ~ dnorm(0,.01)
  betaT ~ dnorm(0,.01)
  betaJ1 ~ dnorm(0,.01)
  betaJ2 ~ dnorm(0,.01)
  betaVar ~ dnorm(0,.01)
  for (k in 1:24) {
    alphaSite[k] ~ dnorm(alpha, tauSite)
  }
  tauSite ~ dgamma(0.001, 0.001)
  
    betaE[1] <- 0
  for (x in 2:3) {
    betaE[x] ~ dnorm(0, .01) 
  }

  r ~ dunif(0.01, 100)
  sdhat ~ dunif(0,5)
  lvar <- pow(sdhat,2)
} # end of model
",fill=TRUE)
sink(NULL)

params <- c("alpha","tauSite","betaT","betaJ1","betaJ2","betaE", "betaVar", "alphaSite")

OrderVariable <- jags.parallel(data = jagsData,
                    inits = NULL,
                    parameters.to.save = params,
                    model.file = args[5],
                    n.chains = 3,
                    n.iter = 24000,
                    n.burnin = 2000,
                    n.thin = 10,
                    n.cluster=3)

saveRDS(OrderVariable,args[3])

png(args[4],width=1920,height=1024)

traceplot(OrderVariable,mfrow=c(3,12))






dev.off()


