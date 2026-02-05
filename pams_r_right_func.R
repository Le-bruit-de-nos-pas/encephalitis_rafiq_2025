library(smacof)
source("BootSmacof.R")

# Read the data from the “csv” format
testdata <- read.csv(file="Cross-sectional.csv", header=TRUE) 

# Label the input variables
colnames(testdata) <- c("OV1","NS2","VA3","LP4","PP5","SR6","VS7","GI8","CF9","NR10","NP11","NW12","VA13","PR14","AS15","ON16","PC17","MW18")
testname <- c("OV1","NS2","VA3","LP4","PP5","SR6","VS7","GI8","CF9","NR10","NP11","NW12","VA13","PR14","AS15","ON16","PC17","MW18")

nprofile <- 3 # number of profile dim

# Choose “ordinal” among “ratio”, “interval”, “ordinal”, “mspline” for scaling
preout <- smacofSym(dist(t(testdata)), nprofile, type="ordinal")
preout$conf   # profile

# Stress value
preout$stress # for the 3-dim solution 

# Profile plot 
par(mfrow=c(1, nprofile), mar=c(2.9,4.1,2,1))
for (i in 1:nprofile) {
    plot(scale(preout$conf[,i]), main=paste("Dim", i, "Scaled profile"), xlab="", ylab="Coordinate", type="b")
    abline(h=0, lty=3)
  }

# To conduct “profileR,” using the “classical” option in PAMS in R.
# However, this analysis can be skipped.
library(profileR)
result0 <- pams(testdata, nprofile)
result0$dimensional.configuration 

for (i in 1:nprofile) plot(scale(preout$conf[,i]), scale(result0$dimensional.configuration[,i]), xlab="Smacof", ylab="ProfileR")

# The results from profileR package and the “classical” option are the same
preout2 <- cmdscale(dist(t(testdata)), k = nprofile)

for (i in 1:nprofile) {
  plot(scale(preout2[,i]), -scale(result0$dimensional.configuration[,i]), 
       xlab="Classical", ylab="ProfileR", main=paste("Corr =", cor(scale(preout2[,i]), -scale(result0$dimensional.configuration[,i]))))
  abline(0,1)
}




direction <- c(-1, 1, 1) # flip the direction of the 1st profile
for (i in 1:nprofile) preout$conf[,i] <- preout$conf[,i] * direction[i]

# profile plot with flipped first profile 
par(mfrow=c(1, nprofile), mar=c(2.9,4.1,2,1))
for (i in 1:nprofile) {
  plot(scale(preout$conf[,i]), main=paste("Dim", i, "Scaled profile"), xlab="", ylab="Coordinate", type="b")
  abline(h=0, lty=3)
}


# Generate an empirical distribution of 2,000 core profiles 
# from the original I x J data set by bootstrapping
set.seed(1)     # to duplicate the result, set the seed as 1.
nprofile <- 3   # the number of core profiles

# With a 3-dimensional solution, we set the seed as 1 
# to conduct bootstrapping with nBoot = 2000, but a user can choose 
# any number (e.g., nBoot = 1000).
args(BootSmacof)

# For example, we include the first ten people for further investigation
# Researchers can specify any participants of interest
participant <- c(1:10) # identify specific participants
empprofile <- BootSmacof(testdata=testdata, participant=participant, mds="smacof", type="ordinal", distance="euclid", 
                           scale=FALSE, nprofile=nprofile, direction=direction, cl=0.95, nBoot=2000, testname=testname, file="profile")


# Save the result. Next time, just load the saved result
#save(empprofile, file = "empprofile.RData")
#load(file = "empprofile.RData")


#Outputs
names(empprofile)

## PART 2: Estimate Correlations between Subscale Mean & Core Profiles 
## and BCa CIs for Core Profile Coordinates

# Estimate how much each core profile is accounted for 
# by the other core profiles to examine collinearity among core profiles
# For example, Dim1 = a*Dim2 + b*Dim3 + error to compute R^2 for Dim1
empprofile$MDSR2 # R^2 for each core profile

# Pairwise latent profile plot  
pairs(empprofile$MDS$conf)

# Compare patterns for the 18 cognitive subscales mean profile and 
# the 1st/2nd/3rd core profiles 

cor(apply(testdata, 2, mean), empprofile$MDS$conf[, 1]) # correlation
cor(apply(testdata, 2, mean), empprofile$MDS$conf[, 2]) # correlation
cor(apply(testdata, 2, mean), empprofile$MDS$conf[, 3]) # correlation

# The bootstrap summary statistics for each core profile
# For example, we chose i = 3 (for the 3rd core profile)
# “BCaLower BCaUpper” for the 95% BCa CI

i <- 1  # for the i^th dimension profile
round(empprofile$MDSsummary[[i]], 3)

# Additional plots for 95% BCa CI of dimension profiles
i <- 3  # for the i^th dimension profile
par(mfrow=c(1,1), mar=c(4.1,4.1,2,1))
plot(empprofile$MDSsummary[[i]]$Ori, type="b", xlab="", ylab="Coordinate", main="Dim_3: 95% BCa CI",
            ylim=range(empprofile$MDSsummary[[i]]$BCaLower, empprofile$MDSsummary[[i]]$BCaUpper))
lines(empprofile$MDSsummary[[i]]$BCaLower, lty=4)
lines(empprofile$MDSsummary[[i]]$BCaUpper, lty=4)
abline(h=0, lty=5, col="red")

# Test a specific jth coordinate in the ith profile for the normality
# since normality cannot be guaranteed for all coordinates
# using SE for test statistics (t-test) would not be recommended
i <- 3 # specify the ith profile.
j <- 12 # specify the jth coordinate of the ith profile.
par(mfrow=c(1,2), mar=c(2.9,4.1,2,1))
hist(empprofile$MDSprofile[[i]][, j], main=paste("Histogram of", testname[j], "coordinate of Dim", i))
qqnorm(empprofile$MDSprofile[[i]][, j], main=paste("Q-Q plot of", testname[j], "coordinate of Dim", i))

## PART 3: Estimating R^2, Person Weights, & Person Correlations

# The mean R^2 for all 1650 participants
round(empprofile$WeightmeanR2, 2) 

# To illustrate, the first 10 participants’ person weights, levels, R^2s, 
# and person correlations were included. Person #9 was included in the paper
# as an example
round(empprofile$Weight[1:10,], 2)

## Part 4: Assessment of A Specific Person (Person #9 for example)

# Plot of standardized person profile and dim profile
# Person #9 profile juxtaposed with the 1st core profile (CP1)
# because this person’s correlation with CP1 was r = .93 
# but the other correlations were r = .13 and r = -.07.
p <- 9 # for the p^th participant
par(mfrow=c(1, nprofile), mar=c(2.9,4.1,2,1))
for (i in 1:nprofile) {
  plot(scale(as.numeric(testdata[p,])), main=paste("Dim", i, "profile"), xlab="", ylab="Coordinate", type="b")
  lines(scale(empprofile$MDS$conf[,i]), lty=2)
  abline(h=0, lty=3)
}


