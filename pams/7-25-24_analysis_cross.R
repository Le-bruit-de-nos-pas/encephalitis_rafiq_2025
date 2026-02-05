library(smacof)
setwd("/Users/sekangkim/Documents/J Sub/PAMSproject/") # the folder where data and function are stored. 

#### Load the R function 
source("BootSmacof.R") # load the R function

#### Load the R function 
#source("newBootSmacof.R") # load the R function

#### Reading data
testdata <- read.csv(file="Cross-sectional.csv", header=FALSE) 
colnames(testdata) <- c("OV1","NS2","VA3","LP4","PP5","SR6","VS7","GI8","CF9","NR10","NP11","NW12","VA13","PR14","AS15","ON16","PC17","MW18")
testname <- c("OV1","NS2","VA3","LP4","PP5","SR6","VS7","GI8","CF9","NR10","NP11","NW12","VA13","PR14","AS15","ON16","PC17","MW18")

# First of all, I am not familiar with the modelling method of MDS area at the conceptual level.
# MDS treats person's scores as variables and analyzes the relationship between those,
# while equation (1) regards person's scores as observations of a variable although person's scores measure distinct characteristics. 
# Is this procedure the conventional approach? I assume such convention, and think that reviewer's comments are weird.
# I doubt whether two reviewers appreciate these two models of person's scores as variables in a way, and as observations in other way.

# Anyway, for a fixed time T, equation (1) will be m = c + X W + epsilon, 
# where epsilon is an error term of iid N(0, sigma^2). It is not residual but error term.
# By regression, hat m is estimated, and residual is obtained as m - hat m.
# And then, using residual, the assumptions of error term will be investigated.

# Here, the terminology "longitudinal profile" is very confusing. 
# Equation (1) is a model for a fixed time T, and reviewer misunderstands this situation.
# Title "Modeling Longitudinal Person Profiles with PAMS in R" means "Modeling Person Profiles with PAMS in R for a time T"

# With respect to multiple regression analysis, we have 18 observations with three independent variables,
# which cannot be statistically validated by data.
# Thus, it is not sufficient to conduct regression analysis and diagnostics, and to deduce any conclusion.
# For example, in case of #100 person, I do not think that Q-Q plot and residual plot provide
# any informative results. However, if based on these plots, one suspects the normality and homogeneous variance of error term.
# If just running the statistical test, the null hypothesis of the normality and homogeneous variance will not 
# be rejected under certain significant level.

# Please run below code, and see the results.

profileOri <- read.csv(file="profileMDSasindependent.csv", header=TRUE) 
nsubject <- nrow(testdata)

# After installing "lmtest" package
library(lmtest)

testres <- NULL
for (i in 1:nsubject) {
  y <- as.numeric(testdata[i,])
  my <- mean(y)
  y <- y - my
  out <- lm(y ~ -1 + Ori1 + Ori2 + Ori3, data.frame(y, profileOri))
  tmp <- shapiro.test(out$residuals)$p.value
  tmp1 <- bptest(y ~ -1 + Ori1 + Ori2 + Ori3, data=data.frame(y, profileOri))$p.value
  testres <- rbind(testres, c(tmp, tmp1))
}
round(testres, 3)

# The first column for the normality test of H_0 : it is normal vs  H_1 : it is not normal, by Shapiro-Wilks.
# A p-value larger than alpha means that H_0 will not be rejected.

# The second column for the constant variance test of H_0 : it is homogeneous vs  H_1 : it is heterogeneous, by Breusch-Pagan test.
# A p-value larger than alpha means that H_0 will not be rejected.

# if the significant level set as alpha,
ahpla <- 0.05
round(testres, 3) < ahpla 
# FALSE means that H_O is not rejected under alpha.

# Q-Q plot for residuals and residual plot of predicted vs residuals
par(mfrow=c(1,2))

# for the i-th person
i <- 100
y <- as.numeric(testdata[i,])
my <- mean(y)
y <- y - my
out <- lm(y ~ -1 + Ori1 + Ori2 + Ori3, data.frame(y, profileOri))

qqnorm(out$residuals); qqline(out$residuals)
plot(out$fitted.values, out$residuals, xlab="predicted", ylab="residual", main="residual plot"); abline(h=0)


