library(smacof)
setwd("/Users/sekangkim/Documents/J Sub/PAMSproject/") 

#### Load the R function 
source("BootSmacof.R") # load the R function

#### Reading data
testdata <- read.csv(file="Longitudinal.csv", header=TRUE) 
colnames(testdata) <- c("E1_1","E1_2","E1_3","E1_4","E1_5","E1_6","E1_7","E1_8","E1_9","E1_10","E1_11","E3_1","E3_2","E3_3","E3_4","E3_5","E3_6","E3_7","E3_8","E3_9","E3_10","E3_11")
testname <- c("E1_1","E1_2","E1_3","E1_4","E1_5","E1_6","E1_7","E1_8","E1_9","E1_10","E1_11","E3_1","E3_2","E3_3","E3_4","E3_5","E3_6","E3_7","E3_8","E3_9","E3_10","E3_11")
#testname <- c("OV1","NS2","VA3","LP4","PP5","SR6","VS7","GI8","CF9","NR10","NP11","NW12","VA13","PR14","AS15","ON16","PC17","MW18")

#### 
# When you conduct the normality test, would you conduct Time 1 (E1_1 – E1_11) and Time 2 (E3_1 –E3_11) separately 
# and then test the Time 1 and Time 2 combined data?

testdata1 <- testdata[, 1:11]
testdata2 <- testdata[, 12:22]

# I just used "MVN" library for multivariate normality test, which I used first time.
# There are five type of tests of H_O : Multivariate Normality vs H_1 : It is not.
# Given significant level alpha, if p-value is smaller alpha, then H_O will be rejected.
# For the Longitudinal data, all cases are rejected with small p-value, which means it is not multivariate normal.

# After installing "MVN" package
library(MVN)

#### Multivariate Normality Test for all data
# H_O : Multivariate Normality vs H_1 : It is not.

# Mardia’s MVN test
result <- mvn(data = testdata, mvnTest = "mardia")
result$multivariateNormality
#Test              Statistic          p value Result
#1 Mardia Skewness 18210.4865944201       0     NO
#2 Mardia Kurtosis   100.36978316286      0     NO

# Henze-Zirkler’s MVN test
result <- mvn(data = testdata, mvnTest = "hz")
result$multivariateNormality
#Test      HZ p value MVN
#1 Henze-Zirkler 1.32284       0  NO

# Royston’s MVN test
result <- mvn(data = testdata, mvnTest = "royston")
result$multivariateNormality
# Test        H p value MVN
# 1 Royston 2872.982       0  NO

# Doornik-Hansen’s MVN test
result <- mvn(data = testdata, mvnTest = "dh")
result$multivariateNormality
#Test        E df      p value MVN
#1 Doornik-Hansen 369.1748 44 5.871267e-53  NO

# Energy test
result <- mvn(data = testdata, mvnTest = "energy")
result$multivariateNormality
#Test Statistic p value MVN
#1 E-statistic  19.02847       0  NO

#### Multivariate Normality Test for Time 1
# H_O : Multivariate Normality vs H_1 : It is not.

# Mardia’s MVN test
result <- mvn(data = testdata1, mvnTest = "mardia")
result$multivariateNormality
#Test       Statistic p value Result
#1 Mardia Skewness 5645.0421163708       0     NO
#2 Mardia Kurtosis 52.463936119043       0     NO
#3             MVN            <NA>    <NA>     NO

# Henze-Zirkler’s MVN test
result <- mvn(data = testdata1, mvnTest = "hz")
result$multivariateNormality
#Test       HZ p value MVN
#1 Henze-Zirkler 2.675676       0  NO

# Royston’s MVN test
result <- mvn(data = testdata1, mvnTest = "royston")
result$multivariateNormality
#Test        H       p value MVN
#1 Royston 1255.793 3.918086e-262  NO

# Doornik-Hansen’s MVN test
result <- mvn(data = testdata1, mvnTest = "dh")
result$multivariateNormality
#Test        E df      p value MVN
#1 Doornik-Hansen 245.6256 22 1.077223e-39  NO

# Energy test
result <- mvn(data = testdata1, mvnTest = "energy")
result$multivariateNormality
#Test Statistic p value MVN
#1 E-statistic  16.75876       0  NO

#### Multivariate Normality Test for Time 2
# H_O : Multivariate Normality vs H_1 : It is not.

# Mardia’s MVN test
result <- mvn(data = testdata2, mvnTest = "mardia")
result$multivariateNormality
#Test        Statistic p value Result
#1 Mardia Skewness 9816.36261767864       0     NO
#2 Mardia Kurtosis 106.438959374944       0     NO
#3             MVN             <NA>    <NA>     NO

# Henze-Zirkler’s MVN test
result <- mvn(data = testdata2, mvnTest = "hz")
result$multivariateNormality
#Test       HZ p value MVN
#1 Henze-Zirkler 8.183358       0  NO

# Royston’s MVN test
result <- mvn(data = testdata2, mvnTest = "royston")
result$multivariateNormality
#Test        H       p value MVN
#1 Royston 1539.341 9.881313e-324  NO

# Doornik-Hansen’s MVN test
result <- mvn(data = testdata2, mvnTest = "dh")
result$multivariateNormality
#Test       E df      p value MVN
#1 Doornik-Hansen 265.846 22 9.598518e-44  NO

# Energy test
result <- mvn(data = testdata2, mvnTest = "energy")
result$multivariateNormality
#Test Statistic p value MVN
#1 E-statistic  43.83538       0  NO

