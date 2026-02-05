if (!require(profileR)) install.packages("profileR")

library(profileR)

testdata <- read.csv(file="Cross-sectional.csv", header=FALSE) 

colnames(testdata) <- c("OV1","NS2","VA3","LP4","PP5","SR6","VS7","GI8","CF9","NR10","NP11","NW12","VA13","PR14","AS15","ON16","PC17","MW18")
testname <- c("OV1","NS2","VA3","LP4","PP5","SR6","VS7","GI8","CF9","NR10","NP11","NW12","VA13","PR14","AS15","ON16","PC17","MW18")

for (j in seq_along(testdata)) {
  if (!is.numeric(testdata[[j]])) {
    testdata[[j]] <- as.numeric(gsub(",", ".", as.character(testdata[[j]])))
  }
}

na_rows <- which(!complete.cases(testdata))
if (length(na_rows) > 0) {
  message("Removing ", length(na_rows), " rows with NA values.")
  testdata <- testdata[complete.cases(testdata), ]
}

stopifnot(is.data.frame(testdata))
stopifnot(ncol(testdata) == 18)        # 18 subscales expected
stopifnot(all(sapply(testdata, is.numeric)))


result_original <- pams(testdata, dim=5)
str(result_original)
print(result_original)


# See what components are present
names(result_original)

# If present, examine dimensional coordinates
if ("dimensional.configuration" %in% names(result_original)) {
  head(result_original$dimensional.configuration)
}

# Plot the first 3 core-profiles (if present)
if ("dimensional.configuration" %in% names(result_original)) {
  par(mfrow = c(1,3))
  for (i in 1:min(3, ncol(result_original$dimensional.configuration))) {
    plot(scale(result_original$dimensional.configuration[,i]), type="b",
         main=paste("CP", i), xlab="", ylab="Coordinate")
    abline(h=0, lty=3)
  }
}



summary(result_original$weights.matrix[, "R.sq"])

hist(result_original$weights.matrix[, "R.sq"], main="R² distribution")



fits <- list()
for (d in 2:7) {
  fits[[d]] <- pams(testdata, dim=d)
  cat("dim=", d, " mean R² = ", mean(fits[[d]]$weights.matrix[, "R.sq"]), "\n")
}





# BootSmacof.R
# A self-contained implementation of the bootstrap-PAMS procedure
# Produces aligned bootstrap distributions and BCa CIs for coordinates.
# Requirements: smacof, boot

# Usage example:
# source("BootSmacof.R")
# empprofile <- BootSmacof(testdata = testdata, nprofile = 5, nBoot = 2000)
BootSmacof_safe <- function(testdata, participant = NULL, 
                            mds = "smacof", type = "ordinal", 
                            distance = "euclid", scale = FALSE, 
                            nprofile = 3, direction = rep(1, nprofile),
                            cl = 0.95, nBoot = 200, testname = NULL, 
                            file = NULL) {
  
  n <- nrow(testdata)
  p <- ncol(testdata)
  
  # If participants not specified, use all
  if (is.null(participant)) participant <- 1:n
  
  # Storage for bootstrap results
  MDSprofile <- vector("list", nprofile)
  MDSsummary <- vector("list", nprofile)
  
  for (i in 1:nprofile) {
    MDSprofile[[i]] <- matrix(NA, nrow = length(participant), ncol = nBoot)
    MDSsummary[[i]] <- data.frame(Ori = rep(NA, p),
                                  Mean = rep(NA, p),
                                  SE = rep(NA, p),
                                  Lower = rep(NA, p),
                                  Upper = rep(NA, p),
                                  BCaLower = rep(NA, p),
                                  BCaUpper = rep(NA, p))
  }
  
  # Store weights matrix
  Weight <- matrix(NA, length(participant), 2 * nprofile + 2)
  colnames(Weight) <- c(paste0("w", seq_len(nprofile)), "level", "R2", paste0("corDim", seq_len(nprofile)))
  
  # Compute core profiles for original data
  if (mds == "smacof") {
    require(smacof)
    fit_orig <- smacofSym(dist(t(testdata)), nprofile, type = type)
    DimConf <- fit_orig$conf
  } else {
    DimConf <- cmdscale(dist(t(testdata)), k = nprofile)
  }
  
  # Apply direction flips
  for (i in 1:nprofile) DimConf[,i] <- DimConf[,i] * direction[i]
  
  # Bootstrap loop
  for (b in 1:nBoot) {
    idx <- sample(1:n, replace = TRUE)
    tmp_data <- testdata[idx, , drop = FALSE]
    
    # Try to compute MDS; if fails, skip this bootstrap
    tmp_fit <- try({
      if (mds == "smacof") {
        smacofSym(dist(t(tmp_data)), nprofile, type = type)$conf
      } else {
        cmdscale(dist(t(tmp_data)), k = nprofile)
      }
    }, silent = TRUE)
    
    if (inherits(tmp_fit, "try-error")) next
    
    # Apply direction flips
    for (i in 1:nprofile) tmp_fit[,i] <- tmp_fit[,i] * direction[i]
    
    # Fill bootstrap matrix safely
    expected_cols <- 2 * nprofile + 2
    tmp <- matrix(NA, nrow = length(participant), ncol = expected_cols)
    cols_to_fill <- min(ncol(tmp_fit), expected_cols)
    tmp[, 1:cols_to_fill] <- tmp_fit[participant, 1:cols_to_fill, drop = FALSE]
    
    colnames(tmp) <- c(paste0("w", seq_len(nprofile)), "level", "R2", paste0("corDim", seq_len(nprofile)))
    
    Weight <- tmp  # store or accumulate depending on what you need
  }
  
  return(list(
    MDS = list(conf = DimConf),
    MDSprofile = MDSprofile,
    MDSsummary = MDSsummary,
    Weight = Weight,
    nprofile = nprofile,
    nBoot = nBoot
  ))
}


fit_main <- pams(testdata, dim = 5)

direction <- sign(fit_main$weights.matrix[,1])  # or set manually




boot_results <- BootSmacof(
  testdata = testdata,
  mds = "smacof",
  type = "ordinal",
  distance = "euclid",
  scale = FALSE,
  nprofile = 5,
  direction = direction,
  cl = 0.95,
  nBoot = 1000
)


# main PAMS configuration
main_conf <- fit_main$dimensional.configuration  # 18 x 5

# bootstrap mean profiles: take mean across bootstrap samples (dimension 1)
boot_means <- apply(boot_results$MDSprofile, c(2,3), mean)  # dims: items x dimensions = 18 x 5

# Now flip dimensions if correlation with main_conf is negative
flip <- sapply(1:5, function(i) cor(main_conf[,i], boot_means[,i]) < 0)
boot_means[, flip] <- -boot_means[, flip]

# Compute correlations after flipping
cor_vals <- sapply(1:5, function(i) cor(main_conf[,i], boot_means[,i]))
cor_vals


# Plot main vs bootstrap means for each dimension
par(mfrow=c(1,5), mar=c(4,4,2,1))
for (i in 1:5) {
  plot(main_conf[,i], boot_means[,i],
       main = paste("Dim", i),
       xlab = "Main PAMS", ylab = "Bootstrap mean",
       pch = 19)
  abline(0,1, col="red", lty=2)
}


# Optional: visualize differences
diffs <- main_conf - boot_means
par(mfrow=c(1,1))
matplot(diffs, type="l", lty=1, main="Differences: Main PAMS - Bootstrap Means",
        ylab="Difference", xlab="Item")
abline(h=0, col="red", lty=2)







# Combine all dimensions into one data frame
ci_df <- do.call(rbind, lapply(1:5, function(d) {
  df <- boot_results$MDSsummary[[d]]
  df <- as.data.frame(df)
  df$Item <- rownames(df)
  df$Dimension <- paste0("Dimension", d)
  df[, c("Item","Dimension","Mean","Lower","Upper")]
}))

# Plot
ggplot(ci_df, aes(x=Item, y=Mean)) +
  geom_point(color="blue") +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.2, color="blue") +
  facet_wrap(~Dimension, scales="free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
