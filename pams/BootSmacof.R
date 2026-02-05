BootSmacof <- function(testdata, participant=NULL, mds=c("smacof", "classical"),  type=c("ratio", "interval", "ordinal", "mspline"), distance=c("euclid", "sqeuclid"), scale=FALSE,
                     nprofile=3, direction=rep(1, nprofile), cl=0.95, nBoot=2000, testname=NULL, file=NULL)
{ 
    ntest <- ncol(testdata) 
    nsubject <- nrow(testdata)
    lalpha <- (1 - cl)/2  
    ualpha <- 1 - lalpha
    
    if (is.null(testname)) testname <- paste0("T", 1:ntest)
    if (ntest < nprofile) stop("The number of profiles must be less than the number of test.")
    if (length(direction) != nprofile) stop("The number of profile direction must be the same to the number of profile.")
    if (!(any(distance == c("euclid", "sqeuclid"))))
        stop("Distance for smacof must be one of euclid, sqeuclid.")
    
    if (scale) testdata <- scale(testdata)
    if (distance == "euclid") 
        distance0 <- dist(t(testdata))
    else if (distance == "sqeuclid") 
      distance0 <- dist(t(testdata))^2     

    # profiles by smacof for original sample
    if (mds == "smacof") {
       MDS <- smacofSym(distance0, ndim = nprofile, type=type)
       for (i in 1:nprofile)
           MDS$conf[,i] <- MDS$conf[,i] * direction[i]      
       profileOri <- MDS$conf
       stressOri <- MDS$stress
    }

    # profiles by classical MDS for original sample
    if (mds == "classical") {
        MDS <- NULL
        MDS$conf <- cmdscale(distance0, k = nprofile)
        MDS$stress <- stressOri <- NULL
        colnames(MDS$conf) <- paste0("D", 1:nprofile)
        for (i in 1:nprofile)
            MDS$conf[,i] <- MDS$conf[,i] * direction[i]  
        profileOri <- MDS$conf
    }
                          
    # empirical distribution of profiles by bootstrapping
    profileBoot <- vector("list", nprofile); stressBoot <- NULL
    for (i in 1:nBoot) {
        # bootstrap sample
        testdataBoot <- testdata[sample(1:nsubject, replace=TRUE), ]
        if (distance == "euclid") 
           distanceB <- dist(t(testdataBoot))
        else if (distance == "sqeuclid")
           distanceB <- dist(t(testdataBoot))^2
        
        # profiles for bootstrap sample
        if (mds == "smacof") {
           tmp0 <- smacofSym(distanceB, ndim = nprofile, type=type)
           tmp1 <- tmp0$conf
           stressBoot <- c(stressBoot, tmp0$stress)
        } else
           tmp1 <- cmdscale(distanceB, k = nprofile)

        tmp1 <- tmp1 %*% diag(sign(diag(cor(tmp1, profileOri)))) 
        for (j in 1:nprofile)
            profileBoot[[j]] <- rbind(profileBoot[[j]], tmp1[,j]) 
    }

    # BCa part is based on "bootBCa" R package Version 1.0
    # Description for "bootBCa" R package Version 1.0
    # Author: S original, from StatLib, by Rob Tibshirani. R port by
    # Friedrich Leisch. Enhancements by David Flater.
    # License: BSD_3_clause + file LICENSE
    
    profileu <- uu <- vector("list", nprofile) 
    stressu <- suu <- NULL
    for (i in 1:nsubject) {
        
        if (distance == "euclid") 
            distance0 <- dist(t(testdata[-i, ]))
        else if (distance == "sqeuclid") 
            distance0 <- dist(t(testdata[-i, ]))^2     
        
        # profiles by smacof 
        if (mds == "smacof") {
            tmp <- smacofSym(distance0, ndim = nprofile, type=type)
            for (i in 1:nprofile)
                tmp$conf[,i] <- tmp$conf[,i] * direction[i]      
            stressu <- c(stressu, tmp$stress)
        }
        
        # profiles by classical MDS 
        if (mds == "classical") {
            tmp <- NULL
            tmp$conf <- cmdscale(distance0, k = nprofile)
            for (i in 1:nprofile)
                tmp$conf[,i] <- tmp$conf[,i] * direction[i]  
        }
        
        for (j in 1:nprofile)
            profileu[[j]] <- rbind(profileu[[j]], tmp$conf[,j]) 
    }   
    
    acc <- NULL
    for (j in 1:nprofile) {
        uu[[j]] <- -sweep(profileu[[j]], 2, apply(profileu[[j]], 2, mean))
        acc <- cbind(acc, apply(uu[[j]], 2, function(x) sum(x*x*x)/(6*(sum(x*x))^1.5)))
    }

    zalpha <- qnorm(c(lalpha, ualpha))
    BCaCI <- vector("list", nprofile)
    for (i in 1:nprofile) {  
        tmp <- NULL
        for(j in 1:ntest) {
            z0 <- qnorm(sum(profileBoot[[i]][,j] < profileOri[,i][j]) / nBoot) 
            tt <- pnorm(z0 + (z0 + zalpha) / (1 - acc[j, i] * (z0 + zalpha))) 
            tmp1 <- quantile(profileBoot[[i]][, j], probs=tt)
            tmp <- rbind(tmp, tmp1)
        }
        BCaCI[[i]] <- tmp
    }

    # Summary Statistics for Profile
    profile <- vector("list", nprofile)
    for (i in 1:nprofile) 
        profile[[i]] <- data.frame(Ori=profileOri[,i], Mean=apply(profileBoot[[i]], 2, mean),
                                   SE=apply(profileBoot[[i]], 2, sd), 
                                   Lower=apply(profileBoot[[i]], 2, quantile, lalpha),
                                   Upper=apply(profileBoot[[i]], 2, quantile, ualpha), 
                                   BCaLower=BCaCI[[i]][,1],
                                   BCaUpper=BCaCI[[i]][,2], row.names=testname)
    
    tmp2 <- profile[[1]]
    if (nprofile > 2)
        for (i in 2:nprofile) 
            tmp2 <- cbind(tmp2, profile[[i]])

    # Summary Statistics for stress
    stresssummary <- NULL
    if (mds == "smacof") {
        suu <- -(stressu - mean(stressu))
        sacc <- sum(suu*suu*suu)/(6*(sum(suu*suu))^1.5)
        z0 <- qnorm(sum(stressBoot < stressOri) / nBoot) 
        tt <- pnorm(z0 + (z0 + zalpha) / (1 - sacc * (z0 + zalpha))) 
        BCaCIstress <- quantile(stressBoot, probs=tt)     
        stresssummary <- data.frame(Ori=stressOri, Mean=mean(stressBoot),
                                    SE=sd(stressBoot), 
                                    Lower=quantile(stressBoot, lalpha),
                                    Upper=quantile(stressBoot, ualpha), 
                                    BCaLower=BCaCIstress[1],
                                    BCaUpper=BCaCIstress[2], row.names=NULL)
    }

    if (!(is.null(file))) {
        cat("Summary Statistics for Stress", file=paste0(file, "MDS.csv"), sep="", "\n")  
        cat(c("Ori", "Mean", "SE", "Lower", "Upper", "BCaLower", "BCaUpper"), file=paste0(file, "MDS.csv"), sep=",", "\n", append=TRUE)
        if (mds == "classical") {
            cat(stressOri, file=paste0(file, "MDS.csv"), sep=",", "\n\n", append=TRUE)             
        }
        if (mds == "smacof") {
            cat(as.numeric(stresssummary), file=paste0(file, "MDS.csv"), sep=",", "\n\n", append=TRUE)           
        }
        
        cat("Summary Statistics for Profile", file=paste0(file, "MDS.csv"), sep="", "\n", append=TRUE)  
        cat(c("Name", paste0(rep(c("Ori", "Mean", "SE", "Lower", "Upper", "BCaLower", "BCaUpper"), nprofile), rep(1:nprofile, each=7))), 
            file=paste0(file, "MDS.csv"), sep=",", "\n", append=TRUE)
        write.table(tmp2, file=paste0(file, "MDS.csv"), sep=",", append=TRUE, col.names=FALSE)    
    }


    # Regression
    formulaall <- as.formula(paste("y ~ -1 + ", paste(paste0("D", 1:nprofile), collapse= "+")))
    
    formulaeach <- NULL
    for (k in 1:nprofile)
      formulaeach <- c(formulaeach, paste(paste0("D", k, " ~ -1+"), paste(paste0("D", setdiff(1:nprofile, k), collapse= "+"))))

    R2 <- rep(0, nprofile)
    outDresiduals <- NULL
    for (j in 1:nprofile) {   
      outD <- lm(formulaeach[[j]], data.frame(profileOri))
      outDresiduals <- cbind(outDresiduals, outD$residuals)
      R2[j] <- summary(outD)$r.squared
    }    
    
    result <- NULL
    pcorr <- rep(0, nprofile)
    for (i in 1:nsubject) {
      y <- as.numeric(testdata[i,])
      my <- mean(y)
      y <- y - my
      out <- lm(formulaall, data.frame(y, profileOri))
      coeffmulti <- out$coefficients
      R2multi <- summary(out)$r.squared 
      
      for (j in 1:nprofile) {   
        outy <- update(out, as.formula(paste0(". ~ . -D", j)))
        #outD <- lm(formulaeach[[j]], data.frame(y, profileOri))
        #outr <- lm(outy$residuals ~ -1 + outD$residuals, data.frame(outy$residuals, outD$residuals))
        pcorr[j] <- cor(outy$residuals, outDresiduals[,j])
      }
      
      result <- rbind(result, c(coeffmulti, my, R2multi, pcorr))
    }
    meanR2 <- mean(result[, nprofile + 2])
    
    colnames(result) <- c(paste0("w", 1:nprofile), "level", "R^2", paste0("corDim", 1:nprofile))
    rownames(result) <- paste0("#", 1:nsubject) 
    if (!(is.null(file))) {
      cat(c("Overall R^2", meanR2), file=paste0(file, "Weight.csv"), sep=",", "\n")  
      cat(c("id", paste0("w", 1:nprofile), "level", "R^2", paste0("corDim", 1:nprofile)), file=paste0(file, "Weight.csv"), sep=",", "\n",  append=TRUE)
      write.table(result, file=paste0(file, "Weight.csv"), sep=",", append=TRUE, col.names=FALSE)    
    }
    
    # Regression Bootstrapping
    formulaallB <- as.formula(paste("y ~ -1 + ", paste(paste0("X", 1:nprofile), collapse= "+")))
    formulaeachB <- NULL
    for (k in 1:nprofile)
      formulaeachB <- c(formulaeachB, paste(paste0("X", k, " ~ -1+"), paste(paste0("X", setdiff(1:nprofile, k), collapse= "+"))))
    
    resultB <- resultBP <- NULL   
    if (!is.null(participant)) {
      sumstat <- rep(0, 5*nprofile)
      index1 <- seq(1, 5*nprofile, by=5)
      index2 <- seq(2, 5*nprofile, by=5)
      index3 <- seq(3, 5*nprofile, by=5)
      index4 <- seq(4, 5*nprofile, by=5)
      index5 <- seq(5, 5*nprofile, by=5)
      
      xBoot <- vector("list", nBoot)
      for (k in 1:nBoot) 
        for (j in 1:nprofile) 
          xBoot[[k]] <- cbind(xBoot[[k]], profileBoot[[j]][k,])
 
      outDresiduals <- vector("list", nBoot)
      for (k in 1:nBoot)       
        for (j in 1:nprofile) {   
          outD <- lm(formulaeachB[[j]], data.frame(y, xBoot[[k]]))
          outDresiduals[[k]] <- cbind(outDresiduals[[k]], outD$residuals)
        } 
      
      for (i in 1:length(participant)) {
        y <- as.numeric(testdata[participant[i],])
        my <- mean(y)
        y <- y - my
        
        tmpweight <- tmppcorr <- NULL
        pcorr <- rep(0, nprofile)
        for (k in 1:nBoot) {
          tmplm <- lm(formulaallB, data.frame(y, xBoot[[k]]))
          tmpweight <- rbind(tmpweight, coefficients(tmplm))
 
          for (j in 1:nprofile) {   
            outy <- update(tmplm, as.formula(paste0(". ~ . -X", j)))     
            pcorr[j] <- cor(outy$residuals, outDresiduals[[k]][,j])  
          }
          tmppcorr <- rbind(tmppcorr, pcorr)
        }
        
        sumstat[index1] <- result[participant[i], 1:nprofile] 
        sumstat[index2] <- apply(tmpweight, 2, mean)
        sumstat[index3] <- apply(tmpweight, 2, sd)
        sumstat[index4] <- apply(tmpweight, 2, quantile, lalpha)
        sumstat[index5] <- apply(tmpweight, 2, quantile, ualpha)     
        resultB <- rbind(resultB, c(sumstat, result[participant[i], -(1:nprofile)]))  
        
        sumstat[index1] <- result[participant[i], -(1:(nprofile+2))] 
        sumstat[index2] <- apply(tmppcorr, 2, mean)
        sumstat[index3] <- apply(tmppcorr, 2, sd)
        sumstat[index4] <- apply(tmppcorr, 2, quantile, lalpha)
        sumstat[index5] <- apply(tmppcorr, 2, quantile, ualpha)    
        resultBP <- rbind(resultBP, sumstat)
      }
      colnames(resultB) <- c(paste0(rep(c("w", "m", "se", "L", "U"), nprofile), rep(1:nprofile, each=5)), 
                             "level", "R^2", paste0("corDim", 1:nprofile))
      rownames(resultB) <- paste0("#", participant)
      if (!(is.null(file))) {
        cat(c("id", paste0(rep(c("w", "m", "se", "L", "U"), nprofile), rep(1:nprofile, each=5)), 
              "level", "R^2", paste0("corDim", 1:nprofile)), file=paste0(file, "WeightB.csv"), sep=",", "\n")
        write.table(resultB, file=paste0(file, "WeightB.csv"), sep=",", append=TRUE, col.names=FALSE)    
      }  
      
      colnames(resultBP) <- paste0(rep(c("corDim", "m", "se", "L", "U"), nprofile), rep(1:nprofile, each=5))
      rownames(resultBP) <- paste0("#", participant)
      if (!(is.null(file))) {
        cat(c("id", paste0(rep(c("corDim", "m", "se", "L", "U"), nprofile), rep(1:nprofile, each=5))), file=paste0(file, "CorB.csv"), sep=",", "\n")
        write.table(resultBP, file=paste0(file, "CorB.csv"), sep=",", append=TRUE, col.names=FALSE)    
      }        
    }  
    
    list(MDS=MDS, MDSsummary=profile, MDSprofile=profileBoot, stresssummary=stresssummary, stressprofile=stressBoot, MDSR2=R2, Weight=result, WeightmeanR2=meanR2, WeightB=resultB, PcorrB=resultBP, nprofile=nprofile, nBoot=nBoot, scale=scale, testname=testname)
} 