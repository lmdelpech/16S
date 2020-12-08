#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: Alban Ramette
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# function to repeatedly subsample Sample by OTU table to obtain more reliable alpha diversity estimates
# calculates number of OTU, chao1, inverse Simpson, absolute and relative singletons, absolute doubletons
#
# input: 
# Data - sample by OTU table (samples are columns)
# n - number of iterations (how many times do you want to subsample)
# sub - library size (how many sequences per sample)
# subTable - should one iteration of rarefied OTU be writen to file (logical)
#
# output:
# iterations - list with values (number of OTU, chao1, inverse Simpson, percentage of OTUs that are
#   absolute and relative singletons, and absolute doubletons) for each iteration
# summaryAlpha - mean of iterations for above indices of rarefyied data set
# summaryHill - mean of iterations for Hill numbers 0, 1, 2 of rarefyied data set
# summaryHillRaw - Hill numbers 0, 1, 2 of original data set
#  
# 
# dependencies:
#  require(vegan)
#=============================================================================
#documentation end

SubsampleNGS <- function(Data, n, sub, subTable = F) {
  
  if(!"vegan" %in% installed.packages()) {
    install.packages("vegan")
  }
  require(vegan)
  
  #removing samples with less than sub sequences
  Data <- Data[, colSums(Data) >= sub]
  
  output <- list(iterations = list(nOTU = matrix(NA, n, ncol(Data)),
                                   chao1 = matrix(NA, n, ncol(Data)),
                                   ace = matrix(NA, n, ncol(Data)),
                                   invS = matrix(NA, n, ncol(Data)),
                                   shannon = matrix(NA, n, ncol(Data)),
                                   SSOabs = matrix(NA, n, ncol(Data)),
                                   SSOrel = matrix(NA, n, ncol(Data)),
                                   DSOabs = matrix(NA, n, ncol(Data))),
                 summaryAlpha = matrix(NA, 8, ncol(Data)),
                 summaryHill = matrix(NA, 3, ncol(Data)),
                 summaryHillRaw = matrix(NA, 3, ncol(Data))
  )
  
  colnames(output$iterations$nOTU) = 
    colnames(output$iterations$chao1) = 
    colnames(output$iterations$ace) = 
    colnames(output$iterations$invS) = 
    colnames(output$iterations$shannon) = 
    colnames(output$iterations$SSOabs) = 
    colnames(output$iterations$SSOrel) = 
    colnames(output$iterations$DSOabs) = 
    colnames(output$summaryAlpha) = 
    colnames(output$summaryHill) = 
    colnames(output$summaryHillRaw) = 
    colnames(Data)
  
  rownames(output$summaryAlpha) <- names(output$iterations)
  rownames(output$summaryHill) <- c("Hill0","Hill1","Hill2")
  rownames(output$summaryHillRaw) <- c("Hill0","Hill1","Hill2")
  
  for(j in 1:n){
    
    print(j)    
    
    #subsampling 
    SampleOTU <-  t(rrarefy(t(Data), sub))
    
    #exclude empty OTUs
    Data_new <- SampleOTU[rowSums(SampleOTU) > 0, ]
    
    if (subTable == T & j == 1) {
      write.table(Data_new, "subTable.txt", quote = F, sep = "\t")
    }
    
    #nOTU
    Data01 <- Data_new
    Data01[Data01 > 0] <- 1
    nOTU <- apply(Data01, 2, sum)
    output$iterations$nOTU[j, ] <- nOTU
    
    #chao1 + ACE
    chao1 <- estimateR(t(Data_new))
    output$iterations$chao1[j, ] <- chao1[2, ]
    output$iterations$ace[j, ] <- chao1[4, ]
    
    #invS + shannon
    invS <- diversity(t(Data_new), "inv")
    output$iterations$invS[j, ] <- invS
    Shannon <- diversity(t(Data_new), "shannon")
    output$iterations$shannon[j, ] <- Shannon
    
    #SSO
    D <- t(Data_new)
    
    if (sum(apply(D, 2, sum) == 1) >= 1) {
      SSOabs.D <- data.frame(D[, apply(D, 2, sum) == 1])
      SSOabsPerSample <- apply(SSOabs.D, 1, sum)
    } else {
      SSOabsPerSample <- rep(0, nrow(D))
    }
    output$iterations$SSOabs[j, ] <- round((SSOabsPerSample / nOTU) * 100, 2)
    
    
    #DSOabs
    if (sum(apply(D, 2, sum) == 2) >= 1) {
      DSOabs.D <- data.frame(D[, apply(D, 2, sum) == 2])
      DSOabsPerSample <- apply(DSOabs.D, 1, sum)
    } else {
      DSOabsPerSample <- rep(0, nrow(D))
    }
    output$iterations$DSOabs[j, ] <- round((DSOabsPerSample / nOTU) * 100, 2)
    
    #SSOrel
    # removing SSOabs
    D1 <- D[, apply(D, 2, sum) > 1] 
    # SSOrel should have at least a sample with one sequence alone. See definition.
    if (sum(apply(D1, 2, function(x) any(x == 1)) == TRUE) >= 1) {
      SSOrel.D <- data.frame(D1[, apply(D1, 2, function(x) any(x == 1)) == TRUE])
      # SSOrel.D contains the subtable with all SSOrel of the dataset
      SSOrelPerSample <- apply(SSOrel.D, 1, function(x) sum(x == 1))  
    } else {
      SSOrelPerSample <- rep(0, nrow(D))
    }
    
    output$iterations$SSOrel[j, ] <- round((SSOrelPerSample / nOTU) * 100, 2)
  }
  
  for (i in 1:8) {
    output$summaryAlpha[i,] <- apply(output$iterations[[i]], 2, mean)
  }
  
  output$summaryHill[1, ] <- output$summaryAlpha["nOTU",]
  output$summaryHill[2, ] <- exp(output$summaryAlpha["shannon",])
  output$summaryHill[3, ] <- output$summaryAlpha["invS",]
  
  output$summaryHillRaw[1, ] <- colSums(decostand(Data, method = "pa"))
  output$summaryHillRaw[2, ] <- exp(diversity(t(Data), "shannon"))
  output$summaryHillRaw[3, ] <- diversity(t(Data), "inv")
  
  return(output) 
}