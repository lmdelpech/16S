#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: 
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# perform pairwise ANOSIM tests
#
# input: 
# M - community matrix
# E - grouping factor
# distance - distance measure to be used in anosim
# padj - p value correction method
# if strata is supplied, the type of permutations for Within() and Plots() also has to be supplied
#
# output:
# list with [[1]] anosim R, [[2]] unadjusted p value, [[3]] adjusted p value in triangular matrix (full)
# 
# dependencies:
# require(vegan)
#=============================================================================
#documentation end


ANOSIMposthoc <- function(
  M,
  E,
  distance = "bray", 
  padj = "fdr",
  strata = NULL,
  type.within = NULL, 
  type.plots = NULL
) {
  E = droplevels(E)
  Mlist = list()
  for(i in 1:length(levels(E))){
    Mlist[[i]] = M[E == levels(E)[i], ]
  }
  
  Elist = list()
  for(i in 1:length(levels(E))){
    Elist[[i]] = as.numeric(E[E == levels(E)[i]])
  }
  
  if(!is.null(strata)) {
    Slist = list()
    for(i in 1:length(levels(E))){
      Slist[[i]] = as.numeric(strata[E == levels(E)[i]])
    }
  }
  
  result = list(
    anosimR = matrix(NA, length(levels(E)), length(levels(E))),
    anosimP = matrix(NA, length(levels(E)), length(levels(E))),
    anosimPadj = matrix(NA, length(levels(E)), length(levels(E))))
  colnames(result$anosimR) = colnames(result$anosimP) = levels(E)
  rownames(result$anosimR) = rownames(result$anosimP) = levels(E)
  for(i in 1:(length(levels(E)) - 1)){
    for(j in (i + 1):length(levels(E))){
      if(is.null(strata)) {
        temp = anosim(
          rbind(Mlist[[i]], Mlist[[j]]),
          c(Elist[[i]], Elist[[j]]),
          distance = distance
        )
      } else {
        temp = anosim(
          rbind(Mlist[[i]], Mlist[[j]]),
          c(Elist[[i]], Elist[[j]]),
          distance = distance,
          permutations = how(
            within = Within(type = type.within),
            plots = Plots(strata = c(Slist[[i]], Slist[[j]]), type = type.plots),
            nperm = 999
          )
        )
      }
      result$anosimR[j,i] = temp$statistic
      result$anosimP[j,i] = temp$signif
    }
  }
  result$anosimPadj = matrix(
    p.adjust(
      as.vector(result$anosimP), 
      method = padj, 
      n = length(which(!is.na(as.vector(result$anosimP))))
    ),
    length(levels(E)),
    length(levels(E))
  )
  colnames(result$anosimPadj) = levels(E)
  rownames(result$anosimPadj) = levels(E)
  return(result)
}


