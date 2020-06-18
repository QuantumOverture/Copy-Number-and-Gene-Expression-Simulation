#set.seed(1)
# Probes -> sequence of dna used measure copy number(in/outfProbes) (spot in microarray)



# Formulae:
# M(c) = log2 [(c · Pt + 2 · (1 − Pt))/2], where Pt = U(0.3, 0.7)
# N(µ = M(2), σ^{2} = S^{2}), where S = U(min = 0.1, max = 0.3) --> Baseline

# Setup functions and values:
CellAdmixture <- function(value){
  Pt <- runif(1,0.3,0.7)
  Result <- log2( (value * Pt + 2 *(1-Pt))/2 )
  return(Result)
}
Variance <- function(value){
  return(runif(1,0.1,0.3)**2)
}




# Copy Number Simulation Types

NormalCell <- function(NumberOfProbes){
return(c(rnorm(NumberOfProbes,CellAdmixture(2),Variance())))
}

GainNarrowMediumAberration1 <- function(NumberOfProbes){
return(c(rnorm(NumberOfProbes,CellAdmixture(6),Variance())))
}

# Get Copy number
GetCopyNumber <- function(Value){
  #return(round(2*(2**Value),0))
  return(2*(2**Value))
}


# Gene Expression Simulation (without adding True/False positives)
LinearExpression <- function(Vector){
  TempVector <- c()
  for (i in Vector){
    CurrCopyNum <- GetCopyNumber(i)
    CurrValue <- 2*CurrCopyNum+2
    if (CurrCopyNum >= 3 | CurrCopyNum == 0 | CurrCopyNum == 1){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else{
      TempVector <- c(TempVector,CurrValue)
    }
  }
  print(TempVector)
  return(TempVector)
}

SigmoidExpression <- function(Vector){
  TempVector <- c()
  for (i in Vector){
    CurrCopyNum <- GetCopyNumber(i)
    el = (CurrCopyNum/2) + 2
    CurrValue <- 8/(1+exp(-1*el))
    if (CurrCopyNum >= 3 | CurrCopyNum == 0 | CurrCopyNum == 1){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else{
      TempVector <- c(TempVector,CurrValue)
    }
  }
  print(TempVector)
  return(TempVector)
}


StepWiseExpression <- function(Vector){
  TempVector <- c()
  for (i in Vector){
    CurrCopyNum <- GetCopyNumber(i)

    if (CurrCopyNum < 1){
      CurrValue <- 2
    }else if (CurrCopyNum <= 2 & CurrCopyNum >= 1){
      CurrValue <- 6
    }else if (CurrCopyNum < 4 & CurrCopyNum > 2){
      CurrValue <- 10
    }else if (CurrCopyNum <= 8 & CurrCopyNum >= 4){
      CurrValue <- 12
    }else{
      CurrValue <- 14
    }


    if (CurrCopyNum >= 3 | CurrCopyNum == 0 | CurrCopyNum == 1){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else{
      TempVector <- c(TempVector,CurrValue)
    }
  }
  print(TempVector)
  return(TempVector)
}

NormalExpression <- function(NumberOfProbes){
  return(c(rnorm(NumberOfProbes,6,sqrt(2))))
}




# https://thepracticalr.wordpress.com/2016/08/30/2-y-axis-plotting/
# https://rpubs.com/riazakhan94/297778
# https://stackoverflow.com/questions/38247907/how-to-set-the-y-range-in-boxplot-graph
plot(c(NormalCell(1000),GainNarrowMediumAberration1(5000),NormalCell(1000)),col = "blue",  ylim=c(-2,3))
par(new = TRUE)
plot(NormalExpression(7000),xaxt = "n", yaxt = "n",
     ylab = "", xlab = "", col = "gray", lty = 2)
axis(side = 4)
