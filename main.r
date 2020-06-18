#set.seed(1)

# Probes -> sequence of dna used measure copy number(in/out genes) (spot in microarray) => change NumberOfGenes



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

NormalCell <- function(NumberOfGenes){
return(c(rnorm(NumberOfGenes,CellAdmixture(2),Variance())))
}

GainNarrowMediumAberration1 <- function(NumberOfGenes){
return(c(rnorm(NumberOfGenes,CellAdmixture(6),Variance())))
}

# Get Copy number
GetCopyNumber <- function(Value){
  return(round(2*(2**Value),0))
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
}

LinearExpression(GainNarrowMediumAberration1(100))


#plot(c(NormalCell(1000),GainNarrowMediumAberration1(5000),NormalCell(1000)))
