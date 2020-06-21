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
  return(TempVector)
}

NormalExpression <- function(NumberOfProbes){
  return(c(rnorm(NumberOfProbes,6,sqrt(2))))
}


GenerateComplementaryGeneExpression <- function(Vector,Mode){
  TempVector = c()
  #if(Mode =="Normal"){
  #
  #  return(TempVector)
  #}


  for(i in Vector){

    if(Mode == "Linear"){

        TempVector = c(TempVector,LinearExpression(c(i)))
  
    }else if(Mode == "Sigmoid"){


        TempVector = c(TempVector,SigmoidExpression(c(i)))

    }else if (Mode=="StepWise"){
        TempVector = c(TempVector,StepWiseExpression(c(i)))
    }
  }

  return(TempVector)

}


CopyErrorCheck <- function(ProbeVector){
  TempVector <- c()
  for(i in 1:length(ProbeVector)){
    if(round(GetCopyNumber(ProbeVector[i])) != 2){
      TempVector <- c(TempVector,c(i))
    }

  }

  return(TempVector)
} 

LossList <- function(ProbeVector){

  TempVector <- c()
  for(i in 1:length(ProbeVector)){
    
    if(round(GetCopyNumber(ProbeVector[i])) < 2){
      TempVector <- c(TempVector,c(i))
    }

  }

  return(TempVector)



}

AmplifictionList <- function(ProbeVector){

  TempVector <- c()
  for(i in 1:length(ProbeVector)){

    if(round(GetCopyNumber(ProbeVector[i])) > 2){
      TempVector <- c(TempVector,c(i))
    }

  }

  return(TempVector)


}


# Ask for user input (which gene should get the type I,II,etc.)
# Hold x values for each Gene Position
GenePositionsTypeI = c() 
GenePositionsTypeII = c()
GenePositionsTypeIII = c()
GenePositionsTypeIV = c()
GenePositionsTypeV = c()

NumberOfTypeI = 3
NumberOfTypeII = 3
NumberOfTypeIII = 3
NumberOfTypeIV = 3
NumberOfTypeV = 3



CopyNumberExpression = c(NormalCell(100),GainNarrowMediumAberration1(500),NormalCell(100))

# To help decide where to put typeI-V genes If User does not specifiy anything in terms of location
ProbeRegionsWithCopyNumberProblems = CopyErrorCheck(CopyNumberExpression)  # For Type I - III
ProbeRegionsWithCopyNumberLoss = LossList(CopyNumberExpression) # For Type IV
ProbeRegionsWithCopyNumberGain = AmplifictionList(CopyNumberExpression) # For Type V



NumberOfPatients = 10
ModelType = "Sigmoid"
TotalNumberOfProbes = length(CopyNumberExpression)

# Setting probabilites of genes types:
#TypeI appears in all patients
#TypeII appears in 75%
NumberofTypeIIs = round(0.75 * NumberOfPatients)
GraphsThatWillHaveTypeIIs = c()
for(i in 1:NumberofTypeIIs){
  GraphsThatWillHaveTypeIIs = c(GraphsThatWillHaveTypeIIs,1)
}
for(i in 1:(NumberOfPatients-NumberofTypeIIs)){
  GraphsThatWillHaveTypeIIs = c(GraphsThatWillHaveTypeIIs,0)
}
GraphsThatWillHaveTypeIIs = sample(GraphsThatWillHaveTypeIIs)

#TypeIII appears in 50%
NumberofTypeIIIs = round(0.50 * NumberOfPatients)
GraphsThatWillHaveTypeIIIs = c()
for(i in 1:NumberofTypeIIIs){
  GraphsThatWillHaveTypeIIIs = c(GraphsThatWillHaveTypeIIIs,1)
}
for(i in 1:(NumberOfPatients-NumberofTypeIIIs)){
  GraphsThatWillHaveTypeIIIs = c(GraphsThatWillHaveTypeIIIs,0)
}
GraphsThatWillHaveTypeIIIs = sample(GraphsThatWillHaveTypeIIIs)



# Still have to add user input version! & Fix Indexes
# https://thomasleeper.com/Rcourse/Tutorials/plotcolors.html
for(k in 0:NumberOfPatients){

      GeneExpression = GenerateComplementaryGeneExpression(CopyNumberExpression,ModelType)
      # Color Change
      IndexCounter = 1
      # InsertTypeI
      IndexestoChangeColorOfTypeI = c()
      for(i in 0:NumberOfTypeI){
        IndexestoChangeColorOfTypeI <- c(IndexestoChangeColorOfTypeI ,c(ProbeRegionsWithCopyNumberProblems[IndexCounter])) 
        IndexCounter =  IndexCounter + 1
      }

      IndexestoChangeColorOfTypeII = c()
      if(GraphsThatWillHaveTypeIIs[k+1] == 1){
        for(i in 0:NumberOfTypeII){
        IndexestoChangeColorOfTypeII <- c(IndexestoChangeColorOfTypeII ,c(ProbeRegionsWithCopyNumberProblems[IndexCounter])) 
        IndexCounter =  IndexCounter + 1
      }
      }
      # https://thepracticalr.wordpress.com/2016/08/30/2-y-axis-plotting/
      # https://rpubs.com/riazakhan94/297778
      # https://stackoverflow.com/questions/38247907/how-to-set-the-y-range-in-boxplot-graph


  
      xcor <- 1:TotalNumberOfProbes
      print("Generating Graphs")
      plot(xcor,CopyNumberExpression,col = "blue",  ylim=c(-2,3))
      par(new = TRUE)
      plot(xcor,GeneExpression ,xaxt = "n", yaxt = "n",
          ylab = "", xlab = "", col=ifelse(xcor ==IndexestoChangeColorOfTypeI,"red","black"), lty = 2,ylim=c(-8,12),
          bg=ifelse(xcor==IndexestoChangeColorOfTypeI,"red","white"))
      axis(side = 4)
      print("Graphs Generated.")
      print("=======")
}
