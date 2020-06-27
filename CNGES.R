# Optional Seed (uncomment the next line and add integer in place of x):
# set.seed(x)

# args hold the command line arguments neccessary for the program
args = commandArgs()
# Model connects copy number to gene expression{Sigmoid,Linear,Stepwise}
Model <- "Linear" #args[3]
# The number of graphs/patients we will generate
NumberOfPatients <- 10 #args[4]
# File that holds information about how many probes there are and what their regions will be(gain,loss,normal and variations - for each if they exist)
ProbeLocationFile <- "./ProbeLocation.txt"#args[5]


# Setup functions and values used in CGH value assignment:
# Formulae:
# Cell admixture is --> M(c) = log2 [(c Â· Pt + 2 Â· (1 â Pt))/2], where Pt = U(0.3, 0.7)
# N(Âµ = M(c), Ï^{2} = S^{2}), where S[is represented by the Variance function below] = U(min = 0.1, max = 0.3)
CellAdmixture <- function(value){
  Pt <- runif(1,0.3,0.7)
  Result <- log2( (value * Pt + 2 *(1-Pt))/2 )
  return(Result)
}
Variance <- function(value){
  return(runif(1,0.1,0.3)**2)
}

# Open the probe file and read the lines from it
ProbeLocation <- readLines(paste(getwd(), ProbeLocationFile, sep = ""))

# Hashtable that holds the different inputs and their corresponding m values
InputTranslationTable <- new.env(hash=T, parent=emptyenv())


# Gain Narrow Medium Aberration1
InputTranslationTable[["GNM1"]] <- 8
# Gain Narrow Medium Aberration2
InputTranslationTable[["GNM2"]] <- 6
# Loss Narrow Medium Aberration
InputTranslationTable[["LNMA"]] <- 1
# Gain Large High Aberration
InputTranslationTable[["GLHA"]] <- 10
# Loss Large High Aberration
InputTranslationTable[["LLHA"]] <- 0
# Gain Broad Low Abberation 
InputTranslationTable[["GBLA"]] <- 3
# Normal Cell
InputTranslationTable[["NOCE"]] <- 2



# Get Copy number
GetCopyNumber <- function(Value){
  return(2*(2**Value))
}


# Gene Expression Simulation(Linear)
LinearExpression <- function(Vector){
  TempVector <- c()
  for (i in Vector){
    # Forumla
    CurrCopyNum <- GetCopyNumber(i)
    CurrValue <- 2*CurrCopyNum+2
    #Abberation
    if (CurrCopyNum >= 3 | CurrCopyNum == 0 | CurrCopyNum == 1){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else{
      # Normal Expression
      TempVector <- c(TempVector,rnorm(1,CurrValue,sqrt(2)))
    }
  }
  return(TempVector)
}

# Gene Expression Simulation(Sigmoid)
SigmoidExpression <- function(Vector){
  TempVector <- c()
  for (i in Vector){
    # Formula
    CurrCopyNum <- GetCopyNumber(i)
    el <- (CurrCopyNum/2) + 2
    CurrValue <- 8/(1+exp(-1*el))
    # Abberation
    if (CurrCopyNum >= 3 | CurrCopyNum == 0 | CurrCopyNum == 1){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else{
      # Normal Expression
      TempVector <- c(TempVector,rnorm(1,CurrValue,sqrt(2)))
    }
  }
  return(TempVector)
}

# Gene Expression Simulation(Step Wise)
StepWiseExpression <- function(Vector){
  TempVector <- c()
  for (i in Vector){
    CurrCopyNum <- GetCopyNumber(i)
    # Formula
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
    
    # Abberation
    if (CurrCopyNum >= 3 | CurrCopyNum == 0 | CurrCopyNum == 1){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else{
      # Normal Expression
      TempVector <- c(TempVector,rnorm(1,CurrValue,sqrt(2)))
    }
  }
  return(TempVector)
}


# Generate corresponding gene expression given a vector with cgh values(each element in which represents a probe)
GenerateComplementaryGeneExpression <- function(Vector,Model){
  TempVector <- c()

  if(Model == "Linear"){
      
    TempVector <- c(TempVector,LinearExpression(Vector))
      
  }else if(Model == "Sigmoid"){
     
    TempVector <- c(TempVector,SigmoidExpression(Vector))
      
  }else if (Model =="StepWise"){
    TempVector <- c(TempVector,StepWiseExpression(Vector))
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






FillTypeIVProbeLocationVector <- function(ProbeRegionsWithCopyNumberLoss,MasterListOfAllGeneTypePositions){
NumberOfTypeIV = 3
Flag <- 3
Index <- 1
InsertTheseXvalues <- c()
# Go through Vector that holds all the probes with copy number problems
# When you find a free spot(i.e one that a user has not specifed and is in the problem region)
# Insert that x value in the file
while(Flag != 0 & Index <= length(CopyNumberExpression)){
  if(ProbeRegionsWithCopyNumberLoss[Index] %in% MasterListOfAllGeneTypePositions){
    Index <- Index + 1
    next
  }else{
    InsertTheseXvalues <- c(InsertTheseXvalues,c(ProbeRegionsWithCopyNumberLoss[Index]))
    MasterListOfAllGeneTypePositions <- c(MasterListOfAllGeneTypePositions, c(ProbeRegionsWithCopyNumberLoss[Index]))
    Flag <- Flag - 1
  }
  
}
return(InsertTheseXvalues)
}

FillTypeVProbeLocationVector <- function(ProbeRegionsWithCopyNumberGain,MasterListOfAllGeneTypePositions){
  NumberOfTypeIV = 3
  Flag <- 3
  Index <- 1
  InsertTheseXvalues <- c()
  # Go through Vector that holds all the probes with copy number problems
  # When you find a free spot(i.e one that a user has not specifed and is in the problem region)
  # Insert that x value in the file
  while(Flag != 0 & Index <= length(CopyNumberExpression)){
    if(ProbeRegionsWithCopyNumberGain[Index] %in% MasterListOfAllGeneTypePositions){
      Index <- Index + 1
      next
    }else{
      InsertTheseXvalues <- c(InsertTheseXvalues,c(ProbeRegionsWithCopyNumberGain[Index]))
      MasterListOfAllGeneTypePositions <- c(MasterListOfAllGeneTypePositions, c(ProbeRegionsWithCopyNumberGain[Index]))
      Flag <- Flag - 1
    }
    
  }
  return(InsertTheseXvalues)
}






# Setting probabilites of genes types:
#TypeI appears in all patients
#TypeII appears in 75%
NumberofTypeIIs <- round(0.75 * NumberOfPatients)
GraphsThatWillHaveTypeIIs <- c()
for(i in 1:NumberofTypeIIs){
  GraphsThatWillHaveTypeIIs <- c(GraphsThatWillHaveTypeIIs,1)
}
for(i in 1:(NumberOfPatients-NumberofTypeIIs)){
  GraphsThatWillHaveTypeIIs <- c(GraphsThatWillHaveTypeIIs,0)
}
GraphsThatWillHaveTypeIIs <- sample(GraphsThatWillHaveTypeIIs)

#TypeIII appears in 50%
NumberofTypeIIIs <- round(0.50 * NumberOfPatients)
GraphsThatWillHaveTypeIIIs <- c()
for(i in 1:NumberofTypeIIIs){
  GraphsThatWillHaveTypeIIIs <- c(GraphsThatWillHaveTypeIIIs,1)
}
for(i in 1:(NumberOfPatients-NumberofTypeIIIs)){
  GraphsThatWillHaveTypeIIIs <- c(GraphsThatWillHaveTypeIIIs,0)
}
GraphsThatWillHaveTypeIIIs <- sample(GraphsThatWillHaveTypeIIIs)











for(k in 1:NumberOfPatients){
  # Get copy number simulation by reading lines and giving each probe its appropriate value
  CopyNumberExpression <- c()
  for(Line in ProbeLocation){
    SplitLine <- unlist(strsplit(Line, split = "\t"))
    CopyNumberExpression <- c(CopyNumberExpression,c(rnorm(as.integer(SplitLine[1]),CellAdmixture(InputTranslationTable[[SplitLine[2]]]),Variance())))
  }  
  
  # Generate Gene Expression
  GeneExpression <- GenerateComplementaryGeneExpression(CopyNumberExpression,Model)
  
  
  # To help decide where to put typeI-V genes If User does not specifiy anything in terms of location
  ProbeRegionsWithCopyNumberProblems <- CopyErrorCheck(CopyNumberExpression)  # For Type I - III
  ProbeRegionsWithCopyNumberLoss <- LossList(CopyNumberExpression) # For Type IV
  ProbeRegionsWithCopyNumberGain <- AmplifictionList(CopyNumberExpression) # For Type V
  
  
  ErrorProbes = CopyErrorCheck(CopyNumberExpression)
  TIProbeLocation <- c(ErrorProbes[1], ErrorProbes[2], ErrorProbes[3])
  TIIProbeLocation <- c(ErrorProbes[4], ErrorProbes[5], ErrorProbes[6])
  TIIIProbeLocation <- c(ErrorProbes[7], ErrorProbes[8], ErrorProbes[9])
  TakenProbes <- c(TIProbeLocation,TIIProbeLocation,TIIIProbeLocation)
  
  if(length(ProbeRegionsWithCopyNumberLoss) != 0){
    TIVProbeLocation <- FillTypeIVProbeLocationVector(ProbeRegionsWithCopyNumberLoss,TakenProbes)
    TakenProbes <- c(TakenProbes,TIVProbeLocation)
  }else{
    TIVProbeLocation <- c()
  }
  
  TVProbeLocation <- FillTypeVProbeLocationVector(ProbeRegionsWithCopyNumberGain,TakenProbes)
  TakenProbes <- c(TakenProbes,TVProbeLocation)
  
  SizeAdj <- c()
  for(i in 1:length(CopyNumberExpression)){
    if(i %in% TakenProbes){
      SizeAdj <- c(SizeAdj,0.6)
    }else{
      SizeAdj <- c(SizeAdj,0.3)
    }
  }
  
  
  
GeneEpxrCol <- c()
BackgroundCol <- c()
for(i in 1:length(CopyNumberExpression)){
  if(i %in% TIProbeLocation){
    GeneEpxrCol <- c(GeneEpxrCol,"red")
    BackgroundCol <- c(BackgroundCol,"red")
  }else if(i %in% TIIProbeLocation && GraphsThatWillHaveTypeIIs[k] == 1){
    GeneEpxrCol <- c(GeneEpxrCol,"yellow")
    BackgroundCol <- c(BackgroundCol,"yellow")
  }else if(i %in% TIIIProbeLocation && GraphsThatWillHaveTypeIIIs[k] == 1){
    GeneEpxrCol <- c(GeneEpxrCol,"orange")
    BackgroundCol <- c(BackgroundCol,"orange")
  }else if(i %in% TIVProbeLocation){
    GeneExpression[i] <-rnorm(1,4,sqrt(0.5))
    GeneEpxrCol <- c(GeneEpxrCol,"purple")
    BackgroundCol <- c(BackgroundCol,"purple")
  }else if(i %in% TVProbeLocation){
    GeneExpression[i] <-rnorm(1,12,sqrt(0.5))
    GeneEpxrCol <- c(GeneEpxrCol,"blue")
    BackgroundCol <- c(BackgroundCol,"blue")
  }else{
    GeneEpxrCol <- c(GeneEpxrCol,"gray")
    BackgroundCol <-c(BackgroundCol,"gray")
  }
}
#plot(xcor,GeneExpression, cex=0.3 ,xaxt = "n", yaxt = "n",
#ylab = "", xlab = "", col=ifelse(xcor %in% TIProbeLocation,"red",ifelse(xcor %in% TIIProbeLocation,"yellow",ifelse(xcor %in% TIIIProbeLocation,"orange",ifelse(xcor %in% TIVProbeLocation,"purple",ifelse(xcor %in% TVProbeLocation,"darkgreen","gray"))))), lty = 2,ylim=c(-8,12))


xcor <- 1:length(CopyNumberExpression)
print("Generating Graphs")
plot(xcor,CopyNumberExpression,col = "blue",  ylim=c(-2,3), cex=0.3)
par(new = TRUE)
plot(xcor,GeneExpression,col=GeneEpxrCol,bg=BackgroundCol,pch=21 ,cex=SizeAdj ,xaxt = "n", yaxt = "n",
     ylab = "", xlab = "",lty = 2,ylim=c(-8,12))

axis(side = 4)
print("Graphs Generated.")
print("=======")
}
