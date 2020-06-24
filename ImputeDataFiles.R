# Setup functions and values:
CellAdmixture <- function(value){
  Pt <- runif(1,0.3,0.7)
  Result <- log2( (value * Pt + 2 *(1-Pt))/2 )
  return(Result)
}
Variance <- function(value){
  return(runif(1,0.1,0.3)**2)
}

NormalCell <- function(NumberOfProbes){
  return(c(rnorm(NumberOfProbes,CellAdmixture(2),Variance())))
}

GainNarrowMediumAberration1 <- function(NumberOfProbes){
  return(c(rnorm(NumberOfProbes,CellAdmixture(6),Variance())))
}

args = commandArgs()

GeneTypeFile <- "./GeneType.txt"#args[3]
ProbeLocationFile <- "./ProbeLocation.txt"#args[4]

ProbeLocation <- readLines(paste(getwd(), ProbeLocationFile, sep = ""))
GeneType <- readLines(paste(getwd(),GeneTypeFile, sep = ""))



CopyNumberExpression <- c()
if(length(ProbeLocation) == 0){
  CopyNumberExpression <- c(NormalCell(10),GainNarrowMediumAberration1(50),NormalCell(10))
}else{
  for(i in ProbeLocation){
    SplitString <- unlist(strsplit(i, split = "\t"))
    CopyNumberExpression <- c(CopyNumberExpression,c(rnorm(as.integer(SplitString[1]),CellAdmixture(as.integer(SplitString[2])),Variance())))
  }
}



GetCopyNumber <- function(Value){
  #return(round(2*(2**Value),0))
  return(2*(2**Value))
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




# To help decide where to put typeI-V genes If User does not specifiy anything in terms of location
ProbeRegionsWithCopyNumberProblems = CopyErrorCheck(CopyNumberExpression)  # For Type I - III
ProbeRegionsWithCopyNumberLoss = LossList(CopyNumberExpression) # For Type IV
ProbeRegionsWithCopyNumberGain = AmplifictionList(CopyNumberExpression) # For Type V




SplitString1 <- as.integer(unlist(strsplit(GeneType[1], split = "\t")))
SplitString2 <- as.integer(unlist(strsplit(GeneType[2], split = "\t")))
SplitString3 <- as.integer(unlist(strsplit(GeneType[3], split = "\t")))
SplitString4 <- as.integer(unlist(strsplit(GeneType[4], split = "\t")))
SplitString5 <- as.integer(unlist(strsplit(GeneType[5], split = "\t")))



NumberOfTypeII <- length(SplitString2)
NumberOfTypeIII <- length(SplitString3)
NumberOfTypeIV <- length(SplitString4)
NumberOfTypeV <- length(SplitString5)


MasterListOfAllGeneTypePositions <- c(SplitString1,SplitString2,SplitString3,SplitString4,SplitString5)





if(length(SplitString1) == 0 || is.na(SplitString1) ){
  NumberOfTypeI <- 3
  Flag <- 3
  Index <- 1
  InsertTheseXvalues <- c()
  # Go through Vector that holds all the probes with copy number problems
  # When you find a free spot(i.e one that a user has not specifed and is in the problem region)
  # Insert that x value in the file
  while(Flag != 0 & Index <= length(CopyNumberExpression)){
    if(ProbeRegionsWithCopyNumberProblems[Index] %in% MasterListOfAllGeneTypePositions){
      Index <- Index + 1
      next
    }else{
      InsertTheseXvalues <- c(InsertTheseXvalues,c(ProbeRegionsWithCopyNumberProblems[Index]))
      MasterListOfAllGeneTypePositions <- c(MasterListOfAllGeneTypePositions, c(ProbeRegionsWithCopyNumberProblems[Index]))
      Flag <- Flag - 1
    }
    
  }
  
  NewString <- paste(InsertTheseXvalues, collapse = "\t")
  sink(GeneTypeFile)
  cat(NewString)
  cat("\n")
  cat(paste(SplitString2, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString3, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString4, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString5, collapse = "\t"))
  cat("\n")
  sink()

  
}

if(length(SplitString2) == 0 || is.na(SplitString2)){
  NumberOfTypeII = 3
  Flag <- 3
  Index <- 1
  InsertTheseXvalues <- c()
  # Go through Vector that holds all the probes with copy number problems
  # When you find a free spot(i.e one that a user has not specifed and is in the problem region)
  # Insert that x value in the file
  while(Flag != 0 & Index <= length(CopyNumberExpression)){
    if(ProbeRegionsWithCopyNumberProblems[Index] %in% MasterListOfAllGeneTypePositions){
      Index <- Index + 1
      next
    }else{
      InsertTheseXvalues <- c(InsertTheseXvalues,c(ProbeRegionsWithCopyNumberProblems[Index]))
      MasterListOfAllGeneTypePositions <- c(MasterListOfAllGeneTypePositions, c(ProbeRegionsWithCopyNumberProblems[Index]))
      Flag <- Flag - 1
    }
    
  }
  
  NewString <- paste(InsertTheseXvalues, collapse = "\t")
  sink(GeneTypeFile)
  cat(paste(SplitString1, collapse = "\t"))
  cat("\n")
  cat(NewString)
  cat("\n")
  cat(paste(SplitString3, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString4, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString5, collapse = "\t"))
  cat("\n")
  sink()
  
  
}

if(length(SplitString3) == 0 || is.na(SplitString3)){
  NumberOfTypeIII = 3
  Flag <- 3
  Index <- 1
  InsertTheseXvalues <- c()
  # Go through Vector that holds all the probes with copy number problems
  # When you find a free spot(i.e one that a user has not specifed and is in the problem region)
  # Insert that x value in the file
  while(Flag != 0 & Index <= length(CopyNumberExpression)){
    if(ProbeRegionsWithCopyNumberProblems[Index] %in% MasterListOfAllGeneTypePositions){
      Index <- Index + 1
      next
    }else{
      InsertTheseXvalues <- c(InsertTheseXvalues,c(ProbeRegionsWithCopyNumberProblems[Index]))
      MasterListOfAllGeneTypePositions <- c(MasterListOfAllGeneTypePositions, c(ProbeRegionsWithCopyNumberProblems[Index]))
      Flag <- Flag - 1
    }
    
  }
  
  NewString <- paste(InsertTheseXvalues, collapse = "\t")
  sink(GeneTypeFile)
  cat(paste(SplitString1, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString2, collapse = "\t"))
  cat("\n")
  cat(NewString)
  cat("\n")
  cat(paste(SplitString4, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString5, collapse = "\t"))
  cat("\n")
  sink()
  
}

if((length(SplitString4) == 0 || is.na(SplitString4)) & length(ProbeRegionsWithCopyNumberLoss) != 0 ){
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
  
  NewString <- paste(InsertTheseXvalues, collapse = "\t")
  sink(GeneTypeFile)
  cat(paste(SplitString1, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString2, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString3, collapse = "\t"))
  cat("\n")
  cat(NewString)
  cat("\n")
  cat(paste(SplitString5, collapse = "\t"))
  cat("\n")
  sink()
  
}


if((length(SplitString5) == 0 || is.na(SplitString5)) & length(ProbeRegionsWithCopyNumberGain) != 0 ){
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
  
  NewString <- paste(InsertTheseXvalues, collapse = "\t")
  sink(GeneTypeFile)
  cat(paste(SplitString1, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString2, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString3, collapse = "\t"))
  cat("\n")
  cat(paste(SplitString4, collapse = "\t"))
  cat("\n")
  cat(NewString)
  cat("\n")
  sink()
  
}
