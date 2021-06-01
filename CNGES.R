# Optional Seed (uncomment the next line and add integer in place of x):
# set.seed(1309)
options(scipen = 999)

# args hold the command line arguments neccessary for the program
args = commandArgs()
# Model connects copy number to gene expression{Sigmoid,Linear,Stepwise}
Model <- args[6] # "StepWise"
if(Model != "StepWise" && Model != "Sigmoid" && Model != "Linear"){
  print("Please enter correct Model [StepWise,Sigmoid,Linear]")
  quit()
}
# The number of graphs/patients we will generate
NumberOfPatients <- as.integer(args[7]) # 3
# File that holds information about how many probes there are and what their regions will be(gain,loss,normal and variations - for each if they exist)
ProbeLocationFile <- args[8] #"ProbeLdocation.txt"
# Give average data
GiveAverageData <- as.logical(args[9]) #as.logical("FALSE")
# Setting for Edira (Auto False)
EdiraDataSet <- as.logical(args[10]) # as.logical("FALSE")
EdiraOffSet <- as.integer(args[11])  # 0
# CSVFormat if not "null" then output in this file
CSVFormat <- args[12]  # "null"

# File to input raw values into
CGHOutputFile <- args[13] #"C:/Users/ismai/Desktop/CopyNumAndGeneExpSimulation/CGH.txt"
GeneExpOutputFile<- args[14] #"C:/Users/ismai/Desktop/CopyNumAndGeneExpSimulation/GeneExp.txt"

GeneExpOutputFileCSV <- args[15] #"C:/Users/ismai/Desktop/CopyNumAndGeneExpSimulation/GeneExp.csv"
CGHOutputFileCSV <- args[16] #"C:/Users/ismai/Desktop/CopyNumAndGeneExpSimulation/CGH.csv"

FalseTruePositivesProbesLocations <- args[17] #"C:/Users/ismai/Desktop/CopyNumAndGeneExpSimulation/Postives.txt"

# Check if a file exists, if not generate one with a default name and notify user
if(file.exists(ProbeLocationFile)){
  # Notify User that the file was successfully read in
  print("File Detected Successfully!")
}else{
  # Intialize Boilerplate text
  BoilerPlateText <- c("ADD LINES FOR EACH SUBSEQUENT REGION OF PROBES BELOW",
                       "EXAMPLE: \"100 NOCE\" is a region with 100 probes with normal copy number expression.",
                       "\n","\n",
                       "================================================================================",
                       "ENTER THE POSITIONS OF EACH OF THE GENE TYPES BELOW",
                       "EXAMPLE: \"Type 1: 1,100,13\" Type 1 genes are places at x positions 1,100 and 13.",
                       "Type 1: ",
                       "Type 2: ",
                       "Type 3: ",
                       "Type 4: ",
                       "Type 5: ")
  # Create file
  file.create(ProbeLocationFile)
  # Add in boiler plate text
  writeLines(BoilerPlateText,ProbeLocationFile)
  # Notify User
  print("Please fill out the generated file!")
  # End Program
  quit()
}

# Setup functions and values used in CGH value assignment:
# Formulae:
# Cell admixture is --> M(c) = log2 [(c * Pt + 2 * (1 ??? Pt))/2], where Pt = U(0.3, 0.7)
# N(mean = M(c), std^{2} = S^{2}), where S[is represented by the Variance function below] = U(min = 0.1, max = 0.3)
CellAdmixture <- function(value){
  Pt <- runif(1,0.3,0.7)
  Result <- log2( (value * Pt + 2 *(1-Pt))/2 )
  return(Result)
}
CellAdmixtureForEdira <- function(value,offset){
  Pt <- runif(1,0.3,0.7)
  Result <- ( (value * Pt + 2 *(1-Pt))/2 ) + offset
  return(Result)
}
STD <- function(value){
  return(runif(1,0.1,0.3))
}

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

# Generate corresponding gene expression given a vector with cgh values(each element in which represents a probe)
GenerateComplementaryGeneExpression <- function(Vector,TakenProbes,Model){
  TempVector <- c()
  if(Model == "Linear"){
    
    TempVector <- c(TempVector,LinearExpression(Vector,TakenProbes))
    
  }else if(Model == "Sigmoid"){
    
    TempVector <- c(TempVector,SigmoidExpression(Vector, TakenProbes))
    
  }else if (Model =="StepWise"){
    TempVector <- c(TempVector,StepWiseExpression(Vector, TakenProbes))
  }else{
    print("MODEL FATAL FAILURE -> No suitable model found")
    quit()
  }
  
  return(TempVector)
  
}


# Gene Expression Simulation(Linear)
LinearExpression <- function(Vector,TakenProbes){
  TempVector <- c()
  # They are generated in sequential order [from type I to V]
  for (i in Vector){
    # Forumla
    CurrCopyNum <- GetCopyNumber(i)
    CurrValue <- 2*CurrCopyNum+2
    # When index is belongs to a gene type
    if (i %in% TakenProbes$T1Loc){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else if(i %in% TakenProbes$T2Loc && runif(1)<0.75){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else if(i %in% TakenProbes$T3Loc && runif(1)<0.50){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else if(i %in% TakenProbes$T4Loc){
      TempVector <- c(TempVector,rnorm(1,4,sqrt(0.5)))
    }else if(i %in% TakenProbes$T5Loc){
      TempVector <- c(TempVector,rnorm(1,12,0.5))
    }else{
      # Normal Expression
      TempVector <- c(TempVector,rnorm(1,6,sqrt(2)))
    }
  }
  return(TempVector)
}



# Gene Expression Simulation(Linear)
SigmoidExpression <- function(Vector,TakenProbes){
  TempVector <- c()
  # They are generated in sequential order [from type I to V]
  for (i in Vector){
    # Forumla
    CurrCopyNum <- GetCopyNumber(i)
    el <- (CurrCopyNum/2) + 2
    CurrValue <- 8/(1+exp(-1*el))
    # When index is belongs to a gene type
    if (i %in% TakenProbes$T1Loc){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else if(i %in% TakenProbes$T2Loc && runif(1)<0.75){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else if(i %in% TakenProbes$T3Loc && runif(1)<0.50){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else if(i %in% TakenProbes$T4Loc){
      TempVector <- c(TempVector,rnorm(1,4,sqrt(0.5)))
    }else if(i %in% TakenProbes$T5Loc){
      TempVector <- c(TempVector,rnorm(1,12,0.5))
    }else{
      # Normal Expression
      TempVector <- c(TempVector,rnorm(1,6,sqrt(2)))
    }
  }
  return(TempVector)
}


# Gene Expression Simulation(Linear)
StepWiseExpression <- function(Vector,TakenProbes){
  TempVector <- c()
  # They are generated in sequential order [from type I to V]
  for (i in Vector){
    # Forumla
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
    # When index is belongs to a gene type
    if (i %in% TakenProbes$T1Loc){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else if(i %in% TakenProbes$T2Loc && runif(1)<0.75){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else if(i %in% TakenProbes$T3Loc && runif(1)<0.50){
      TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
    }else if(i %in% TakenProbes$T4Loc){
      TempVector <- c(TempVector,rnorm(1,4,0.5))
    }else if(i %in% TakenProbes$T5Loc){
      TempVector <- c(TempVector,rnorm(1,12,0.5))
    }else{
      # Normal Expression
      TempVector <- c(TempVector,rnorm(1,6,sqrt(2)))
    }
  }
  return(TempVector)
}




# Open the probe file and read the lines from it
RawFileData <- readLines(ProbeLocationFile)
GeneLocation <- tail(RawFileData,n=5)
ProbeLocation <- RawFileData[3:(length(RawFileData)-8)]

# Convert file data(about gene type locations) into usable vectors
Type1Location <- unlist(strsplit(substr(GeneLocation[1],9,nchar(GeneLocation[1])), split = "\t"))
if(length(Type1Location) == 0){
  Type1Location <- -1
}

Type2Location <- unlist(strsplit(substr(GeneLocation[2],9,nchar(GeneLocation[2])), split = "\t"))
if(length(Type2Location) == 0){
  Type2Location <- -1
}

Type3Location <- unlist(strsplit(substr(GeneLocation[3],9,nchar(GeneLocation[3])), split = "\t"))
if(length(Type3Location) == 0){
  Type3Location <- -1
}

Type4Location <- unlist(strsplit(substr(GeneLocation[4],9,nchar(GeneLocation[4])), split = "\t"))
if(length(Type4Location) == 0){
  Type4Location <- -1
}

Type5Location <- unlist(strsplit(substr(GeneLocation[5],9,nchar(GeneLocation[5])), split = "\t"))
if(length(Type5Location) == 0){
  Type5Location <- -1
}

TakenProbes <- data.frame("T1Loc" = Type1Location,"T2Loc" = Type2Location,"T3Loc" = Type3Location,"T4Loc" = Type4Location,"T5Loc" = Type5Location, stringsAsFactors=F)



ProbeSetForAverage <- c()
GeneSetForAverage <- c()
ProbeSetForFull <- c()
GeneSetForFull <- c()
# Generate a graph for each simulated patient
PatientCGHSTD <-STD()
for(k in 1:NumberOfPatients){
  # Get copy number simulation by reading lines and giving each probe its appropriate value
  CopyNumberExpression <- c()
  for(Line in ProbeLocation){
    # Split tab-delimited line into two(the first value is the number of probes with a specific CGH value
    # and the second value is the actual CGH ratio value assigned to those probes)
    SplitLine <- unlist(strsplit(Line, split = "\t"))
    # add the appropriate amount and type of genes
    if(!EdiraDataSet){
      CopyNumberExpression <- c(CopyNumberExpression,c(rnorm(as.integer(SplitLine[1]),CellAdmixture(InputTranslationTable[[SplitLine[2]]]),PatientCGHSTD)))
    }else{
      CopyNumberExpression <- c(CopyNumberExpression,c(rnorm(as.integer(SplitLine[1]),CellAdmixtureForEdira(InputTranslationTable[[SplitLine[2]]],EdiraOffSet),PatientCGHSTD)))
    }
  }  
  
  if(length(ProbeSetForAverage) == 0 ){
    ProbeSetForAverage <- numeric(length(CopyNumberExpression))
    GeneSetForAverage <- numeric(length(CopyNumberExpression))
  }
  # Insert Gene expression values (base,without gene types)
  GeneExpression <- GenerateComplementaryGeneExpression(CopyNumberExpression,TakenProbes,Model)
  
  SizeAdj <- c()
  # Size adjustment for normal and typeI-V genes
  for(i in 1:length(CopyNumberExpression)){
    # if a probe is a typeI-V gene then make it double the size of a normal gene expression dot
    if(any(TakenProbes == i)){
      SizeAdj <- c(SizeAdj,0.6)
    }else{
      SizeAdj <- c(SizeAdj,0.3)
    }
  }
  
  
  GeneEpxrCol <- c()
  BackgroundCol <- c()
  # Loop sets the color of typeI-V genes differently(in order to make them easier to view) 
  # & changes the values of the type IV and V as mandated by the source paper
  for(i in 1:length(CopyNumberExpression)){
    # If a probe is type 1 then make it a red dot
    if(i %in% TakenProbes$T1Loc){
      GeneEpxrCol <- c(GeneEpxrCol,"red")
      BackgroundCol <- c(BackgroundCol,"red")
    }else if(i %in% TakenProbes$T2Loc ){
      # If a probe is type 2 then make it a yellow dot
      GeneEpxrCol <- c(GeneEpxrCol,"yellow")
      BackgroundCol <- c(BackgroundCol,"yellow")
    }else if(i %in% TakenProbes$T3Loc){
      # If a probe is type 1 then make it a orange dot
      GeneEpxrCol <- c(GeneEpxrCol,"orange")
      BackgroundCol <- c(BackgroundCol,"orange")
    }else if(i %in% TakenProbes$T4Loc){
      # If a probe is type 4 then make it a purple dot and change it's value(according to the source paper)
      GeneEpxrCol <- c(GeneEpxrCol,"purple")
      BackgroundCol <- c(BackgroundCol,"purple")
    }else if(i %in% TakenProbes$T5Loc){
      # If a probe is type 5 then make it a purple dot and change it's value(according to the source paper)
      GeneEpxrCol <- c(GeneEpxrCol,"blue")
      BackgroundCol <- c(BackgroundCol,"blue")
    }else{
      # If a probe is a regular gene expression value then make it gray
      GeneEpxrCol <- c(GeneEpxrCol,"gray")
      BackgroundCol <-c(BackgroundCol,"gray")
    }
  }
  
  if(GiveAverageData){
    for(i in 1:length(CopyNumberExpression)){
      ProbeSetForAverage[i] <- ProbeSetForAverage[i] + CopyNumberExpression[i]
      GeneSetForAverage[i] <- GeneSetForAverage[i] + GeneExpression[i]
    }
  }else{
    ProbeSetForFull <- c(ProbeSetForFull,CopyNumberExpression)
    GeneSetForFull <- c(GeneSetForFull,GeneExpression)
  }
  # Create the x-coordinates for the graph
  xcor <- 1:length(CopyNumberExpression)
  print("Generating Graphs")
  # Plot the first part of the graph(cgh ratio/blue dots)
  if(!EdiraDataSet){
    plot(xcor,CopyNumberExpression,col = "blue",ylim=c(-2,3), cex=0.3 ,ylab="CGH Ratio Values", xlab="Probes/Genes")
  }else{
    plot(xcor,CopyNumberExpression,col = "blue",ylim=c(-1,3), cex=0.3 ,ylab="CGH Ratio Values", xlab="Probes/Genes")
  }
  # Give the graph a title that denotes the patient number
  
  title(paste("Patient ",as.character(k)))
  # Adjust the size of the plot and allow the second part(gene expression) to be added on to the first(CGH ratio values)
  par(new = TRUE,mar=c(5.1,4.1,4.1,4.6))
  # Plot the second part of the graph(gene expression/gray dots)
  plot(xcor,GeneExpression,col=GeneEpxrCol,bg=BackgroundCol,pch=21 ,cex=SizeAdj ,xaxt = "n", yaxt = "n",
       ylab = "", xlab = "",lty = 2,ylim=c(-8,12))
  # Set a y-axis for the gene expression portion of the graph
  axis(side = 4)
  # Add the gene expression y-axis text
  mtext("Gene Expression Values",side=4,line=2.5)
  
  print("Graphs Generated.")
  print("=======")
  
  # Add marker to seperate patients at end
  
}



# Output Data into these files
ProbeSetForAverage <- ProbeSetForAverage/length(ProbeSetForAverage)
GeneSetForAverage <- GeneSetForAverage/length(GeneSetForAverage)


if(GiveAverageData){
  writeLines(as.character(ProbeSetForAverage),CGHOutputFile)
  writeLines(as.character(GeneSetForAverage),GeneExpOutputFile)
  
}else if(CSVFormat != "null"){
  ResultLines <- c()
  for(i in 1:length(CopyNumberExpression)){
    ResultLines <- c(ResultLines,paste(as.character(CopyNumberExpression[i]),as.character(GeneExpression[i]),sep = ","))
  }
  writeLines(ResultLines,CSVFormat)
}else{
  writeLines(as.character(c(paste("Num. of Patients:",as.character(NumberOfPatients)),
                            paste("Probers Per Patient:",as.character(length(CopyNumberExpression))),
                            ProbeSetForFull)),CGHOutputFile)
  writeLines(as.character(c(paste("Num. of Patients:",as.character(NumberOfPatients)),
                            paste("Probers Per Patient:",as.character(length(CopyNumberExpression))),
                            GeneSetForFull)),GeneExpOutputFile)
  
  # Prep matrix for data entry of the gene data
  GenePrintMatrix <- matrix(ncol=NumberOfPatients,nrow=length(CopyNumberExpression))

  StartFrame <- 1
  for(i in 1:NumberOfPatients){
    
    GenePrintMatrix[,i] <- GeneSetForFull[StartFrame:(i*length(CopyNumberExpression))]
    StartFrame <- StartFrame + length(CopyNumberExpression)
    
  }
  # Write to CSV (Gene)
  write.table(x=GenePrintMatrix,file=GeneExpOutputFileCSV,sep=",",row.names = FALSE,col.names = FALSE)
  
  
  
  
  
  # Prep matrix for data entry of the CGH data
  CGHPrintMatrix <- matrix(ncol=NumberOfPatients,nrow=length(CopyNumberExpression))
  
  StartFrame <- 1
  for(i in 1:NumberOfPatients){
    
    CGHPrintMatrix[,i] <- ProbeSetForFull[StartFrame:(i*length(CopyNumberExpression))]
    StartFrame <- StartFrame + length(CopyNumberExpression)
    
  }
  # Write to CSV (GGH)
  write.table(x=CGHPrintMatrix,file=CGHOutputFileCSV,sep=",",row.names = FALSE,col.names = FALSE)
  

}
writeLines(as.character(c(
                          paste("T1:",length(TakenProbes$T1Loc)),
                          paste("T2:",length(TakenProbes$T2Loc)),
                          paste("T3:",length(TakenProbes$T3Loc)),
                          paste("T4:",length(TakenProbes$T4Loc)),
                          paste("T5:",length(TakenProbes$T5Loc)),
                          TakenProbes$T1Loc,
                          TakenProbes$T2Loc,
                          TakenProbes$T3Loc,
                          TakenProbes$T4Loc,
                          TakenProbes$T5Loc))
           ,FalseTruePositivesProbesLocations)

