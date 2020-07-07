# Optional Seed (uncomment the next line and add integer in place of x):
# set.seed(x)

# args hold the command line arguments neccessary for the program
args = commandArgs()
# Model connects copy number to gene expression{Sigmoid,Linear,Stepwise}
Model <- "Linear"#args[4]
# The number of graphs/patients we will generate
NumberOfPatients <- 10#as.integer(args[5])
# File that holds information about how many probes there are and what their regions will be(gain,loss,normal and variations - for each if they exist)
ProbeLocationFile <- "/ProbeLocation.txt"#args[6]


# Setup functions and values used in CGH value assignment:
# Formulae:
# Cell admixture is --> M(c) = log2 [(c * Pt + 2 * (1 ??? Pt))/2], where Pt = U(0.3, 0.7)
# N(mean = M(c), std^{2} = S^{2}), where S[is represented by the Variance function below] = U(min = 0.1, max = 0.3)
CellAdmixture <- function(value){
  Pt <- runif(1,0.3,0.7)
  Result <- log2( (value * Pt + 2 *(1-Pt))/2 )
  return(Result)
}
Variance <- function(value){
  return(runif(1,0.1,0.3))
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
LinearExpression <- function(Vector,TakenProbes){
  TempVector <- c()
  # There are a total of 15 genes(false/true positive), 3 for each group(which there are 5 of)
  # They are generated in sequential order [from type I to V]
  TypeIndex <-1
  for (i in Vector){
    # Forumla
    CurrCopyNum <- GetCopyNumber(i)
    CurrValue <- 2*CurrCopyNum+2
    # When index is belongs to a gene type
    if (i %in% TakenProbes){
      # Type 1
      if(TypeIndex <= 3){
        TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
      }else if(TypeIndex <= 6 && runif(1)<0.75){
        # Type 2
        TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
      }else if(TypeIndex <= 9 && runif(1)<0.50){
        # Type 3
        TempVector <- c(TempVector,runif(1,CurrValue,CurrValue+0.5))
      }else if(TypeIndex <= 12){
        # Type 4
        TempVector <- c(TempVector,rnorm(1,4,sqrt(0.5)))
      }else if(TypeIndex <= 15){
        # Type 5
        TempVector <- c(TempVector,rnorm(1,12,sqrt(0.5)))
      }else{
        # Normal Expression - due to the type I-III not having the proper runif value
        TempVector <- c(TempVector,rnorm(1,CurrValue,sqrt(2)))
      }
      TypeIndex <- TypeIndex + 1
    }else{
      # Normal Expression
      CurrValue <- 6
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
    # Normal Expression
    TempVector <- c(TempVector,rnorm(1,CurrValue,sqrt(2)))
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
GenerateComplementaryGeneExpression <- function(Vector,TakenProbes,Model){
  TempVector <- c()
  
  if(Model == "Linear"){
    
    TempVector <- c(TempVector,LinearExpression(Vector,TakenProbes))
    
  }else if(Model == "Sigmoid"){
    
    TempVector <- c(TempVector,SigmoidExpression(Vector))
    
  }else if (Model =="StepWise"){
    TempVector <- c(TempVector,StepWiseExpression(Vector))
  }
  
  return(TempVector)
  
}


# Return the probe regions(the x values in the graph) which have copy number abberations
CopyErrorCheck <- function(ProbeVector){
  TempVector <- c()
  for(i in 1:length(ProbeVector)){
    if(round(GetCopyNumber(ProbeVector[i])) != 2){
      TempVector <- c(TempVector,c(i))
    }
    
  }
  
  return(TempVector)
} 


  

  















# Generate a graph for each simulated patient
for(k in 1:NumberOfPatients){
  # Get copy number simulation by reading lines and giving each probe its appropriate value
  CopyNumberExpression <- c()
  for(Line in ProbeLocation){
    # Split tab-delimited line into two(the first value is the number of probes with a specific CGH value
    # and the second value is the actual CGH ratio value assigned to those probes)
    SplitLine <- unlist(strsplit(Line, split = "\t"))
    # add the appropriate amount and type of genes
    CopyNumberExpression <- c(CopyNumberExpression,c(rnorm(as.integer(SplitLine[1]),CellAdmixture(InputTranslationTable[[SplitLine[2]]]),Variance())))
  }  
  
  # Vector that hold the x-values of probes that have copy number abberations
  ErrorProbes = CopyErrorCheck(CopyNumberExpression)
  # The first 9 are always assigned to types I-III
  TIProbeLocation <- c(ErrorProbes[1], ErrorProbes[2], ErrorProbes[3])
  TIIProbeLocation <- c(ErrorProbes[4], ErrorProbes[5], ErrorProbes[6])
  TIIIProbeLocation <- c(ErrorProbes[7], ErrorProbes[8], ErrorProbes[9])
  TIVProbeLocation <- c(ErrorProbes[10], ErrorProbes[11], ErrorProbes[12])
  TVProbeLocation <- c(ErrorProbes[12], ErrorProbes[13], ErrorProbes[14])
  TakenProbes <- c(TIProbeLocation,TIIProbeLocation,TIIIProbeLocation,TIVProbeLocation,TVProbeLocation)
  
  
  
  
  # Generate Gene Expression
  GeneExpression <- GenerateComplementaryGeneExpression(CopyNumberExpression,TakenProbes,Model)
  
  
  # To help decide where to put typeI-V genes If User does not specifiy anything in terms of location
  ProbeRegionsWithCopyNumberProblems <- CopyErrorCheck(CopyNumberExpression)  # For Type I - III

  
 
  SizeAdj <- c()
  # Size adjustment for normal and typeI-V genes
  for(i in 1:length(CopyNumberExpression)){
    # if a probe is a typeI-V gene then make it double the size of a normal gene expression dot
    if(i %in% TakenProbes){
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
    if(i %in% TIProbeLocation){

 
      GeneEpxrCol <- c(GeneEpxrCol,"red")
      BackgroundCol <- c(BackgroundCol,"red")
    }else if(i %in% TIIProbeLocation ){
      # If a probe is type 2 then make it a yellow dot
      GeneEpxrCol <- c(GeneEpxrCol,"yellow")
      BackgroundCol <- c(BackgroundCol,"yellow")
    }else if(i %in% TIIIProbeLocation){
      # If a probe is type 1 then make it a orange dot
      GeneEpxrCol <- c(GeneEpxrCol,"orange")
      BackgroundCol <- c(BackgroundCol,"orange")
    }else if(i %in% TIVProbeLocation){
      # If a probe is type 4 then make it a purple dot and change it's value(according to the source paper)
      GeneEpxrCol <- c(GeneEpxrCol,"purple")
      BackgroundCol <- c(BackgroundCol,"purple")
    }else if(i %in% TVProbeLocation){
      # If a probe is type 5 then make it a purple dot and change it's value(according to the source paper)
      GeneEpxrCol <- c(GeneEpxrCol,"blue")
      BackgroundCol <- c(BackgroundCol,"blue")
    }else{
      # If a probe is a regular gene expression value then make it gray
      GeneEpxrCol <- c(GeneEpxrCol,"gray")
      BackgroundCol <-c(BackgroundCol,"gray")
    }
  }
  
  # Create the x-coordinates for the graph
  xcor <- 1:length(CopyNumberExpression)
  print("Generating Graphs")
  # Plot the first part of the graph(cgh ratio/blue dots)
  plot(xcor,CopyNumberExpression,col = "blue",ylim=c(-2,3), cex=0.3 ,ylab="CGH Ratio Values", xlab="Probes/Genes")
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
}


