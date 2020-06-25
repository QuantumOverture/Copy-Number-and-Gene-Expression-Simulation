# Optional Seed (uncomment the next line and add integer in place of x):
# set.seed(x)
library(list)
# args hold the command line arguments neccessary for the program
args = commandArgs()
# Model connects copy number to gene expression{Sigmoid,Linear,Stepwise}
Model <- "Linear" #args[3]
# The number of graphs/patients we will generate
NumberOfPatients <- 1 #args[4]
# File that holds information about how many probes there are and what their regions will be(gain,loss,normal and variations - for each if they exist)
ProbeLocationFile <- "./ProbeLocation.txt"#args[5]


# Setup functions and values used in CGH value assignment:
# Formulae:
# Cell admixture is --> M(c) = log2 [(c · Pt + 2 · (1 − Pt))/2], where Pt = U(0.3, 0.7)
# N(µ = M(c), σ^{2} = S^{2}), where S[is represented by the Variance function below] = U(min = 0.1, max = 0.3)
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
InputTranslationTable = list()
# Gain Narrow Medium Aberration1
InputTranslationTable("GNM1") = 8
# Gain Narrow Medium Aberration2
InputTranslationTable("GNM2") = 6
# Loss Narrow Medium Aberration
InputTranslationTable("LNMA") <- 1
# Gain Large High Aberration
InputTranslationTable("GLHA") <- 10
# Loss Large High Aberration
InputTranslationTable("LLHA") <- 0
# Gain Broad Low Abberation 
InputTranslationTable("GBLA") <- 3
# Normal Cell
InputTranslationTable("NOCE") <- 2

CopyNumberExpression <- c()
for(Line in ProbeLocation){
  SplitLine <- unlist(strsplit(Line, split = "\t"))
  CopyNumberExpression <- c(CopyNumberExpression,c(rnorm(as.integer(SplitLine[1]),CellAdmixture(InputTranslationTable[SplitLine[2]]),Variance())))
}
# Left here
print(CopyNumberExpression)
