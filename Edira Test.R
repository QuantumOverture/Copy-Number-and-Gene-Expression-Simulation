dir = "C:/Users/ismai/Desktop/CopyNumAndGeneExpSimulation/"

# Load in AMLData:
# load("C:/Users/ismai/Downloads/edira_1.2.3/edira/data/AMLdata.RData")

# Load Edira function:
# lazyLoad("C:/Users/ismai/Downloads/edira_1.2.3/edira/R/edira")

args = commandArgs()

# Fix Formating Issues ya foo'
CGHFile <- "CGH.txt" #args[1]
GeneFile <- "GeneExp.txt" #args[2]
PositiveFile <- "Postives.txt" #args[3]

# Only works for non=average output for patients
CGHValues <- readLines(CGHFile)
GeneValues <- readLines(GeneFile)
Positives <- as.numeric(readLines(PositiveFile))

# Get Num Of Patients & Probes per patient and remove from vectors
NumOfPatients <- unlist(strsplit(CGHValues[1],split=" "))
NumOfPatients <- as.numeric(NumOfPatients[4])
ProbesPerPatient <- unlist(strsplit(CGHValues[2],split=" "))
ProbesPerPatient <- as.numeric(ProbesPerPatient[4])

CGHValues <- CGHValues[3:length(CGHValues)]
GeneValues <- GeneValues[3:length(GeneValues)]

# Intialize the GE$data data frame
GE_mine <- c()
GE_mine$data <- c()

# Intialize and fill data in matrix that will become GE_mine$data
PatientGeneData <- matrix(nrow=ProbesPerPatient,ncol=NumOfPatients)
colnames(PatientData) <- seq(1,NumOfPatients,1)

for(i in 1:NumOfPatients){
  for(j in 1:ProbesPerPatient){
    PatientGeneData[j,i] <- GeneValues[j + (ProbesPerPatient*(i-1))]
  }
}

# Add matrix/patient data(in correct format) to the dataframe
mode(PatientGeneData) = "numeric"
GE_mine$data <- data.frame(PatientGeneData)
colnames(GE_mine$data) <- seq(1,NumOfPatients,1)


# Intialize the GE$info data frame
GE_mine$info <- c()
GE_mine$info$chromosome <- c()
GE_mine$info$position <- seq(1,ProbesPerPatient,1)



GE_mine_ref <- 


# Current Tasks:
# What is the signal conversion for CGH things
# What should chromosome(in GE) look like







results <- edira(GE,CN_segmented,GE_ref,CN_ref_segmented)


