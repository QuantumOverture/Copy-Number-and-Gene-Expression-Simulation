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
Positives$RawLines <- readLines(PositiveFile)

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
colnames(PatientGeneData) <- seq(1,NumOfPatients,1)

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
GE_mine$info$chromosome <- rep.int(4,ProbesPerPatient)
GE_mine$info$position <- seq(1,ProbesPerPatient,1)


# Prepare GE_ref
PatientGeneREFData <- matrix(nrow=ProbesPerPatient,ncol=2)
colnames(PatientGeneREFData) <- c(1,2)
GE_mine_ref <- c()
GE_mine_ref$data <- c()
for(i in 1:2){
  for(j in 1:ProbesPerPatient){
    PatientGeneREFData[j,i] <- rnorm(1,6,sqrt(2))
  }
}
# Add data to GE_ref dataframe "data"
mode(PatientGeneREFData) = "numeric"
GE_mine_ref$data <- data.frame(PatientGeneREFData)
colnames(GE_mine_ref$data) <- c(1,2)

# Intialize the GE$info data frame
GE_mine_ref$info <- c()
GE_mine_ref$info$chromosome <- rep.int(4,ProbesPerPatient)
GE_mine_ref$info$position <- seq(1,ProbesPerPatient,1)

#=========================

# Intialize the GE$data data frame
CN_mine <- c()
CN_mine$data <- c()

# Intialize and fill data in matrix that will become GE_mine$data
PatientCNData <- matrix(nrow=ProbesPerPatient,ncol=NumOfPatients)
colnames(PatientGeneData) <- seq(1,NumOfPatients,1)

for(i in 1:NumOfPatients){
  for(j in 1:ProbesPerPatient){
    PatientCNData[j,i] <- CGHValues[j + (ProbesPerPatient*(i-1))]
  }
}

# Add matrix/patient data(in correct format) to the dataframe
mode(PatientCNData) = "numeric"
CN_mine$data <- data.frame(PatientCNData)
colnames(CN_mine$data) <- seq(1,NumOfPatients,1)


# Intialize the GE$info data frame
CN_mine$info <- c()
CN_mine$info$chromosome <- rep.int(4,ProbesPerPatient)
CN_mine$info$position <- seq(1,ProbesPerPatient,1)




# Intialize the GE$data data frame
CN_ref_mine <- c()
CN_ref_mine$data <- c()

# Intialize and fill data in matrix that will become GE_mine$data
PatientCNREFData <- matrix(nrow=ProbesPerPatient,ncol=2)
colnames(PatientCNREFData) <- c(1,2)
# Change when neccessary
offset <- 10

for(i in 1:2){
  for(j in 1:ProbesPerPatient){
    Pt <- runif(1,0.3,0.7)
    PatientCNREFData[j,i] <- ((2 * Pt + 2 *(1-Pt))/2 ) + offset
  }
}

# Add matrix/patient data(in correct format) to the dataframe
mode(PatientCNREFData) = "numeric"
CN_ref_mine$data <- data.frame(PatientCNREFData)
colnames(CN_ref_mine$data) <- c(1,2)


# Intialize the GE$info data frame
CN_ref_mine$info <- c()
CN_ref_mine$info$chromosome <- rep.int(4,ProbesPerPatient)
CN_ref_mine$info$position <- seq(1,ProbesPerPatient,1)



# Check segment function for any missed info

CN_mine_ref_segmented <- segmentCN(CN_ref_mine, GE_mine_ref$info, chr=4:7, algorithm="CBS")

CN_mine_segmented <- segmentCN(CN_mine, GE_mine$info, chr=4:7, algorithm="CBS")
# Current Tasks:
# What is the signal conversion for CGH things
# What should chromosome(in GE) look like



for(i in 1:length(Positives$RawLines)){
  if(i == 1){
    
  }
  
} 



results <- edira(GE,CN_segmented,GE_ref,CN_ref_segmented)




