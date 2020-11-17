require("TDA")
require("stringr")
options(scipen=999)
setwd("C:/Users/ismai/Desktop/Script4Update/Research/TAaCGH")
# Command Line Arguements
args <- commandArgs()
dataSet <- "horlings" #args[4]
homDim <- 2 #as.numeric(args[5])
partNum <- 101 #as.numeric(args[6])
epsIncr <- 0.01 #round(as.numeric(args[7]), 3)
action <- "sect" #args[8]


decimalSpotsOfIncr <- nchar(strsplit(as.character(epsIncr),"[.]")[[1]][2])

#### TO DO LIST ####

# Add build_individual_profile for my simulated data

####################


#################IMPORTANT INTERNAL FUNCTIONS#######################

build_individual_profile <- function(CGHdata , patient , start , finish){
  tProfile <- CGHdata[patient,(start+1):(finish + 1)]
  return(tProfile)
}

build_cloud <- function(CGHprofile , dim){
  tempPoint <- c()
  tCloud <- c()
  # Go through the number of points for that chromosome arm
  for(i in 1:ncol(CGHprofile)){
    # Have as many coordinates as dimensions
    for(j in 1:dim){
        # This ensures proper wraparound.
        if(i+j-1 <= ncol(CGHprofile)){
          tempPoint <- c( tempPoint,round(as.numeric(CGHprofile[ 1,(i+j-1)]),digits=6) )
        }else{
          tempPoint <- c( tempPoint,round(as.numeric(CGHprofile[ 1,( (i+j-1) %% ncol(CGHprofile) )]),digits=6) )
        }
      #tempPoint <- c( tempPoint,round(as.numeric(CGHprofile[ 1,( (i+j-1) %% ncol(CGHprofile) )]),digits=6) )
    }
  tCloud <- rbind(tCloud,tempPoint)
  tempPoint <- c()
  }
  
return(tCloud)
}


Maximum_Distance_Between_Points <- function(cloud, cloudDim){
  MaxDist = -1
  for(Row in 1:(length(cloud)/cloudDim)){
    for(SecondRow in (Row):(length(cloud)/cloudDim)){
      if(SecondRow == Row){
        next
      }
      TempDist = 0
      for(dim in 1:cloudDim){
        TempDist = TempDist + (cloud[Row,dim]-cloud[SecondRow,dim])**2
      }
      
      if(sqrt(TempDist) > MaxDist){
        MaxDist = sqrt(TempDist)
      }
    }
  }
  
  return(MaxDist)
}

# Start here



ClearAndSortBarCode <- function(BarCode,EpsLimit,epsIncr){
  Result <- c()#data.frame("Dimension"=c(1),"Birth"=c(),"Death"=c("-"))
  Dim0 <- c()
  Dim1 <- data.frame("Dimension"=c(1),"Birth"=c(1),"Death"=c(1))
  for(Row in 1:nrow(BarCode)){

    Dimension = as.numeric(BarCode[Row,"dimension"])
    Birth = trunc(as.numeric(BarCode[Row,"Birth"])*10^decimalSpotsOfIncr)/10^decimalSpotsOfIncr
    Death = trunc(as.numeric(BarCode[Row,"Death"])*10^decimalSpotsOfIncr)/10^decimalSpotsOfIncr
    if(Birth != Death){
      if(Dimension == 1){
        Dim1<- rbind(Dim1,c(Dimension,Birth,Death))

      }else{
        Dim0 <- c(Dim0,Death)
      }
      
    }
  
  }
  
  
  Dim1 <- Dim1[-1,]
  if(nrow(Dim1)>1){
  Dim1$"Birth" <- as.numeric(as.character(Dim1$"Birth"))
  Dim1 <- Dim1[order(Dim1$"Birth"),]
  }
  Dim0 <- sort(Dim0)
  
  # Set to correct Interval
  Interval = 1
  for(Index in 1:(length(Dim0)-1)){
    if(Dim0[Index] <= (Interval * epsIncr)){
      Dim0[Index] = Interval * epsIncr 
    }else{
      while(Dim0[Index] > (Interval * epsIncr)){
        Interval = Interval + 1
      }
      Dim0[Index] = Interval * epsIncr 
    }
  }
  
  if(nrow(Dim1) !=  0){
  Interval = 1
  for(row in 1:nrow(Dim1)){
    if(Dim1[row,2] <= (Interval * epsIncr)){
      Dim1[row,2] = Interval * epsIncr 
    }else{
      while(Dim1[row,2] > (Interval * epsIncr)){
        Interval = Interval + 1
      }
      Dim1[row,2] = Interval * epsIncr 
    }
  }
  
  Interval = 1
  for(row in 1:nrow(Dim1)){
    if(Dim1[row,3] <= (Interval * epsIncr)){
      Dim1[row,3] = Interval * epsIncr 
    }else{
      while(Dim1[row,3] > (Interval * epsIncr)){
        Interval = Interval + 1
      }
      Dim1[row,3] = Interval * epsIncr 
    }
  }
  }
  # Betti 0 first
  for(num in Dim0){
    if(num != EpsLimit){
      Result <- rbind(Result,c(0,0,num))
    }else{
      Result <- rbind(Result,c(0,0,"inf"))
    }
  }
  
  
  
  # Betti 1 second
  if(nrow(Dim1)>1){
    for(row in 1:nrow(Dim1)){
      Result <- rbind(Result,c(1,Dim1[row,2],Dim1[row,3]))
    }
  }
  
  return(Result)
}


# Assuming Result ==> Betti Chunk Of BarCode
CreateJagLine <- function(Result){
  # Count Betti 0
  TempVal = "0"
  NumOfBetti0 = 1
  Betti0List <- c()
  while(TempVal != "1" && NumOfBetti0 <= nrow(Result)){
    TempVal = Result[NumOfBetti0,1]
    NumOfBetti0 = NumOfBetti0 + 1
    
  }
  
  MaxEps = as.numeric(Result[NumOfBetti0-2,3])
  StartEps = epsIncr
  Counts <- c(NumOfBetti0-1)
  TotalAlive = NumOfBetti0-1
  while(StartEps < MaxEps){
    NumOccurence = sum(as.character(StartEps) == Result[,3])
    if(NumOccurence > 0){
      TotalAlive = TotalAlive - sum(as.character(StartEps) == Result[,3])
      Counts <- c(Counts,TotalAlive)
    }else{
      Counts <- c(Counts,tail(Counts,n=1))
    }
    StartEps = StartEps + epsIncr 
  }
  Counts <- c(Counts,1)

}

####################################################################

# Setup correct file paths
setwd(file.path(getwd(),".."))
dictFile <- paste(dataSet,"_",action,"_dict_",as.character(partNum),".txt",sep="")
dictPath <- file.path(getwd(),"Data",dataSet,action,dictFile)

dataFile <- paste(dataSet,"_data.txt",sep = "")
dataPath <- file.path(getwd(),"Data",dataSet,dataFile)

# Read data and dictionary and remove unneccessary rows
inputList <- read.table(dataPath, sep = '\t',header = FALSE,stringsAsFactors=FALSE) 
inputList <- tail(inputList,-2)
dictList <- read.table(dictPath, sep = '\t',header = FALSE,stringsAsFactors=FALSE)
dictList <- dictList[-1,]

# Using dimension 2 as window
cloudDim <- 2

# Result Path
reasultsPath <- file.path(getwd(),"..","Results",dataSet,action,paste(as.character(cloudDim),"D",sep = ""),"Homology")


for(i in 1:nrow(dictList)){
  chr = as.numeric(dictList[i,1])
  arm = dictList[i,2]
  
  beg = as.numeric(dictList[i,3])
  end = as.numeric(dictList[i,4])
  
  seg = as.numeric(dictList[i,6])

  
  
  for(j in 1:nrow(inputList)){
    profile = build_individual_profile(inputList,j,beg,end)
    cloud = build_cloud(profile, cloudDim)
    
    MaxDist = Maximum_Distance_Between_Points(cloud, cloudDim) * 0.9
    MaxDIM = 1 # 1 = Loops/holes
    MaxSCA = 100 #  
    
    
    Diag <- ripsDiag(X = cloud, maxdimension = MaxDIM, maxscale =MaxSCA,
                     library = "Dionysus", printProgress = TRUE)
    
    
    BarCode <- Diag[[1]]

    ClearedAndSortedBarCode <- ClearAndSortBarCode(BarCode,MaxSCA,epsIncr)
    # Start for betti 0
    
    CreateJagLine(ClearedAndSortedBarCode)
    # Make sure cloud graph is working properly
    
  }
}







