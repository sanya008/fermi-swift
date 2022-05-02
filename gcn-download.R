library(checkpoint)
checkpoint("2022-01-01")

#install.packages("stringr")                        # Install stringr package
library("stringr")                                 # Load stringr package



drive <- "e:"
dir <- "/Research/Fermi_Swift_comp/__rev7/Swift-Fermi-comp/gcns/"

# Importing gcn # list

infile <- paste(drive,dir,"gcn.txt",sep="")
gcn <- read.table(infile,header=F,sep=",")
names(gcn) <- c("gcn_number","nul","grb_id")
gcn$nul <- NULL


# Importing BAT data

infile <- paste(drive,dir,"Swift_BAT.txt",sep="")

swg <- read.table(infile,header=T)


infile <- paste(drive,dir,"../BAT_T90.txt",sep="")
batT90e <- read.table(infile,header=T)


# Downloading gcns

#for(i in gcn$V1){
#  download.file(paste("https://gcn.gsfc.nasa.gov/gcn3/",i,".gcn3",sep=""),paste(dir,i,".gcn3",sep=""))
#}

#filelist <- list.files(path="e:/Research/Fermi_Swift_comp/__rev7/",pattern=".gcn3")



#T90=0
#T90err=0

for (j in c(1:length(gcn$gcn_number))){
  textforparse <- readChar(paste(drive,"/Research/Fermi_Swift_comp/__rev7/",
                                 gcn$gcn_number[j],".gcn3",sep=""),file.info(paste(drive,"/Research/Fermi_Swift_comp/__rev7/",gcn$gcn_number[j],".gcn3",sep=""))$size)
  textforparse <- gsub("[\r\n]", "", textforparse)
  textforparse <- gsub("\\+/\\-", "+-", textforparse)
  loc <- str_locate_all(pattern = "T90", textforparse)[[1]][1]
#  loc <- str_locate_all(pattern = "T90 [\\(]", textforparse)[[1]][1]
  loc_end <- str_locate_all(pattern = "sec", substring(textforparse,loc,nchar(textforparse)))[[1]][1]
 if (loc_end > 50 | is.na(loc_end)) {  
   T90[j] <- NA
   T90err[j] <- NA
   } else {
     #print(c(loc,loc+loc_end))
      loc <- loc + str_locate_all(pattern = "s", substring(textforparse,loc,loc+loc_end))[[1]][1]
      T90[j] <- as.numeric(gsub(".*?([0-9.]+).*", "\\1", substring(textforparse,loc,nchar(textforparse))))   
      loc <- loc + str_locate_all(pattern = "[\\+][\\-]", substring(textforparse,loc,nchar(textforparse)))[[1]][1]
      T90err[j] <- as.numeric(gsub(".*?([0-9.]+).*", "\\1", substring(textforparse,loc,nchar(textforparse))))   
      if (is.na(T90err[j])){T90[j]=NA}
      }}

gcn$T90 <- as.numeric(T90[1:length(T90)])
gcn$T90err <- as.numeric(T90err[1:length(T90err)])

gcn$T90[gcn$gcn_number==8019] <- 120
gcn$T90err[gcn$gcn_number==8019] <- 75

outfile <- drive
outfile <- paste(drive,dir,"../swift_t90_error.txt",sep="")

write.table(gcn,outfile,quote=F,row.names=F)


#test <- merge(swg,batT90e)
#plot(log10(test$T90),log10(test$T90cat))



