#
# title: Comparison of Fermi and Swift GRB Data
# author: Balázs, L.G., Pinter, S.
# date (as of rev 7): 2022.05.02
#

library(checkpoint)
checkpoint("2022-01-01")

point_to_mm <- 0.3527777778
fontsize=18
library(ggplot2)
library(cowplot)
library(celestial)
library(GGally)
library(ggExtra)
library(ggforce)
library(stringr)



##### Basic Data Input 1 Swift-BAT #####

# Specify drive and directory
drive <- "e:"
dir <- "/Research/Fermi_Swift_comp/__rev7/Swift-Fermi-comp/"

# Importing BAT data
infile <- paste(drive,dir,"Swift_BAT.txt",sep="")
swg <- read.table(infile,header=T)

# Setting fluence to erg*cm^-2
swg$flu <- swg$flu/1e7

# Importing BAT T90 error
infile <- paste(drive,dir,"swift_t90_error.txt",sep="")
swgT90e <- read.table(infile,header=T,sep=" ")
colnames(swgT90e)[3] <- "T90gcn"
colnames(swgT90e)[2] <- "Name"

# Merge swg & batT90e
swg <- merge(swg,swgT90e,by="Name",all.x = T)

# Checking T90 differences between gcns and database
which(swg$T90!=swg$T90gcn)

# Importing Swift trigger time data
infile <- paste(drive,dir,"Swift_Trtime.txt",sep="")
swTr <- read.table(infile,header=T)

# Merge swg & trigger time
swg <- merge(swg,swTr,by="Name",all.x = T)

#swg<-swg[!is.na(swg$TrTime), ]


# Dropping the pre-Fermi time observations
swg <- swg[as.numeric(substr(swg$Name,1,6)) > 80714,]



# Year, month and day of trigger
syear <- substring(swg$Name,1,2)
syear <- paste("20",syear,sep="")
smonth <-substring(swg$Name,3,4)
sday <-  substring(swg$Name,5,6)
sdate <- paste(syear,"-",smonth,"-",sday,sep="")

# Swift satellite launch date
start <- "2004-11-20 17:16:01"

# Date of trigger
swTrdate <- paste(sdate,swg$TrTime,sep=" ")

# Elapsed time of Trigger from the start in days
swTrtime <-as.numeric(difftime(swTrdate,start,units="day"))

# Converting Ra, Dec into Descartes coordinates
sx <- cos(swg$RA*pi/180)*cos(swg$Dec*pi/180)
sy <- sin(swg$RA*pi/180)*cos(swg$Dec*pi/180)
sz <- sin(swg$Dec*pi/180)

# Creating a position-time dataframe 
sxyzt <- data.frame(sx,sy,sz,swTrtime)

ggplot(swg,aes(RA,Dec)) +
  theme_bw() +
  geom_point(data=swg,size=0.9, shape=16,color='blue') + 
  coord_map(projection="aitoff",clip = "off", 
            orientation=c(90,180,0), 
            xlim=c(0,360))+#, ylim=c(0,0)) + 
  scale_y_continuous(breaks=(-3:3)*30,limits=c(-90,90)) + 
  scale_x_continuous(breaks=(0:4)*90,limits=c(0,360),expand = c(1,0)) +
  theme(legend.text = element_text(size = fontsize*0.8), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        strip.text = element_text(size = fontsize*0.8, colour='black'),
        axis.ticks = element_blank(),#element_line(colour = 'black'),
        axis.text = element_blank(),
        panel.border = element_blank()) + #element_text(size = fontsize*0.8, colour='black')) + 
  labs(subtitle="", y="Declination", x="Right Ascension", title="")

#outfile <- paste(drive,dir,"figs/fig01a.pdf",sep="")
#ggsave(outfile,width = 6, height = 4,dpi=300)


##### Basic Data Input 1 Fermi-GBM #####

infile <- paste(drive,dir,"GBM21.txt",sep="")
gbm <- read.table(infile,header=T,sep="|")

# Dropping observations after 2020/06/31
gbm <- gbm[as.numeric(substr(gbm$Fname,4,9)) < 200631,]

# Ellapsed time of trigger from the Swift start in days
feTrtime <-as.numeric(difftime(gbm$trigger_time,start,units="day"))

# Converting Ra, Dec into Descartes coordinates
gbm$ra  <- hms2deg(as.character(gbm$ra),sep=" ")
gbm$dec <- dms2deg(as.character(gbm$dec),sep=" ") 

fx <- cos(gbm$ra*pi/180)*cos(gbm$dec*pi/180)
fy <- sin(gbm$ra*pi/180)*cos(gbm$dec*pi/180)
fz <- sin(gbm$dec*pi/180)

# Creating a position-time dataframe 
fxyzt <- data.frame(fx,fy,fz,feTrtime)

ggplot(gbm,aes(ra,dec)) +
  theme_bw() +
  geom_point(data=gbm,size=0.9, shape=16,color = "darkred") + 
  coord_map(projection="aitoff",clip = "off", 
            orientation=c(90,180,0), 
            xlim=c(0,360))+#, ylim=c(0,0)) + 
  scale_y_continuous(breaks=(-3:3)*30,limits=c(-90,90)) + 
  scale_x_continuous(breaks=(0:4)*90,limits=c(0,360),expand = c(1,0)) +
  theme(legend.text = element_text(size = fontsize*0.8), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        strip.text = element_text(size = fontsize*0.8, colour='black'),
        axis.ticks = element_blank(),#element_line(colour = 'black'),
        axis.text = element_blank(),
        panel.border = element_blank()) + #element_text(size = fontsize*0.8, colour='black')) + 
  labs(subtitle="", y="Declination", x="Right Ascension", title="")

#outfile <- paste(drive,dir,"figs/fig01b.pdf",sep="")
#ggsave(outfile,width = 6, height = 4,dpi=300)


##### Computing nearest neighbour distance #####

library(FNN)

fsknn1 <- get.knnx(fxyzt,sxyzt,k=1)
nndist <- fsknn1$nn.dist
nnind <-  fsknn1$nn.index
nnres <- data.frame(nndist,nnind)

# Plotting histogram of nndist
ggplot(nnres[nndist<10,], aes(x=log10(nndist))) +
  geom_histogram(bins=40,fill="#00bfc4",alpha=0.5,colour="black") + ggtitle("NN distances between Swift and Fermi GRBs") + theme_bw() + theme(plot.title = element_text(size=12,hjust = 0.5)) +  geom_vline(aes(xintercept= -2.5), color="red",             linetype="dashed") + geom_text(x=-4.5, y=200, label="Swift-Fermi coincidences", color="red",size=4.0) +
  labs(subtitle="", y="count", 
       x="lg(NN distances) [days]", 
       title="") + 
  theme_bw() + 
  theme(legend.text = element_text(size = fontsize*0.8), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        legend.justification=c(1,0),
        strip.text = element_text(size = fontsize*0.8, colour='black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text = element_text(size = fontsize*0.8, colour='black')     )

outfile <- paste(drive,dir,"figs/fig01-nndist.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)


# Numbers of Swift and Fermi GRBs
Swno <- dim(swg)[1]
Feno <- dim(gbm)[1]
#
# Serial numbers for jointly detected Swift and Fermi GRBs
#
FeSwind <- nnind[log10(nndist) < -2.5]
SwFeind <- c(1:Swno)[log10(nndist) < -2.5]
noFeSw <- length(FeSwind)
noSwFe <- length(SwFeind)

# Creating data frame for Swift and Fermi "couples"
BATcoup <- swg[SwFeind,]
GBMcoup <- gbm[FeSwind,2:32]
SwFe_couples <- data.frame(BATcoup,GBMcoup)

# Writing SwFe_couples into a file
outfile <- paste(drive,dir,"SwFe_couples.txt",sep="")
write.table(SwFe_couples,outfile,quote=F,row.names=F)

# Testing the Sw-Fe matching
Swtest <- as.numeric(substring(BATcoup$Name,1,6))
Fetest <- as.numeric(substring(GBMcoup$Fname,4,9))
Swtest-Fetest

# Identifying the "widows"
Swidow <- rep("yes",Swno)
Swidow[SwFeind] <- rep("no",noSwFe)

Fewidow <-rep("yes",Feno)
Fewidow[FeSwind] <- rep("no",noFeSw)

# Creating BATwidow and GBMwidow frames
BATwidow <- swg[Swidow=="yes",]
GBMwidow <- gbm[Fewidow=="yes",]

# Writing frames created
outfile <- paste(drive,dir,"BATwidows.txt",sep="")
write.table(BATwidow,outfile,quote=F,row.names=F)
outfile <- paste(drive,dir,"GBMwidows.txt",sep="")
write.table(GBMwidow,outfile,quote=F,row.names=F)


##### Comparing Swift "couples" and "widows" #####


Swcouples <- SwFe_couples[,1:16]
ncoup <- dim(Swcouples)[1]
status <- rep("couples",ncoup)
Sptype <- Swcouples[,11]
Swcouples <-data.frame(Swcouples,status,Sptype)

nwid <- dim(BATwidow)[1]
indwid <- sample(1:nwid,ncoup)
status <- rep("widows",ncoup)
Sptype <- BATwidow[indwid,11]
Swidows <- data.frame(BATwidow[indwid,],status,Sptype)

#Swidows <- Swidows[as.numeric(substr(Swidows$Name,1,6)) > 80714,]

Swcowid <- rbind(Swcouples,Swidows)

# Displaying data in pairplot
status <- Swcowid$status
Sptype <- Swcowid$Sptype

# Computing and plotting logarithmic variables
lSwcowid <- data.frame(log10(Swcowid[,c(5,6,8,10 )]),status,Sptype)
colnames(lSwcowid)[2]<- "fluence"
colnames(lSwcowid)[3]<- "peak flux"
colnames(lSwcowid)[4]<- "phi"

ggpairs(lSwcowid,legend=1,columns=1:4,  ggplot2::aes(colour=status),  
        lower = list(continuous = wrap("points", size=0.3,shape=19,alpha=0.5)), 
        diag = list(continuous = wrap("densityDiag", alpha=0.5)), 
        upper = list(continuous = wrap("cor", size = 2.5))) + 
  theme_bw() + removeGrid() + 
  ggtitle('Distribution of BAT "couples" and "widows"') + 
  theme(legend.position="bottom",legend.title = element_blank()) + 
  theme(plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

outfile <- paste(drive,dir,"figs/fig06-dist-bat.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)



# Testing differences between Swift "couples" and "widows"
T90coup <- Swcowid$T90[status=="couples"]
T90wid <- Swcowid$T90[status=="widows"]
ks.test(T90coup,T90wid)
#
flucoup <- Swcowid$flu[status=="couples"]
fluwid <- Swcowid$flu[status=="widows"]
ks.test(flucoup,fluwid)
#
P1coup <- Swcowid$P1[status=="couples"]
P1wid <- Swcowid$P1[status=="widows"]
ks.test(P1coup,P1wid)
#
Phicoup <- Swcowid$Phi[status=="couples"]
Phiwid <- Swcowid$Phi[status=="widows"]
ks.test(Phicoup,Phiwid)

# Computing cross correlations between status and Sp type
tb <- table(status,Sptype)
tb
#
summary(tb)
#
# Preparing Data for barplot
#
colnames <- rep(names(tb[1,]),dim(tb)[1])
value   <- as.numeric(tb)
mstatus <- rep(row.names(tb),dim(tb)[2])
# Sptyp   <- colnames[c(1,4,2,5,3,6)]
Sptyp   <- colnames[c(1,3,2,4)]
sd <-data.frame(value, mstatus,Sptyp)

# Displaying barplot
ggplot(sd, aes(fill=mstatus, y=value, x=Sptyp)) + 
  geom_bar(position="dodge", stat="identity") + ggtitle("Distribution of Spectral types") + theme_bw() +
  scale_x_discrete(name="Spectral type") + scale_y_continuous(name="Number of GRBs") + theme(plot.title = element_text(size=12,hjust = 0.5))
#ggsave("F:/Research/Fermi_Swift_comp/newfigs/06.pdf",width = 6, height = 4)




# Outliers' diagnostics of BAT data
# Using boxplots
bpg1 <- ggplot(lSwcowid,aes(T90)) + geom_boxplot(fill="blue",alpha=0.5,colour="black") + theme_bw() + coord_flip()+ theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+ theme(axis.title.y = element_text(size=fontsize*1.2),axis.text.y = element_text(size=fontsize*0.8))
bpg2 <- ggplot(lSwcowid,aes(fluence)) + geom_boxplot(fill="blue",alpha=0.5,colour="black") + theme_bw() + coord_flip()+ theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+ theme(axis.title.y = element_text(size=fontsize*1.2),axis.text.y = element_text(size=fontsize*0.8))
bpg3 <- ggplot(lSwcowid,aes(`peak flux`)) + geom_boxplot(fill="blue",alpha=0.5,colour="black") + theme_bw() + coord_flip()+ theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+ theme(axis.title.y = element_text(size=fontsize*1.2),axis.text.y = element_text(size=fontsize*0.8))
bpg4 <- ggplot(lSwcowid,aes(phi)) +geom_boxplot(fill="blue",alpha=0.5,colour="black") + theme_bw() + coord_flip() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()) + theme(axis.title.y = element_text(size=fontsize*1.2),axis.text.y = element_text(size=fontsize*0.8))

plot_grid(bpg1,bpg2,bpg3,bpg4,ncol=4,nrow=1)

outfile <- paste(drive,dir,"figs/fig03a-outlier-bat.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)



swp1 <- boxplot.stats(lSwcowid$T90) 
swp2 <- boxplot.stats(lSwcowid$fluence)
swp3 <- boxplot.stats(lSwcowid$`peak flux`)
swp4 <- boxplot.stats(lSwcowid$phi)

inlierg1 <- lSwcowid$T90 > swp1$stats[1]-0.5 & lSwcowid$T90 < swp1$stats[5]
inlierg2 <- lSwcowid$fluence > swp2$stats[1] & lSwcowid$fluence < swp2$stats[5]
inlierg3 <- lSwcowid$`peak flux`  > swp3$stats[1] & lSwcowid$`peak flux`  < swp3$stats[5]
inlierg4 <- lSwcowid$phi > swp4$stats[1] & lSwcowid$phi < swp4$stats[5]
inlier   <- inlierg1 & inlierg2 & inlierg3 & inlierg4

# Identify outliers as NA
inlier <- ifelse(inlier,inlier,NA)

# Summary of boxplot's outlier diagnostics 
summary(data.frame(inlierg1,inlierg2,inlierg3,inlierg4, inlier))



library(MASS)

nlSwcowid <- na.omit(lSwcowid[inlier,])
Swlda <- lda(nlSwcowid[,1:4],nlSwcowid$status)
Swlda

# Computing best discriminating variable
LD1 <- Swlda$scaling
LD1 <- as.matrix(nlSwcowid[,1:4]) %*% LD1

# Testing the discriminant power of LD1, using ANOVA
df <-data.frame(nlSwcowid,LD1)
dfcor <- df[,c(5,7)]
aov.res <- aov(LD1 ~ status,data=dfcor)
summary(aov.res)

# Displaying discriminant power of LD1
pLD1 <- ggplot(df, aes(LD1, fill = status)) + 
  geom_density(alpha = 0.5) + ggtitle("Greatest couples-widows difference") +
  theme_bw()+ theme(plot.title = element_text(size=11,hjust = 0),legend.title = element_blank())

# Correlating observed variables with LD1
corLD1 <- cor(nlSwcowid[,1:4],LD1)
cordf <- data.frame(value=as.numeric(corLD1),variables=row.names(corLD1))
corLD1

# Displaying correlations p=0.05 corr =0.075
pcor <- ggplot(cordf, aes(x=variables,y=corLD1)) + 
  geom_bar(stat="identity",fill="#00bfc4",alpha=0.5,colour="black") + 
  ggtitle("Correlation of LD1 with variables") + labs(y="correlation") +
  theme_bw() + ylim(-1,1) + 
  geom_hline(yintercept=c(-0.075,0.075),linetype="dashed", color = "red", size=1) + 
  geom_text(x= 2.5, y= 0.7, label="- -: p-value = 0.05", color="red",size=4.5) + 
  theme(plot.title = element_text(size=11,hjust = 0.5))

plot_grid(pLD1,pcor,nrow=1,ncol=2)

outfile <- paste(drive,dir,"figs/fig08-sep-bat.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)




# Displaying dependence between observed and LD1 variables
pT90 <- ggplot(df,aes(x=LD1,y= T90,color=status)) + geom_point(size=0.7,alpha=0.5) + theme_bw() + ggtitle("LDA in T90")+ theme(plot.title = element_text(size=12,hjust = 0.5),legend.title = element_blank())
pflu <- ggplot(df,aes(x=LD1,y= fluence,color=status)) + geom_point(size=0.7,alpha=0.5) + theme_bw() + ggtitle("LDA in fluence")+ theme(plot.title = element_text(size=12,hjust = 0.5),legend.title = element_blank())
pP1 <- ggplot(df,aes(x=LD1,y=peak.flux,color=status)) + geom_point(size=0.7,alpha=0.5) + theme_bw() + ggtitle("LDA in peak flux")+ theme(plot.title = element_text(size=12,hjust = 0.5),legend.title = element_blank()) + labs(y="peak flux")
pPhi <- ggplot(df,aes(x=LD1,y= phi,color=status)) + geom_point(size=0.7,alpha=0.5) + theme_bw() + ggtitle("LDA in photon index")+ theme(plot.title = element_text(size=12,hjust = 0.5),legend.title = element_blank())

plot_grid(pT90,pflu,pP1,pPhi,nrow=2,ncol=2)

outfile <- paste(drive,dir,"figs/fig09-lda-bat.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)



# Displaying frequency distribution of variables
pT90 <- ggplot(df,aes(x= T90,fill=status)) +   geom_density(size=0.7,alpha=0.5) + theme_bw() + ggtitle("T90 distribution")+ theme(plot.title = element_text(size=12,hjust = 0.5))
pflu <- ggplot(df,aes(x= fluence,fill=status)) + geom_density(size=0.7,alpha=0.5) + theme_bw() + ggtitle("fluence distribution")+ theme(plot.title = element_text(size=12,hjust = 0.5))
pP1 <- ggplot(df,aes(x=peak.flux,fill=status)) + geom_density(size=0.7,alpha=0.5) + theme_bw() + ggtitle("peak flux distribution")+ theme(plot.title = element_text(size=12,hjust = 0.5))
pPhi <- ggplot(df,aes(x= phi,fill=status)) + geom_density(size=0.7,alpha=0.5) + theme_bw() + ggtitle("phi distribution")+ theme(plot.title = element_text(size=12,hjust = 0.5))

plot_grid(pT90,pflu,pP1,pPhi,nrow=2,ncol=2)

#ggsave("E:/Research/Fermi_Swift_comp/________rev6/10.pdf",width = 6, height = 4)





##### Comparing Fermi "couples" and "widows" #####

Fecouples <- SwFe_couples[,c(13+4,17+4,20+4,22+4,33+4,43+4)]
ncoup <- dim(Fecouples)[1]
status <- rep("couples",ncoup)
Fecouples <-data.frame(Fecouples,status)

Fewidows <- GBMwidow[,c(2,6,9,11,22,32)]
nwid <- dim(GBMwidow)[1]
indwid <-sample(1:nwid,ncoup)
status <-rep("widows",ncoup)
Fewidows <-data.frame(Fewidows[indwid,],status)

# Joining Fermi couples and resampled widows data
Fecowid <- rbind(Fecouples,Fewidows)
Fecowid[603,5] <- NA

# REname selectd variables
names(Fecowid)[2:6] <- c("T90","fluence","peak flux","alpha","model")

# Displaying data in pairplot
alpha  <- Fecowid$alpha
status <- Fecowid$status
#model  <- as.numeric(Fecowid$model)
model  <- Fecowid$model
#
# Rename models
#
model <- ifelse(model=="                       ",NA,model)
model <- ifelse(model=="flnc_band              ","band",model)
model <- ifelse(model=="flnc_comp              ","comp",model)
model <- ifelse(model=="flnc_plaw              ","plaw",model)
model <- ifelse(model=="flnc_sbpl              ","sbpl",model)
#
# Computing and plotting logarithmic variables
#
lFecowid <- data.frame(log10(Fecowid[,2:4]),alpha,status,model )
names(lFecowid)[3] <- "peak flux"


#library(GGally)
#p = ggpairs(iris)
#p[2,1] = p[2,1] + scale_y_continuous(labels = scales::percent_format())

ggpairs(lFecowid,legend=1,columns=1:4, ggplot2::aes(colour=status), 
        lower = list(continuous = wrap("points", size=0.3,shape=19,alpha=0.5)), 
        diag = list(continuous = wrap("densityDiag", alpha=0.5)), 
        upper = list(continuous = wrap("cor", size = 2.5))) +
  theme_bw() + removeGrid() + 
  ggtitle('Distribution of GBM "couples" and "widows"') + 
  theme(legend.position="bottom",legend.title = element_blank())+ 
  theme(plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

outfile <- paste(drive,dir,"figs/fig07-dist-gbm.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)



# Testing differences between Swift "couples" and "widows"
FT90coup <- Fecowid$T90[status=="couples"]
FT90wid  <- Fecowid$T90[status=="widows"]
ks.test(FT90coup,FT90wid)

Fflucoup <- Fecowid$fluence[status=="couples"]
Ffluwid  <- Fecowid$fluence[status=="widows"]
ks.test(Fflucoup,Ffluwid)

FP1coup <- Fecowid$`peak flux`[status=="couples"]
FP1wid  <- Fecowid$`peak flux`[status=="widows"]
ks.test(FP1coup,FP1wid)

alphacoup <- Fecowid$alpha[status=="couples"]
alphawid  <- Fecowid$alpha[status=="widows"]
ks.test(alphacoup,alphawid)

# Computing cross correlations between status and Sp type
tb <- table(status,model)
tb

summary(tb)

# Computing cross correlations between status and Sp type

# Preparing Data for barplot

colnames <- rep(names(tb[1,]),dim(tb)[1])
value <- as.numeric(tb)
mstatus <- rep(row.names(tb),dim(tb)[2])
Spmod   <- colnames[c(1,5,2,6,3,7,4,8)]
sd <- data.frame(value, mstatus,Spmod)

# Displaying barplot

ggplot(sd, aes(fill=mstatus, y=value, x=Spmod)) + 
  geom_bar(position="dodge", stat="identity") + ggtitle("Distribution of Spectral types") + theme_bw()+ theme(plot.title = element_text(size=12,hjust = 0.5))

#ggsave("F:/Research/Fermi_Swift_comp/newfigs/12.pdf",width = 6, height = 4)




# Outliers' diagnostics of GBM data
bpf1 <- ggplot(lFecowid,aes(T90)) + geom_boxplot(fill="darkred",alpha=0.5,colour="black") + theme_bw() +   coord_flip() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()) + theme(axis.title.y = element_text(size=fontsize*1.2),axis.text.y = element_text(size=fontsize*0.8))
bpf2 <- ggplot(lFecowid,aes(fluence)) + geom_boxplot(fill="darkred",alpha=0.5,colour="black") + theme_bw() +   coord_flip() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+ theme(axis.title.y = element_text(size=fontsize*1.2),axis.text.y = element_text(size=fontsize*0.8))
bpf3 <- ggplot(lFecowid,aes(`peak flux`)) + geom_boxplot(fill="darkred",alpha=0.5,colour="black") + theme_bw() +   coord_flip() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+ theme(axis.title.y = element_text(size=fontsize*1.2),axis.text.y = element_text(size=fontsize*0.8))
bpf4 <- ggplot(lFecowid,aes(alpha)) +geom_boxplot(fill="darkred",alpha=0.5,colour="black") + theme_bw() + coord_flip() + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+ theme(axis.title.y = element_text(size=fontsize*1.2),axis.text.y = element_text(size=fontsize*0.8))

plot_grid(bpf1,bpf2,bpf3,bpf4,ncol=4,nrow=1)

outfile <- paste(drive,dir,"figs/fig03b-outlier-gbm.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)

#ggsave("E:/Research/Fermi_Swift_comp/________rev6/fig05b-outlier-gbm.pdf",width = 6, height = 4)

#
fep1 <- boxplot.stats(lFecowid$T90) 
fep2 <- boxplot.stats(lFecowid$fluence)
fep3 <- boxplot.stats(lFecowid$`peak flux`)
fep4 <- boxplot.stats(lFecowid$alpha)
#
inlierf1 <- lFecowid$T90 > fep1$stats[1] & lFecowid$T90 < fep1$stats[5]
inlierf2 <- lFecowid$fluence > fep2$stats[1] & lFecowid$fluence < fep2$stats[5]
inlierf3 <- lFecowid$`peak flux`  > fep3$stats[1] & lFecowid$`peak flux`  < fep3$stats[5]
inlierf4 <- lFecowid$alpha > fep4$stats[1] & lFecowid$alpha < fep4$stats[5]
inlierf   <- inlierf1 & inlierf2 & inlierf3 & inlierf4
#
# Identify outliers as NA
# 
inlierf <- ifelse(inlierf,inlierf,NA)
#
# Summary of boxplot's outlier diagnostics 
# 
summary(data.frame(inlierf1,inlierf2,inlierf3,inlierf4, inlierf))
#


#
# Performing LDA
#
nlFecowid <- na.omit(lFecowid[inlierf,])
Felda <- lda(nlFecowid[,1:4],nlFecowid$status)
Felda
#
# Computing best discriminating variable
#
FLD1 <- Felda$scaling
FLD1 <- as.matrix(nlFecowid[,1:4]) %*% FLD1
#
# Testing the discriminant power of FLD1, using ANOVA
#
df <-data.frame(nlFecowid,FLD1)
dfcor <- df[,c(5,7)]
aov.res <- aov(FLD1 ~ status,data=dfcor)
summary(aov.res)
#
# Displaying discriminant power of FLD1
#
pFLD1 <- ggplot(df, aes(FLD1, fill = status))  + geom_density(alpha = 0.5) + 
  ggtitle("Greatest couples-widows difference") +theme_bw()+ 
  theme(plot.title = element_text(size=11,hjust = 0.0),legend.title = element_blank()) + 
  labs(x="LD1")
#
# Correlating observed variables with LD1
#
corFLD1 <- cor(nlFecowid[,1:4],FLD1,method="pearson")
cordf <- data.frame(value=as.numeric(corFLD1),variables=row.names(corFLD1))
corFLD1
#
# Displaying correlations p=0.05 corr =0.075
# 
pFcor <- ggplot(cordf, aes(x=variables,y=corFLD1)) + 
  labs(y="correlation")  + geom_bar(stat="identity",fill="#00bfc4",alpha=0.5,colour="black") + 
  ggtitle("Correlation of LD1 with variables") +theme_bw() + ylim(-1,1) + 
  geom_hline(yintercept=c(-0.075,0.075),linetype="dashed",
             color = "red", size=1) + 
  geom_text(x= 2.5, y= 0.5, label="- -: p-value = 0.05", color="red",size=4.5)+ 
  theme(plot.title = element_text(size=11,hjust = 0.5))
#
plot_grid(pFLD1,pFcor,nrow=1,ncol=2)

outfile <- paste(drive,dir,"figs/fig10-sep-gbm.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)



# Displaying dependence between observed and LD1 variables

pFT90 <- ggplot(df,aes(x=FLD1,y= T90,color=status)) + geom_point(size=0.7,alpha=0.5) + theme_bw() + ggtitle("LDA in T90")+ theme(plot.title = element_text(size=12,hjust = 0.5),legend.title = element_blank()) + labs(x="LD1")
pFflu <- ggplot(df,aes(x=FLD1,y= fluence,color=status)) + geom_point(size=0.7,alpha=0.5) + theme_bw() + ggtitle("LDA in fluence")+ theme(plot.title = element_text(size=12,hjust = 0.5),legend.title = element_blank()) + labs(x="LD1")
pFP1 <- ggplot(df,aes(x=FLD1,y= peak.flux,color=status)) + geom_point(size=0.7,alpha=0.5) + theme_bw() + ggtitle("LDA in peak flux")+ theme(plot.title = element_text(size=12,hjust = 0.5),legend.title = element_blank()) + labs(x="LD1",y="peak flux")
palpha <- ggplot(df,aes(x=FLD1,y= alpha,color=status)) + geom_point(size=0.7,alpha=0.5) + theme_bw() + ggtitle("LDA in alpha")+ theme(plot.title = element_text(size=12,hjust = 0.5),legend.title = element_blank()) + labs(x="LD1")

plot_grid(pFT90,pFflu,pFP1,palpha,nrow=2,ncol=2)

outfile <- paste(drive,dir,"figs/fig11-lda-gbm.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)



# Displaying frequency distributions of variables

pFT90 <-   ggplot(df,aes(x=T90,fill=status)) +   geom_density(size=0.7,alpha=0.4) + 
  theme_bw() + ggtitle("FT90 distribution")+ 
  theme(plot.title = element_text(size=12,hjust = 0.5))
#
pFflu <- ggplot(df,aes(x= fluence,fill=status)) + geom_density(size=0.7,alpha=0.4) + theme_bw() + ggtitle("Fflu distribution")+ theme(plot.title = element_text(size=12,hjust = 0.5))
#
pFP1 <- ggplot(df,aes(x= peak.flux,fill=status)) + geom_density(size=0.7,alpha=0.4) + theme_bw() + ggtitle("FP1 distribution")+ theme(plot.title = element_text(size=12,hjust = 0.5))
#
palpha <- ggplot(df,aes(x= alpha,fill=status)) + geom_density(size=0.7,alpha=0.4) + theme_bw() + ggtitle("LDA in alpha")+ theme(plot.title = element_text(size=12,hjust = 0.5))
#
plot_grid(pFT90,pFflu,pFP1,palpha,nrow=2,ncol=2)

#ggsave("F:/Research/Fermi_Swift_comp/newfigs/16.pdf",width = 6, height = 4)



##### CC #####


library(CCA)
library(CCP)


SwFecoup <- SwFe_couples[,c(5,6,8,10,17+4,20+4,22+4,33+4)]

names(SwFecoup) <- c("Swift T90","Swift fluence","Swift peak flux","phi","Fermi T90","Fermi fluence","Fermi peak flux","alph")

nSwFecoup <- na.omit(SwFecoup)
nSwFecoup <- SwFecoup 

inSw <- data.frame(log10(nSwFecoup[,c(1:3)]),phi=nSwFecoup[,4])

inFe <- data.frame(log10(nSwFecoup[,c(5:7)]),alpha=nSwFecoup[,8])

# Perform canonical correlations
#
cc.res <- cc(inSw,inFe)
#
# Retrieve canonical scores for BAT (xscores) and GBM (yscores)
#
U <- cc.res$scores$xscores
V <- cc.res$scores$yscores
#
# Estimate signficance with Wilks Lambda
#
rho <- cc.res$cor
N <- dim(inSw)[1]
p <- dim(inSw)[2]
q <- dim(inFe)[2] 
#
p.asym(rho, N, p, q, tstat = "Wilks")


# Display results

# Plot U-BAT, pairwise

UBAT <- data.frame(U,inSw)
names(UBAT)[1:8] <- c("U1","U2","U3","U4","T90","fluence","peak flux","phi")
ggpairs(UBAT, columns = c(1:8),
        lower = list(continuous = wrap(ggally_density, size =0.1, color = "blue")),
        upper = list(continuous=wrap(ggally_cor,size=3,color="black")),
        diag = list(continuous = wrap(ggally_barDiag,bins=20,fill="lightblue",color = "blue")))  + theme_bw() + removeGrid() + ggtitle("BAT & U data (outliers excluded)")+ theme(plot.title = element_text(size=12,hjust = 0.5),axis.text = element_text(size=7))

outfile <- paste(drive,dir,"figs/fig12-bat-u.pdf",sep="")
ggsave(outfile,width = 9, height = 6,dpi=300)


# Plot V-Bat, pairwise 

VBAT <- data.frame(V,inSw)
names(VBAT)[1:8] <- c("V1","V2","V3","V4","T90","fluence","peak flux","phi")
ggpairs(VBAT, columns = c(1:8),
        lower = list(continuous = wrap(ggally_density, size =0.1, color = "blue")),
        upper = list(continuous=wrap(ggally_cor,size=3,color="black")),
        diag = list(continuous = wrap(ggally_barDiag,bins=20,fill="lightblue",color = "blue")))  + theme_bw() + removeGrid() + ggtitle("BAT & V data (outliers excluded)")+ theme(plot.title = element_text(size=12,hjust = 0.5),axis.text = element_text(size=7))

outfile <- paste(drive,dir,"figs/fig13-bat-v.pdf",sep="")
ggsave(outfile,width = 9, height = 6,dpi=300)



# Plot U-GBM pairwise 

UGBM <- data.frame(U,inFe)
names(UGBM)[1:8] <- c("U1","U2","U3","U4","T90","fluence","peak flux","alpha")
ggpairs(UGBM, columns = c(1:8),
        lower = list(continuous = wrap(ggally_density, size =0.1, color = "darkred")),
        upper = list(continuous=wrap(ggally_cor,size=3,color="black")),
        diag = list(continuous = wrap(ggally_barDiag,bins=20,fill="pink",color = "darkred")))  + theme_bw() + removeGrid() + ggtitle("GBM & U data (outliers excluded)")+ theme(plot.title = element_text(size=12,hjust = 0.5),axis.text = element_text(size=7))

outfile <- paste(drive,dir,"figs/fig14-gbm-u.pdf",sep="")
ggsave(outfile,width = 9, height = 6,dpi=300)



# Plot V-GBM, pairwise 

VGBM <- data.frame(V,inFe)
names(VGBM)[1:8] <- c("V1","V2","V3","V4","T90","fluence","peak flux","alpha")
ggpairs(VGBM, columns = c(1:8),
        lower = list(continuous = wrap(ggally_density, size =0.1, color = "darkred")),
        upper = list(continuous=wrap(ggally_cor,size=3,color="black")),
        diag = list(continuous = wrap(ggally_barDiag,bins=20,fill="pink",color = "darkred")))  + theme_bw() + removeGrid() + ggtitle("GBM & V data (outliers excluded)")+ theme(plot.title = element_text(size=12,hjust = 0.5),axis.text = element_text(size=7))

outfile <- paste(drive,dir,"figs/fig15-gbm-v.pdf",sep="")
ggsave(outfile,width = 9, height = 6,dpi=300)









# Color display of correlations

library(corrplot)

# Compute structure martices for u and v values

UV  <- data.frame(U,V)
cxuv <- cor(inSw,UV,use="pairwise.complete")
cyuv <- cor(inFe,UV,use="pairwise.complete")

# Cross correlations of U and V with BAT and GBM data

colnames(cxuv) <- c("U1","U2","U3","U4","V1","V2","V3","V4")
rownames(cxuv) <- c("Sw: T90","Sw: fluence","Sw: peak flux","Sw: phi")
colnames(cyuv) <- c("U1","U2","U3","U4","V1","V2","V3","V4")
rownames(cyuv) <- c("Fe: T90","Fe: fluence","Fe: peak flux","Fe: alpha")
cxyuv <- round(rbind(cxuv,cyuv),3)
cxyuv



# Estimate significance

psig <- 0.003 # 3sigma

colx <- dim(inSw)[2]
coly <- dim(inFe)[2]
coluv <- dim(UV)[2]
pmatx <- matrix(rep(0,colx*coluv),colx,coluv)
pmaty <- matrix(rep(0,coly*coluv),coly,coluv)

for(i in 1:colx){
  for(j in 1:coluv){pmatx[i,j] <- cor.test(inSw[,i],UV[,j])$p.value
  cxuv[i,j] <- ifelse(pmatx[i,j]>psig,0,cxuv[i,j])}
}

for(i in 1:coly){
  for(j in 1:coluv){pmaty[i,j] <- cor.test(inFe[,i],UV[,j])$p.value
  cyuv[i,j] <- ifelse(pmaty[i,j]>psig,0,cyuv[i,j])}
}

# Structure matrix resulted

cxyuv <- round(rbind(cxuv,cyuv),3)
cxyuv


# Colored display of the structure matrix

outfile <- paste(drive,dir,"figs/fig16-corrplot.pdf",sep="")
#ggsave(outfile,width = 6.5, dpi=300)


pdf(file = outfile,width = 6.5)
corrplot(0-cxyuv[,c(1,5,2,6,3,7)],is.cor=T,tl.cex=1,cl.align.text = "c",cl.pos="b",cl.cex=1,cl.lim=c(-1,1),cl.length=5)
dev.off()






##### Compare Swift and Fermi T90 distributions #####

T90 <- c(inSw[,1],inFe[,1])
satellite <- c(rep("Swift",length(inSw[,1])),rep("Fermi",length(inFe[,1])))
dfT90 <- data.frame(T90,satellite)

# Display jointly Swift and Fermi T90 distributions

ggplot(dfT90,aes(x = T90,color = satellite, fill = satellite)) +  
  geom_density(size=0.7,alpha=0.4) + theme_bw() + 
  ggtitle("Fermi-Swift T90 distributions")+ 
  theme(plot.title = element_text(size=12,hjust = 0.5),legend.title = element_blank())+
  labs(x="lg(T90)") +
  scale_fill_manual(values=c("darkred","darkblue")) + 
  scale_color_manual(values=c("darkred","darkblue"))


outfile <- paste(drive,dir,"figs/fig18-t90-dist.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)






# Perform Mclust on Fermi and Swift T90 data

library(mclust)

Fclust <- Mclust((inFe[,1]))
inSw2<-inSw[!is.na(inSw$Swift.T90), ]
Sclust <- Mclust((inSw2[,1]))


clust_nums <- c(1,2,3,4,5,6,7,8,9)

clusters <- cbind(clust_nums,Fclust$BIC[,1],Fclust$BIC[,2],Sclust$BIC[,1],Sclust$BIC[,2])
colnames(clusters) <- c("No. of clusters","Fermi E","Fermi V","Swift E","Swift V")
rownames(clusters) <- clust_nums
clusters = as.data.frame(clusters)
#clusters$`Swift V`[7] <- (clusters$`Swift V`[6]+clusters$`Swift V`[8])/2

library(reshape)
clusters_to_plot <- melt(clusters,id.vars = "No. of clusters")
#ggplot(clusters_to_plot,aes(x="No. of clusters",y=value,col=variable)) + geom_line()
colnames(clusters_to_plot)[2] <- "Model"
colnames(clusters_to_plot)[3] <- "BIC"
ggplot(clusters_to_plot,aes(x=`No. of clusters`, y=BIC,col= Model)) + geom_line(linetype = "dashed") + theme_bw() + geom_point() + scale_color_manual(values=c("red3","#F8766D", "blue3","#00bfc4")) + 
  scale_x_continuous(name="Number of components", breaks=seq(1, 9, by = 1))

outfile <- paste(drive,dir,"figs/fig17-bic.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)




##### log-log scale of observed variables #####



#SwFe_couples

swp1 <- boxplot.stats(SwFe_couples$T90) 
swp2 <- boxplot.stats(SwFe_couples$flu)
swp3 <- boxplot.stats(SwFe_couples$P1)
swp4 <- boxplot.stats(SwFe_couples$Phi)

inlierg1 <- SwFe_couples$T90 > swp1$stats[1]-0.5 & SwFe_couples$T90 < swp1$stats[5]
inlierg2 <- SwFe_couples$flu > swp2$stats[1] & SwFe_couples$flu < swp2$stats[5]
inlierg3 <- SwFe_couples$P1  > swp3$stats[1] & SwFe_couples$P1  < swp3$stats[5]
inlierg4 <- SwFe_couples$Phi > swp4$stats[1] & SwFe_couples$Phi < swp4$stats[5]
inlier_s   <- inlierg1 & inlierg2 & inlierg3 & inlierg4

# Identify outliers as NA
inlier_s <- ifelse(inlier_s,inlier_s,NA)




fep1 <- boxplot.stats(SwFe_couples$t90) 
fep2 <- boxplot.stats(SwFe_couples$fluence)
fep3 <- boxplot.stats(SwFe_couples$flux_1024)
fep4 <- boxplot.stats(SwFe_couples$flnc_band_alpha)
#
inlierf1 <- SwFe_couples$t90 > fep1$stats[1] & SwFe_couples$t90 < fep1$stats[5]
inlierf2 <- SwFe_couples$fluence > fep2$stats[1] & SwFe_couples$fluence < fep2$stats[5]
inlierf3 <- SwFe_couples$flux_1024  > fep3$stats[1] & SwFe_couples$flux_1024  < fep3$stats[5]
inlierf4 <- SwFe_couples$flnc_band_alpha > fep4$stats[1] & SwFe_couples$flnc_band_alpha < fep4$stats[5]
inlier_f   <- inlierf1 & inlierf2 & inlierf3 & inlierf4
#
# Identify outliers as NA
# 
inlier_f <- ifelse(inlier_f,inlier_f,FALSE)
#
# Summary of boxplot's outlier diagnostics 
# 
#summary(data.frame(inlierf1,inlierf2,inlierf3,inlierf4, inlier_f))

inlier_all <- inlierf1 & inlierf2 & inlierf3 & inlierf4 & inlierg1 & inlierg2 & inlierg3 & inlierg4
inlier_all <- ifelse(inlier_all,inlier_all,FALSE)

SwFe_couples_in <- na.omit(SwFe_couples[inlier_all,])
SwFe_couples_out <- na.omit(SwFe_couples[!inlier_all,])
SwFe_couples_all <- na.omit(SwFe_couples)



ST90 <- SwFe_couples_in$T90
FT90 <- SwFe_couples_in$t90
#
Sflu <- SwFe_couples_in$flu
Fflu <- SwFe_couples_in$fluence
#
Speak <- SwFe_couples_in$P1
Fpeak <- SwFe_couples_in$flux_1024
#
# Plotting jointly detected variables
#
lmpT <- lm(log10(SwFe_couples_in$t90)~log10(SwFe_couples_in$T90))
pT <- ggplot(SwFe_couples_in,aes(log10(T90),log10(t90)))  +
  stat_smooth(method = "lm", col = "red",alpha=0.5,level=0.997,fullrange=T,size=0.5) +
  stat_smooth(data=SwFe_couples_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange=T,size=0.5) + 
  geom_point(shape=19, color="#00bfc4",size=1) + 
  geom_point(data=SwFe_couples_out,shape=19, color="magenta",size=1) + 
  theme_bw() + ggtitle("Swift-Fermi T90 durations") + 
  geom_abline(intercept = 0, slope = 1,alpha=0.5, color="blue",  linetype="dashed", size=1) + 
  annotate(geom = 'text', label = '- - -: y = x', x = -Inf, y = Inf, hjust = -0.1, vjust = 2, color="blue", size=4) + 
  annotate(geom = 'text', label=paste("slope = ",round(lmpT$coefficients[2],2),'\u00B1',round(coef(summary(lmpT))[2,"Std. Error"],2),sep=""),adj=1.1,x = Inf, y = -Inf, vjust = -1.7, color="red", size=4) +
  annotate(geom = 'text', label=paste("intercept = ",round(lmpT$coefficients[1],2),'\u00B1',round(coef(summary(lmpT))[1,"Std. Error"],2),sep=""),adj=1.1,x = Inf, y = -Inf, vjust = -0.5, color="red", size=4) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab(expression(paste("Swift lg(",T[90],")",sep=""))) + 
  ylab(expression(paste("Fermi lg(",T[90],")",sep="")))

lmpF <- lm(log10(SwFe_couples_in$fluence)~log10(SwFe_couples_in$flu))
pF <- ggplot(SwFe_couples_in,aes(log10(flu),log10(fluence))) +
  stat_smooth(method = "lm", col = "red",alpha=0.5,level=0.997,fullrange=T,size=0.5) +
  stat_smooth(data=SwFe_couples_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange=T,size=0.5) + 
  geom_point(shape=19, color="#00bfc4",size=1) + 
  geom_point(data=SwFe_couples_out,shape=19, color="magenta",size=1) + 
  theme_bw() + ggtitle("Swift-Fermi fluences") + 
  geom_abline(intercept = 0, slope = 1, color="blue",alpha=0.5,  linetype="dashed", size=1) + 
  annotate(geom = 'text', label='- - -: y = x', x = -Inf, y = Inf, 
           hjust = -0.2, vjust = 1.75, color="blue", size=4) + 
  annotate(geom = 'text', label=paste("slope = ",round(lmpF$coefficients[2],2),'\u00B1',round(coef(summary(lmpF))[2,"Std. Error"],2),sep=""),adj=1.1,x = Inf, y = -Inf, vjust = -1.7, color="red", size=4) +
  annotate(geom = 'text', label=paste("intercept = ",round(lmpF$coefficients[1],2),'\u00B1',round(coef(summary(lmpF))[1,"Std. Error"],2),sep=""),adj=1.1,x = Inf, y = -Inf, vjust = -0.5, color="red", size=4) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-8.0, -3) + ylim(-8.0, -3) + 
  xlab(expression(paste("Swift lg(fluence)",sep=""))) + 
  ylab(expression(paste("Fermi lg(fluence)",sep=""))) 

lmpP <- lm(log10(SwFe_couples_in$flux_1024)~log10(SwFe_couples_in$P1))
pP <- ggplot(SwFe_couples_in,aes(log10(P1),log10(flux_1024))) +
  stat_smooth(method = "lm", col = "red",alpha=0.5,level=0.997,fullrange=T,size=0.5) +
  stat_smooth(data=SwFe_couples_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange=T,size=0.5) + 
  geom_point(shape=19, color="#00bfc4",size=1) +  
  geom_point(data=SwFe_couples_out,shape=19, color="magenta",size=1) + 
  theme_bw() + ggtitle("Swift-Fermi 1024 ms peak flux") + 
  geom_abline(intercept=0,slope=1,color="blue",linetype="dashed",size=1,alpha=0.5) + 
  annotate(geom = 'text', label=paste("slope = ",round(lmpP$coefficients[2],2),'\u00B1',round(coef(summary(lmpP))[2,"Std. Error"],2),sep=""),adj=1.1,x = Inf, y = -Inf, vjust = -1.7, color="red", size=4) +
  annotate(geom = 'text', label=paste("intercept = ",round(lmpP$coefficients[1],2),'\u00B1',round(coef(summary(lmpP))[1,"Std. Error"],2),sep=""),adj=1.1,x = Inf, y = -Inf, vjust = -0.5, color="red", size=4) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  annotate(geom = 'text', label = '- - -: y = x',x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, color="blue", size=4)  + 
  xlim(-1, 2.5) + ylim(-1, 2.5) + 
  xlab("Swift lg(peak flux)") + 
  ylab("Fermi lg(peak flux)") 

plot_grid(pT,pF,pP,ncol=1,nrow=3)

outfile <- paste(drive,dir,"figs/fig02-log-log.pdf",sep="")
ggsave(outfile,width = 6, height = 6,dpi=300)



###### SWIFT errors #################


# Compute logarithmic variables

lT90c <- log10(swg$T90gcn)
lT90 <-  log10(swg$T90)
lflu <-  log10(swg$flu)
lP1  <-  log10(swg$P1)

lT90e <- log10(swg$T90err)
lfle  <- log10(swg$fle)
lp1e  <- log10(swg$p1e)

lswg <- data.frame(Name=paste("GRB",swg$Name,sep=""),lT90c,lT90,lflu,lP1,Phi=swg$Phi,lT90e,lfle,lp1e,Phie=swg$Phie)

swp1 <- boxplot.stats(lswg$lT90) 
swp2 <- boxplot.stats(lswg$lflu)
swp3 <- boxplot.stats(lswg$lP1)
swp4 <- boxplot.stats(lswg$Phi)

inlierg1 <- lswg$lT90 > swp1$stats[1]-0.5 & lswg$lT90 < swp1$stats[5]
inlierg2 <- lswg$lflu > swp2$stats[1] & lswg$lflu < swp2$stats[5]
inlierg3 <- lswg$lP1  > swp3$stats[1] & lswg$lP1  < swp3$stats[5]
inlierg4 <- lswg$Phi > swp4$stats[1] & lswg$Phi < swp4$stats[5]
inlier   <- inlierg1 & inlierg2 & inlierg3 & inlierg4

# Identify outliers as NA
#inlier <- ifelse(inlier,inlier,FALSE)
#
# Summary of boxplot's outlier diagnostics 
# 
#summary(data.frame(inlierg1,inlierg2,inlierg3,inlierg4, inlier))
#
# Removing NAs from the data
#
nlswg <- na.omit(lswg[inlier,])# & (lfle > 0)
nlswg_out <- na.omit(lswg[!inlier,])# & (lfle > 0)
nlswg_all <- na.omit(lswg[,])#(lfle > 0)
#nlswg <- na.omit(lswg[,])#(lfle > 0)

pT90 <- ggplot(nlswg,aes(lT90,lT90e)) +
  geom_point(data=nlswg_out,shape=19, color="darkred",size=0.3)+ 
  geom_point(data=nlswg,shape=19, color="blue",size=0.3) + theme_bw() + 
  ggtitle(expression(paste("Swift ",T[90]," vs error",sep=""))) + 
  theme(plot.title = element_text(hjust = 0.5,size = fontsize*0.8),legend.text = element_text(size = fontsize*0.4), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        strip.text = element_text(size = fontsize*0.8, colour='black')) + 
  #   geom_abline(intercept = -0.8, slope = 1, color="red",  linetype="dashed", size=1) +
  xlab(expression(paste("lg(",T[90],")",sep=""))) + 
  ylab(expression(paste("lg(",T[90]," error)",sep=""))) +
  stat_smooth(data=nlswg_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange = TRUE) +
  stat_smooth(method = "lm", col = "red",alpha=0.5,level=0.997,fullrange = TRUE) 
#
pflu <- ggplot(nlswg,aes(lflu,lfle)) +
  geom_point(data=nlswg_out,shape=19, color="darkred",size=0.3)+ 
  geom_point(data=nlswg,shape=19, color="blue",size=0.3) + theme_bw() + 
  ggtitle("Swift fluence vs error") + 
  theme(plot.title = element_text(hjust = 0.5,size = fontsize*0.8),legend.text = element_text(size = fontsize*0.4), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        strip.text = element_text(size = fontsize*0.8, colour='black')) + 
  #   geom_abline(intercept = -0.8, slope = 0.8, color="red",  linetype="dashed", size=1) +
  xlab("lg(fluence)") + ylab("lg(fluence error)") +
  stat_smooth(data=nlswg_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange = TRUE) +
  stat_smooth(method = "lm", col = "red",alpha=0.5,level=0.997,fullrange = TRUE)
#
pP1 <- ggplot(nlswg,aes(lP1,lp1e)) +
  geom_point(data=nlswg_out,shape=19, color="darkred",size=0.3)+ 
  geom_point(data=nlswg,shape=19, color="blue",size=0.3) + theme_bw() + 
  ggtitle("Swift peak flux vs error") + 
  theme(plot.title = element_text(hjust = 0.5,size = fontsize*0.8),legend.text = element_text(size = fontsize*0.4), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        strip.text = element_text(size = fontsize*0.8, colour='black')) + 
  #   geom_abline(intercept = -0.7, slope = 0.5, color="red",  linetype="dashed", size=1)+
  xlab("lg(peak flux)") + ylab("lg(peak flux error)") +
  stat_smooth(data=nlswg_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange = TRUE) +
  stat_smooth(method = "lm", col = "red",alpha=0.5,level=0.997,fullrange = TRUE)
#
pPhi <- ggplot(nlswg,aes(Phi,log10(Phie))) +
  geom_point(data=nlswg_out,shape=19, color="darkred",size=0.3)+ 
  geom_point(data=nlswg,shape=19, color="blue",size=0.3) + theme_bw() + 
  ggtitle("Swift photon index vs error") + 
  theme(plot.title = element_text(hjust = 0.5,size = fontsize*0.8),legend.text = element_text(size = fontsize*0.4), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        strip.text = element_text(size = fontsize*0.8, colour='black')) + 
  #   geom_hline(yintercept = 0.2, color="red",  linetype="dashed", size=1)+ xlab("Phind") +
  xlab("photon index") + ylab("photon index error")  +
  stat_smooth(data=nlswg_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange = TRUE) +
  stat_smooth(method = "lm", col = "red",alpha=0.5,level=0.997,fullrange = TRUE)
# 
plot_grid(pT90 ,pflu,pP1,pPhi,ncol=2,nrow=2) 

outfile <- paste(drive,dir,"figs/fig04-swift-error.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)



###### Fermi error ####################

# Compute logarithmic variables

lFT90 <- log10(gbm$t90)
lFflu <- log10(gbm$fluence)
lFP1  <-log10(gbm$flux_1024)
alpha <- gbm$flnc_band_alpha

lFT90e <- log10(gbm$t90_error)
lFflue <- log10(gbm$fluence_error)
lFP1e  <- log10(gbm$flux_1024_error)
alphae <- gbm$flnc_band_alpha_pos_err

lgbm <- data.frame(Fname=gbm$Fname,lFT90,lFflu,lFP1,alpha,lFT90e,lFflue,lFP1e,alphae)

# display data in matrix plot

#pairs(lgbm,pch=19,col="darkred",main="Fermi GBM GRBs")
#
# Outliers' diagnostics of GBM data
#
# Using boxplots
#
#bpf1 <- ggplot(lgbm,aes(lFT90)) + geom_boxplot(col="black",fill="darkred") + theme_bw() +   coord_flip()
#bpf2 <- ggplot(lgbm,aes(lFflu)) + geom_boxplot(col="black",fill="darkred") + theme_bw() +   coord_flip()
#bpf3 <- ggplot(lgbm,aes(lFP1)) + geom_boxplot(col="black",fill="darkred") + theme_bw() +   coord_flip()
#bpf4 <- ggplot(lgbm,aes(alpha)) +geom_boxplot(col="black",fill="darkred") + theme_bw() +   coord_flip() 
#
#plot_grid(bpf1,bpf2,bpf3,bpf4,ncol=4,nrow=1)
#
#
fep1 <- boxplot.stats(lgbm$lFT90) 
fep2 <- boxplot.stats(lgbm$lFflu)
fep3 <- boxplot.stats(lgbm$lFP1)
fep4 <- boxplot.stats(lgbm$alpha)
#
inlierf1 <- lgbm$lFT90 > fep1$stats[1] & lgbm$lFT90 < fep1$stats[5]
inlierf2 <- lgbm$lFflu > fep2$stats[1] & lgbm$lFflu < fep2$stats[5]
inlierf3 <- lgbm$lFP1  > fep3$stats[1] & lgbm$lFP1  < fep3$stats[5]
inlierf4 <- lgbm$alpha > fep4$stats[1] & lgbm$alpha < fep4$stats[5]
inlierf   <- inlierf1 & inlierf2 & inlierf3 & inlierf4
#
# Identify outliers as NA
# 
#inlierf <- ifelse(inlierf,inlierf,NA)
#
# Summary of boxplot's outlier diagnostics 
# 
#summary(data.frame(inlierf1,inlierf2,inlierf3,inlierf4, inlierf))
#
# Removing NAs from the data
#
nlgbm <- na.omit(lgbm[inlierf,])
nlgbm_out <- na.omit(lgbm[!inlier,])
nlgbm_all <- na.omit(lgbm)


pFT90 <- ggplot(nlgbm,aes(lFT90,lFT90e)) +
  geom_point(data=nlgbm_out,shape=19,size=0.3,  color="blue")+ 
  geom_point(shape=19, color="darkred",size=0.3) + theme_bw() + 
  ggtitle(expression(paste("Fermi ",T[90]," vs error",sep=""))) + 
  theme(plot.title = element_text(hjust = 0.5,size = fontsize*0.8),legend.text = element_text(size = fontsize*0.4), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        strip.text = element_text(size = fontsize*0.8, colour='black')) + 
  #   geom_abline(intercept = -0.5, slope = 0.6, color="black",  linetype="dashed", size=1) +
  xlab(expression(paste("lg(",T[90],")",sep=""))) + ylab(expression(paste("lg(",T[90]," error)",sep=""))) +
  stat_smooth(data=nlgbm_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange = TRUE) +
  stat_smooth(method = "lm", col = "blue",alpha=0.5,level=0.997,fullrange = TRUE)

pFflu <- ggplot(nlgbm,aes(lFflu,lFflue)) +
  geom_point(data=nlgbm_out,shape=19, color="blue",size=0.3)+ 
  geom_point(shape=19, color="darkred",size=0.3) + theme_bw() + 
  ggtitle("Fermi fluence vs error") + 
  theme(plot.title = element_text(hjust = 0.5,size = fontsize*0.8),legend.text = element_text(size = fontsize*0.4), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        strip.text = element_text(size = fontsize*0.8, colour='black')) + 
  #   geom_abline(intercept = -4.9, slope = 0.45, color="black",  linetype="dashed", size=1) +
  xlab("lg(fluence)") + ylab("lg(fluence error)") +
  stat_smooth(data=nlgbm_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange = TRUE) +
  stat_smooth(method = "lm", col = "blue",alpha=0.5,level=0.997,fullrange = TRUE)
#
pFP1 <- ggplot(nlgbm,aes(lFP1,lFP1e)) +
  geom_point(data=nlgbm_out,shape=19, color="blue",size=0.3)+ 
  geom_point(shape=19, color="darkred",size=0.3) + theme_bw() + 
  ggtitle("Fermi peak flux vs error") + 
  theme(plot.title = element_text(hjust = 0.5,size = fontsize*0.8),legend.text = element_text(size = fontsize*0.4), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        strip.text = element_text(size = fontsize*0.8, colour='black')) + 
  #   geom_abline(intercept = -0.8, slope = 0.35, color="black",  linetype="dashed", size=1) +
  xlab("lg(peak flux)") + ylab("lg(peak flux error)") +
  stat_smooth(data=nlgbm_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange = TRUE) +
  stat_smooth(method = "lm", col = "blue",alpha=0.5,level=0.997,fullrange = TRUE)
#
pFalph <- ggplot(nlgbm,aes((alpha),log10(alphae))) +
  geom_point(data=nlgbm_out,shape=19, color="blue",size=0.3) + 
  geom_point(shape=19, color="darkred",size=0.3) + theme_bw() + 
  ggtitle("Fermi alpha vs error") + 
  theme(plot.title = element_text(hjust = 0.5,size = fontsize*0.8),legend.text = element_text(size = fontsize*0.4), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        strip.text = element_text(size = fontsize*0.8, colour='black')) + 
  #   geom_hline(yintercept = 0.2, color="black",  linetype="dashed", size=1) + 
  xlab("band alpha") + ylab("alpha error")  +
  stat_smooth(data=nlgbm_all,method = "lm", col = "green",alpha=0.5,level=0.997,fullrange = TRUE) +
  stat_smooth(method = "lm", col = "blue",alpha=0.5,level=0.997,fullrange = TRUE) + xlim(-2,1)

plot_grid(pFT90 ,pFflu,pFP1,pFalph,ncol=2,nrow=2) 
  
outfile <- paste(drive,dir,"figs/fig05-fermi-error.pdf",sep="")
ggsave(outfile,width = 6, height = 4,dpi=300)





