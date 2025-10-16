#
#
############################
### USING THIS SCRIPT... ###
############################
#
#
#
#####
# 1 #
#####
###
### You need an ANALYZE format file from ImageJ
### if the files are "2841_OD_R_1_0__.hdr" and "2841_OD_R_1_0__.img", 
### then the next line should read
#FILENAME="2841_OD_R_1_0__"
###
###

FILENAME="forregStack"



########################################################################################
########################################################################################
########################################################################################
####################
### REQUIREMENTS ###
####################
#
#####
# 1 #
#####
# You'll need to install the AnalyzefMRI package, if you haven't already.
# To see if you have it, try running the this line of code:
library(AnalyzeFMRI)
# If there was no error, you're good to go.  
#
# Otherwise,
# to install AnalyzefMRI package, go to Packages>>Install Packages>>
# and pick a location, then scroll around to find "AnalyzefMRI"
# (you'll need to be connected to the internet to do this).
#
#
#####
# 2 #
#####
# You'll need to install the RNiftyReg package, if you haven't already.
# To see if you have it, try running the this line of code:
library(RNiftyReg)
# If there was no error, you're good to go.  
#
# Otherwise,
# to install AnalyzefMRI package, go to Packages>>Install Packages>>
# and pick a location, then scroll around to find "AnalyzefMRI"
# (you'll need to be connected to the internet to do this).
#

DATA<-f.read.analyze.volume(paste(FILENAME,".hdr",sep=""))
if(length(dim(DATA))==4) DATA=DATA[,,,1]


TARGET=DATA[,,1]
## The mean will be used as the initial target
for(z in 2:dim(DATA)[3]) TARGET=TARGET+DATA[,,z]
TARGET=TARGET/dim(DATA)[3]

REG.DATA=DATA
REG.DATA[,,]<-NA
for(y in 1:dim(DATA)[3])
 {REG.DATA[,,y]=niftyreg.linear(DATA[,,y], TARGET, scope = "rigid", init = NULL, 
                sourceMask = NULL, targetMask = NULL, symmetric = FALSE, 
                nLevels = 3, maxIterations = 5, estimateOnly = FALSE)$image}

REG.DATA[which(is.nan(REG.DATA))]=NA

#REG.DATAdisplay=REG.DATA
#REG.DATAdisplay[which(is.na(REG.DATA))]=0
#f.write.analyze(REG.DATAdisplay[,,],paste(FILENAME,"_REG_step1",sep=""),size="float")
#rm(REG.DATAdisplay)


######################################################
##
## for the next step, we've hopefully
## taken care of most of the rotational aspects. 
##
## I've noticed some side-to-side jitter of the retina; we'll try to blunt this by just using the
## (STANDARD DEVIATION) signal within each "coulmn" of the image, which makes the bright spot for the vascular bundle "pop" visually
## (the bright cluster of blood vessels at the optic N head is a nice landmark)

## we'll only use the central 50%
BORDERS=round(c(dim(REG.DATA)[1]*1/4,dim(REG.DATA)[1]*3/4))
CENTRAL.HALF=REG.DATA[BORDERS[1]:BORDERS[2],,]

## the list of mean profiles we'll get...
PROFILE.LIST=matrix(,length(BORDERS[1]:BORDERS[2]),dim(REG.DATA)[3])

## and we need to use tapply
MEANx=function(x) mean(na.rm=T,x)
SDx=function(x) sd(na.rm=T,x)
CENTRAL.HALF.factor=CENTRAL.HALF[,,1]
for(x in 1:nrow(CENTRAL.HALF.factor)) CENTRAL.HALF.factor[x,]=x
CENTRAL.HALF.factor=as.factor(CENTRAL.HALF.factor)

#for(z in 1:dim(CENTRAL.HALF)[3]) PROFILE.LIST[,z]=as.numeric(tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,MEANx))
#for(z in 1:dim(CENTRAL.HALF)[3]) PROFILE.LIST[,z]=as.numeric(tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,SDx)/tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,MEANx))
for(z in 1:dim(CENTRAL.HALF)[3]) PROFILE.LIST[,z]=as.numeric(tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,SDx))

## as you can see, the data are noisey:
show.slice1=1
show.slice2=2
plot(PROFILE.LIST[,show.slice1],pch=19)
matpoints(PROFILE.LIST[,show.slice2],pch=21,col="red")

## anf we'll use the mean as a target:
MEAN.PROFILE=PROFILE.LIST[,1:2]
MEAN.PROFILE[,]<-NA
for(x in 1:nrow(MEAN.PROFILE)) MEAN.PROFILE[x,2]=mean(na.rm=T,PROFILE.LIST[x,])

plot(MEAN.PROFILE[,2],type="l")
matpoints(PROFILE.LIST[,show.slice1],pch=19,col="black")
matpoints(PROFILE.LIST[,show.slice2],pch=21,col="red")

rm(show.slice1,show.slice2)


## now, how to slide one past another to minimize variance?
## allow movements of +/- 30 pixels
## a negative number is to the left. 

EMPTY=matrix(,nrow(PROFILE.LIST),61)

LIST.OF.MOVEMENTS=matrix(,dim(REG.DATA)[3],2)
for(z in 1:ncol(PROFILE.LIST))
  {BUILD=EMPTY;
   STRIP=PROFILE.LIST[,z];
   BUILD[,31]=STRIP;
   for(x in seq(-30,-1,1)) BUILD[1:(length(STRIP)+x),(x+31)]=STRIP[((-1*x)+1):length(STRIP)];
   for(x in seq(1,30,1)) BUILD[(x+1):length(STRIP),(x+31)]=STRIP[1:(length(STRIP)-x)];
   ## now, compare correlation coefficients to the MEAN.PROFILE;
   BEST=c(0,0);
   for(x in seq(-30,30,1)) {A=cor.test(MEAN.PROFILE[,2],BUILD[,(x+31)])[4]$est; if(A>BEST[1]) {BEST[1]=A;BEST[2]=x}}
   LIST.OF.MOVEMENTS[z,2]=BEST[2]}

## now, apply those corrections:

REG.DATA2=REG.DATA
REG.DATA2[,,]<-NA
for(z in 1:nrow(LIST.OF.MOVEMENTS))
 {move=LIST.OF.MOVEMENTS[z,2];
  if(move==0) REG.DATA2[,,z]=REG.DATA[,,z];
  if(move<0) REG.DATA2[1:(dim(REG.DATA2)[1]+move),,z]=REG.DATA[((-1*move)+1):dim(REG.DATA2)[1],,z];
  if(move>0) REG.DATA2[(move+1):dim(REG.DATA2)[1],,z]=REG.DATA[1:(dim(REG.DATA2)[1]-move),,z]}

#image(REG.DATA2[,,1])

#REG.DATAdisplay=REG.DATA2
#REG.DATAdisplay[which(is.na(REG.DATA2))]=0
#f.write.analyze(REG.DATAdisplay[,,],paste(FILENAME,"_REG_step2",sep=""),size="float")
#rm(REG.DATAdisplay)


######################################################
##
## Now, we're going to break the retina into several "panels":
## One panel will be the central 1/4 of the image (6/16 to 10/16ths). The other 12 "panels" will be on each side of the center,
## representing 1/16th of the retina. 
##
## The panels will be shifted vertially to try and match the mean.
## Whatever the proposed vertical movement, we'll submit the whole thing to a smooth spline across the image
## so judge how we should move each "column" of pixels vertically.


TARGET=REG.DATA2[,,1]
to.average.across.z=matrix(,dim(REG.DATA2)[2],dim(REG.DATA2)[3])
for(x in 1:nrow(to.average.across.z)) to.average.across.z[x,]=x;
to.average.across.z=as.factor(to.average.across.z)
for(x in 1:nrow(TARGET)) TARGET[x,]=as.numeric(tapply(REG.DATA2[x,,],to.average.across.z,MEANx));

PANEL.BORDERS=matrix(,13,2)
#PANEL.BORDERS[1,1]=1
#PANEL.BORDERS[13,2]=dim(REG.DATA)[1]
for(x in 1:6)
 {PANEL.BORDERS[x,]=c(1+round(dim(REG.DATA)[1]*(x-1)/16),round(dim(REG.DATA)[1]*(x)/16))}
for(x in 8:13)
 {PANEL.BORDERS[x,]=c(1+round(dim(REG.DATA)[1]*(x+2)/16),round(dim(REG.DATA)[1]*(x+3)/16))}
PANEL.BORDERS[7,1]=PANEL.BORDERS[6,2]+1
PANEL.BORDERS[7,2]=PANEL.BORDERS[8,1]-1
PANEL.BORDERS=cbind(PANEL.BORDERS,PANEL.BORDERS[,1])
PANEL.BORDERS[,3]<-(PANEL.BORDERS[,1]+PANEL.BORDERS[,2])/2
colnames(PANEL.BORDERS)<-c("start","end","mean")


#
#

## we want a matrix of vertical moves for each panel, for each image:
# +1 will be "move up by 1 to match target", -1 will be to move down by 1. 
# check for any distances +/-15 pixels
PANEL.MOVES=matrix(,13,dim(REG.DATA2)[3])

#####################################################################################################
## first panel:
PANEL.NUMBER=1

PANEL=REG.DATA2[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],,]

PROFILE.LIST=matrix(,dim(REG.DATA)[2],dim(REG.DATA)[3])

PANEL.factor=PANEL[,,1]
for(x in 1:ncol(PANEL.factor)) PANEL.factor[,x]=x
PANEL.factor=as.factor(PANEL.factor)

for(z in 1:dim(PANEL)[3]) PROFILE.LIST[,z]=as.numeric(tapply(PANEL[,,z],PANEL.factor,MEANx))
## note:this is reading from the bottom to the top
PROFILE.TARGET=PROFILE.LIST[,1:2]
PROFILE.TARGET[,]<-NA
PROFILE.TARGET[,2]=as.numeric(tapply(TARGET[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],],PANEL.factor,MEANx))

EMPTY=matrix(,dim(REG.DATA2)[2],31)

for(z in 1:ncol(PANEL.MOVES))
  {BUILD=EMPTY;
   STRIP=PROFILE.LIST[,z];
   BUILD[,16]=STRIP;
   for(x in seq(-15,-1,1)) BUILD[1:(length(STRIP)+x),(x+16)]=STRIP[((-1*x)+1):length(STRIP)];
   for(x in seq(1,15,1)) BUILD[(x+1):length(STRIP),(x+16)]=STRIP[1:(length(STRIP)-x)];
   ## now, compare correlation coefficients to the MEAN.PROFILE;
   BEST=c(0,0);
   for(x in seq(-15,15,1)) {A=cor.test(PROFILE.TARGET[,2],BUILD[,(x+16)])[4]$est; if(A>BEST[1]) {BEST[1]=A;BEST[2]=x}};
   PANEL.MOVES[PANEL.NUMBER,z]=BEST[2]}


#####################################################################################################
## next panels:

for(PANEL.NUMBER in 2:13)
 {PANEL=REG.DATA2[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],,];
  PROFILE.LIST=matrix(,dim(REG.DATA)[2],dim(REG.DATA)[3]);
  PANEL.factor=PANEL[,,1];
  for(x in 1:ncol(PANEL.factor)) PANEL.factor[,x]=x;
  PANEL.factor=as.factor(PANEL.factor);
  for(z in 1:dim(PANEL)[3]) PROFILE.LIST[,z]=as.numeric(tapply(PANEL[,,z],PANEL.factor,MEANx));
  ## note:this is reading from the bottom to the top
  PROFILE.TARGET=PROFILE.LIST[,1:2];
  PROFILE.TARGET[,]<-NA;
  PROFILE.TARGET[,2]=as.numeric(tapply(TARGET[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],],PANEL.factor,MEANx));
  EMPTY=matrix(,dim(REG.DATA2)[2],31);
  for(z in 1:ncol(PANEL.MOVES))
    {BUILD=EMPTY;
     STRIP=PROFILE.LIST[,z];
     BUILD[,16]=STRIP;
     for(x in seq(-15,-1,1)) BUILD[1:(length(STRIP)+x),(x+16)]=STRIP[((-1*x)+1):length(STRIP)];
     for(x in seq(1,15,1)) BUILD[(x+1):length(STRIP),(x+16)]=STRIP[1:(length(STRIP)-x)];
     ## now, compare correlation coefficients to the MEAN.PROFILE;
     BEST=c(0,0);
     for(x in seq(-15,15,1)) {A=cor.test(PROFILE.TARGET[,2],BUILD[,(x+16)])[4]$est; if(A>BEST[1]) {BEST[1]=A;BEST[2]=x}};
     PANEL.MOVES[PANEL.NUMBER,z]=BEST[2]}}


#####################################################################################################
## now, use splines to estimate how much to move each "column" of pixels:


## image by image in the stack

REG.DATA3=REG.DATA2
REG.DATA3[,,]<-NA

for(z in 1:dim(REG.DATA3)[3])
 {make.spline=cbind(PANEL.BORDERS[,3],PANEL.MOVES[,z]);
  SPLINE=smooth.spline(make.spline[,1],make.spline[,2], df=7);
  plot(make.spline[,1],make.spline[,2],ylim=c(-10,10));
  predict.range=seq(1,(dim(REG.DATA3)[1]),1);
  SPLINE.PREDICTED=cbind(predict.range,predict(SPLINE,predict.range)$y);
  matlines(SPLINE.PREDICTED[,1],SPLINE.PREDICTED[,2]);
  ## once we round these, they'll be ready to use!
  SPLINE.PREDICTED[,2]=round(SPLINE.PREDICTED[,2])
  matlines(SPLINE.PREDICTED[,1],SPLINE.PREDICTED[,2],col="red");
  for(x in 1:nrow(SPLINE.PREDICTED))
   {move=SPLINE.PREDICTED[x,2];
    if(move==0) REG.DATA3[x,,z]=REG.DATA2[x,,z];
    if(move>0) REG.DATA3[x,(move+1):(dim(REG.DATA3)[2]),z]=REG.DATA2[x,1:(dim(REG.DATA3)[2]-move),z];
    if(move<0) REG.DATA3[x,1:(dim(REG.DATA3)[2]+move),z]=REG.DATA2[x,((-1*move)+1):(dim(REG.DATA3)[2]),z]}}


########################################################
########################################################
########################################################
########################################################  repeat most steps:
########################################################
########################################################


REG.DATA=REG.DATA3


######################################################
##
## for the next step, we've hopefully
## taken care of most of the rotational aspects. 
##
## I've noticed some side-to-side jitter of the retina; we'll try to blunt this by just using the
## (STANDARD DEVIATION) signal within each "coulmn" of the image, which makes the bright spot for the vascular bundle "pop" visually
## (the bright cluster of blood vessels at the optic N head is a nice landmark)

## we'll only use the central 50%
BORDERS=round(c(dim(REG.DATA)[1]*1/4,dim(REG.DATA)[1]*3/4))
CENTRAL.HALF=REG.DATA[BORDERS[1]:BORDERS[2],,]

## the list of mean profiles we'll get...
PROFILE.LIST=matrix(,length(BORDERS[1]:BORDERS[2]),dim(REG.DATA)[3])

## and we need to use tapply
MEANx=function(x) mean(na.rm=T,x)
SDx=function(x) sd(na.rm=T,x)
CENTRAL.HALF.factor=CENTRAL.HALF[,,1]
for(x in 1:nrow(CENTRAL.HALF.factor)) CENTRAL.HALF.factor[x,]=x
CENTRAL.HALF.factor=as.factor(CENTRAL.HALF.factor)

#for(z in 1:dim(CENTRAL.HALF)[3]) PROFILE.LIST[,z]=as.numeric(tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,MEANx))
#for(z in 1:dim(CENTRAL.HALF)[3]) PROFILE.LIST[,z]=as.numeric(tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,SDx)/tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,MEANx))
for(z in 1:dim(CENTRAL.HALF)[3]) PROFILE.LIST[,z]=as.numeric(tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,SDx))

## as you can see, the data are noisey:
show.slice1=1
show.slice2=2
plot(PROFILE.LIST[,show.slice1],pch=19)
matpoints(PROFILE.LIST[,show.slice2],pch=21,col="red")

## anf we'll use the mean as a target:
MEAN.PROFILE=PROFILE.LIST[,1:2]
MEAN.PROFILE[,]<-NA
for(x in 1:nrow(MEAN.PROFILE)) MEAN.PROFILE[x,2]=mean(na.rm=T,PROFILE.LIST[x,])

plot(MEAN.PROFILE[,2],type="l")
matpoints(PROFILE.LIST[,show.slice1],pch=19,col="black")
matpoints(PROFILE.LIST[,show.slice2],pch=21,col="red")

rm(show.slice1,show.slice2)


## now, how to slide one past another to minimize variance?
## allow movements of +/- 30 pixels
## a negative number is to the left. 

EMPTY=matrix(,nrow(PROFILE.LIST),61)

LIST.OF.MOVEMENTS=matrix(,dim(REG.DATA)[3],2)
for(z in 1:ncol(PROFILE.LIST))
  {BUILD=EMPTY;
   STRIP=PROFILE.LIST[,z];
   BUILD[,31]=STRIP;
   for(x in seq(-30,-1,1)) BUILD[1:(length(STRIP)+x),(x+31)]=STRIP[((-1*x)+1):length(STRIP)];
   for(x in seq(1,30,1)) BUILD[(x+1):length(STRIP),(x+31)]=STRIP[1:(length(STRIP)-x)];
   ## now, compare correlation coefficients to the MEAN.PROFILE;
   BEST=c(0,0);
   for(x in seq(-30,30,1)) {A=cor.test(MEAN.PROFILE[,2],BUILD[,(x+31)])[4]$est; if(A>BEST[1]) {BEST[1]=A;BEST[2]=x}}
   LIST.OF.MOVEMENTS[z,2]=BEST[2]}

## now, apply those corrections:

REG.DATA2=REG.DATA
REG.DATA2[,,]<-NA
for(z in 1:nrow(LIST.OF.MOVEMENTS))
 {move=LIST.OF.MOVEMENTS[z,2];
  if(move==0) REG.DATA2[,,z]=REG.DATA[,,z];
  if(move<0) REG.DATA2[1:(dim(REG.DATA2)[1]+move),,z]=REG.DATA[((-1*move)+1):dim(REG.DATA2)[1],,z];
  if(move>0) REG.DATA2[(move+1):dim(REG.DATA2)[1],,z]=REG.DATA[1:(dim(REG.DATA2)[1]-move),,z]}

#image(REG.DATA2[,,1])

#REG.DATAdisplay=REG.DATA2
#REG.DATAdisplay[which(is.na(REG.DATA2))]=0
#f.write.analyze(REG.DATAdisplay[,,],paste(FILENAME,"_REG_step2",sep=""),size="float")
#rm(REG.DATAdisplay)


######################################################
##
## Now, we're going to break the retina into several "panels":
## One panel will be the central 1/4 of the image (6/16 to 10/16ths). The other 12 "panels" will be on each side of the center,
## representing 1/16th of the retina. 
##
## The panels will be shifted vertially to try and match the mean.
## Whatever the proposed vertical movement, we'll submit the whole thing to a smooth spline across the image
## so judge how we should move each "column" of pixels vertically.


TARGET=REG.DATA2[,,1]
to.average.across.z=matrix(,dim(REG.DATA2)[2],dim(REG.DATA2)[3])
for(x in 1:nrow(to.average.across.z)) to.average.across.z[x,]=x;
to.average.across.z=as.factor(to.average.across.z)
for(x in 1:nrow(TARGET)) TARGET[x,]=as.numeric(tapply(REG.DATA2[x,,],to.average.across.z,MEANx));

PANEL.BORDERS=matrix(,13,2)
#PANEL.BORDERS[1,1]=1
#PANEL.BORDERS[13,2]=dim(REG.DATA)[1]
for(x in 1:6)
 {PANEL.BORDERS[x,]=c(1+round(dim(REG.DATA)[1]*(x-1)/16),round(dim(REG.DATA)[1]*(x)/16))}
for(x in 8:13)
 {PANEL.BORDERS[x,]=c(1+round(dim(REG.DATA)[1]*(x+2)/16),round(dim(REG.DATA)[1]*(x+3)/16))}
PANEL.BORDERS[7,1]=PANEL.BORDERS[6,2]+1
PANEL.BORDERS[7,2]=PANEL.BORDERS[8,1]-1
PANEL.BORDERS=cbind(PANEL.BORDERS,PANEL.BORDERS[,1])
PANEL.BORDERS[,3]<-(PANEL.BORDERS[,1]+PANEL.BORDERS[,2])/2
colnames(PANEL.BORDERS)<-c("start","end","mean")


#
#

## we want a matrix of vertical moves for each panel, for each image:
# +1 will be "move up by 1 to match target", -1 will be to move down by 1. 
# check for any distances +/-15 pixels
PANEL.MOVES=matrix(,13,dim(REG.DATA2)[3])

#####################################################################################################
## first panel:
PANEL.NUMBER=1

PANEL=REG.DATA2[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],,]

PROFILE.LIST=matrix(,dim(REG.DATA)[2],dim(REG.DATA)[3])

PANEL.factor=PANEL[,,1]
for(x in 1:ncol(PANEL.factor)) PANEL.factor[,x]=x
PANEL.factor=as.factor(PANEL.factor)

for(z in 1:dim(PANEL)[3]) PROFILE.LIST[,z]=as.numeric(tapply(PANEL[,,z],PANEL.factor,MEANx))
## note:this is reading from the bottom to the top
PROFILE.TARGET=PROFILE.LIST[,1:2]
PROFILE.TARGET[,]<-NA
PROFILE.TARGET[,2]=as.numeric(tapply(TARGET[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],],PANEL.factor,MEANx))

EMPTY=matrix(,dim(REG.DATA2)[2],31)

for(z in 1:ncol(PANEL.MOVES))
  {BUILD=EMPTY;
   STRIP=PROFILE.LIST[,z];
   BUILD[,16]=STRIP;
   for(x in seq(-15,-1,1)) BUILD[1:(length(STRIP)+x),(x+16)]=STRIP[((-1*x)+1):length(STRIP)];
   for(x in seq(1,15,1)) BUILD[(x+1):length(STRIP),(x+16)]=STRIP[1:(length(STRIP)-x)];
   ## now, compare correlation coefficients to the MEAN.PROFILE;
   BEST=c(0,0);
   for(x in seq(-15,15,1)) {A=cor.test(PROFILE.TARGET[,2],BUILD[,(x+16)])[4]$est; if(A>BEST[1]) {BEST[1]=A;BEST[2]=x}};
   PANEL.MOVES[PANEL.NUMBER,z]=BEST[2]}


#####################################################################################################
## next panels:

for(PANEL.NUMBER in 2:13)
 {PANEL=REG.DATA2[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],,];
  PROFILE.LIST=matrix(,dim(REG.DATA)[2],dim(REG.DATA)[3]);
  PANEL.factor=PANEL[,,1];
  for(x in 1:ncol(PANEL.factor)) PANEL.factor[,x]=x;
  PANEL.factor=as.factor(PANEL.factor);
  for(z in 1:dim(PANEL)[3]) PROFILE.LIST[,z]=as.numeric(tapply(PANEL[,,z],PANEL.factor,MEANx));
  ## note:this is reading from the bottom to the top
  PROFILE.TARGET=PROFILE.LIST[,1:2];
  PROFILE.TARGET[,]<-NA;
  PROFILE.TARGET[,2]=as.numeric(tapply(TARGET[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],],PANEL.factor,MEANx));
  EMPTY=matrix(,dim(REG.DATA2)[2],31);
  for(z in 1:ncol(PANEL.MOVES))
    {BUILD=EMPTY;
     STRIP=PROFILE.LIST[,z];
     BUILD[,16]=STRIP;
     for(x in seq(-15,-1,1)) BUILD[1:(length(STRIP)+x),(x+16)]=STRIP[((-1*x)+1):length(STRIP)];
     for(x in seq(1,15,1)) BUILD[(x+1):length(STRIP),(x+16)]=STRIP[1:(length(STRIP)-x)];
     ## now, compare correlation coefficients to the MEAN.PROFILE;
     BEST=c(0,0);
     for(x in seq(-15,15,1)) {A=cor.test(PROFILE.TARGET[,2],BUILD[,(x+16)])[4]$est; if(A>BEST[1]) {BEST[1]=A;BEST[2]=x}};
     PANEL.MOVES[PANEL.NUMBER,z]=BEST[2]}}


#####################################################################################################
## now, use splines to estimate how much to move each "column" of pixels:


## image by image in the stack

REG.DATA3=REG.DATA2
REG.DATA3[,,]<-NA

for(z in 1:dim(REG.DATA3)[3])
 {make.spline=cbind(PANEL.BORDERS[,3],PANEL.MOVES[,z]);
  SPLINE=smooth.spline(make.spline[,1],make.spline[,2], df=7);
  plot(make.spline[,1],make.spline[,2],ylim=c(-10,10));
  predict.range=seq(1,(dim(REG.DATA3)[1]),1);
  SPLINE.PREDICTED=cbind(predict.range,predict(SPLINE,predict.range)$y);
  matlines(SPLINE.PREDICTED[,1],SPLINE.PREDICTED[,2]);
  ## once we round these, they'll be ready to use!
  SPLINE.PREDICTED[,2]=round(SPLINE.PREDICTED[,2])
  matlines(SPLINE.PREDICTED[,1],SPLINE.PREDICTED[,2],col="red");
  for(x in 1:nrow(SPLINE.PREDICTED))
   {move=SPLINE.PREDICTED[x,2];
    if(move==0) REG.DATA3[x,,z]=REG.DATA2[x,,z];
    if(move>0) REG.DATA3[x,(move+1):(dim(REG.DATA3)[2]),z]=REG.DATA2[x,1:(dim(REG.DATA3)[2]-move),z];
    if(move<0) REG.DATA3[x,1:(dim(REG.DATA3)[2]+move),z]=REG.DATA2[x,((-1*move)+1):(dim(REG.DATA3)[2]),z]}}

########################################################
########################################################
########################################################
########################################################  repeat most steps:
########################################################
########################################################


REG.DATA=REG.DATA3


######################################################
##
## for the next step, we've hopefully
## taken care of most of the rotational aspects. 
##
## I've noticed some side-to-side jitter of the retina; we'll try to blunt this by just using the
## (STANDARD DEVIATION) signal within each "coulmn" of the image, which makes the bright spot for the vascular bundle "pop" visually
## (the bright cluster of blood vessels at the optic N head is a nice landmark)

## we'll only use the central 50%
BORDERS=round(c(dim(REG.DATA)[1]*1/4,dim(REG.DATA)[1]*3/4))
CENTRAL.HALF=REG.DATA[BORDERS[1]:BORDERS[2],,]

## the list of mean profiles we'll get...
PROFILE.LIST=matrix(,length(BORDERS[1]:BORDERS[2]),dim(REG.DATA)[3])

## and we need to use tapply
MEANx=function(x) mean(na.rm=T,x)
SDx=function(x) sd(na.rm=T,x)
CENTRAL.HALF.factor=CENTRAL.HALF[,,1]
for(x in 1:nrow(CENTRAL.HALF.factor)) CENTRAL.HALF.factor[x,]=x
CENTRAL.HALF.factor=as.factor(CENTRAL.HALF.factor)

#for(z in 1:dim(CENTRAL.HALF)[3]) PROFILE.LIST[,z]=as.numeric(tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,MEANx))
#for(z in 1:dim(CENTRAL.HALF)[3]) PROFILE.LIST[,z]=as.numeric(tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,SDx)/tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,MEANx))
for(z in 1:dim(CENTRAL.HALF)[3]) PROFILE.LIST[,z]=as.numeric(tapply(CENTRAL.HALF[,,z],CENTRAL.HALF.factor,SDx))

## as you can see, the data are noisey:
show.slice1=1
show.slice2=2
plot(PROFILE.LIST[,show.slice1],pch=19)
matpoints(PROFILE.LIST[,show.slice2],pch=21,col="red")

## anf we'll use the mean as a target:
MEAN.PROFILE=PROFILE.LIST[,1:2]
MEAN.PROFILE[,]<-NA
for(x in 1:nrow(MEAN.PROFILE)) MEAN.PROFILE[x,2]=mean(na.rm=T,PROFILE.LIST[x,])

plot(MEAN.PROFILE[,2],type="l")
matpoints(PROFILE.LIST[,show.slice1],pch=19,col="black")
matpoints(PROFILE.LIST[,show.slice2],pch=21,col="red")

rm(show.slice1,show.slice2)


## now, how to slide one past another to minimize variance?
## allow movements of +/- 30 pixels
## a negative number is to the left. 

EMPTY=matrix(,nrow(PROFILE.LIST),61)

LIST.OF.MOVEMENTS=matrix(,dim(REG.DATA)[3],2)
for(z in 1:ncol(PROFILE.LIST))
  {BUILD=EMPTY;
   STRIP=PROFILE.LIST[,z];
   BUILD[,31]=STRIP;
   for(x in seq(-30,-1,1)) BUILD[1:(length(STRIP)+x),(x+31)]=STRIP[((-1*x)+1):length(STRIP)];
   for(x in seq(1,30,1)) BUILD[(x+1):length(STRIP),(x+31)]=STRIP[1:(length(STRIP)-x)];
   ## now, compare correlation coefficients to the MEAN.PROFILE;
   BEST=c(0,0);
   for(x in seq(-30,30,1)) {A=cor.test(MEAN.PROFILE[,2],BUILD[,(x+31)])[4]$est; if(A>BEST[1]) {BEST[1]=A;BEST[2]=x}}
   LIST.OF.MOVEMENTS[z,2]=BEST[2]}

## now, apply those corrections:

REG.DATA2=REG.DATA
REG.DATA2[,,]<-NA
for(z in 1:nrow(LIST.OF.MOVEMENTS))
 {move=LIST.OF.MOVEMENTS[z,2];
  if(move==0) REG.DATA2[,,z]=REG.DATA[,,z];
  if(move<0) REG.DATA2[1:(dim(REG.DATA2)[1]+move),,z]=REG.DATA[((-1*move)+1):dim(REG.DATA2)[1],,z];
  if(move>0) REG.DATA2[(move+1):dim(REG.DATA2)[1],,z]=REG.DATA[1:(dim(REG.DATA2)[1]-move),,z]}

#image(REG.DATA2[,,1])

#REG.DATAdisplay=REG.DATA2
#REG.DATAdisplay[which(is.na(REG.DATA2))]=0
#f.write.analyze(REG.DATAdisplay[,,],paste(FILENAME,"_REG_step2",sep=""),size="float")
#rm(REG.DATAdisplay)


######################################################
##
## Now, we're going to break the retina into several "panels":
## One panel will be the central 1/4 of the image (6/16 to 10/16ths). The other 12 "panels" will be on each side of the center,
## representing 1/16th of the retina. 
##
## The panels will be shifted vertially to try and match the mean.
## Whatever the proposed vertical movement, we'll submit the whole thing to a smooth spline across the image
## so judge how we should move each "column" of pixels vertically.


TARGET=REG.DATA2[,,1]
to.average.across.z=matrix(,dim(REG.DATA2)[2],dim(REG.DATA2)[3])
for(x in 1:nrow(to.average.across.z)) to.average.across.z[x,]=x;
to.average.across.z=as.factor(to.average.across.z)
for(x in 1:nrow(TARGET)) TARGET[x,]=as.numeric(tapply(REG.DATA2[x,,],to.average.across.z,MEANx));

PANEL.BORDERS=matrix(,13,2)
#PANEL.BORDERS[1,1]=1
#PANEL.BORDERS[13,2]=dim(REG.DATA)[1]
for(x in 1:6)
 {PANEL.BORDERS[x,]=c(1+round(dim(REG.DATA)[1]*(x-1)/16),round(dim(REG.DATA)[1]*(x)/16))}
for(x in 8:13)
 {PANEL.BORDERS[x,]=c(1+round(dim(REG.DATA)[1]*(x+2)/16),round(dim(REG.DATA)[1]*(x+3)/16))}
PANEL.BORDERS[7,1]=PANEL.BORDERS[6,2]+1
PANEL.BORDERS[7,2]=PANEL.BORDERS[8,1]-1
PANEL.BORDERS=cbind(PANEL.BORDERS,PANEL.BORDERS[,1])
PANEL.BORDERS[,3]<-(PANEL.BORDERS[,1]+PANEL.BORDERS[,2])/2
colnames(PANEL.BORDERS)<-c("start","end","mean")


#
#

## we want a matrix of vertical moves for each panel, for each image:
# +1 will be "move up by 1 to match target", -1 will be to move down by 1. 
# check for any distances +/-15 pixels
PANEL.MOVES=matrix(,13,dim(REG.DATA2)[3])

#####################################################################################################
## first panel:
PANEL.NUMBER=1

PANEL=REG.DATA2[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],,]

PROFILE.LIST=matrix(,dim(REG.DATA)[2],dim(REG.DATA)[3])

PANEL.factor=PANEL[,,1]
for(x in 1:ncol(PANEL.factor)) PANEL.factor[,x]=x
PANEL.factor=as.factor(PANEL.factor)

for(z in 1:dim(PANEL)[3]) PROFILE.LIST[,z]=as.numeric(tapply(PANEL[,,z],PANEL.factor,MEANx))
## note:this is reading from the bottom to the top
PROFILE.TARGET=PROFILE.LIST[,1:2]
PROFILE.TARGET[,]<-NA
PROFILE.TARGET[,2]=as.numeric(tapply(TARGET[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],],PANEL.factor,MEANx))

EMPTY=matrix(,dim(REG.DATA2)[2],31)

for(z in 1:ncol(PANEL.MOVES))
  {BUILD=EMPTY;
   STRIP=PROFILE.LIST[,z];
   BUILD[,16]=STRIP;
   for(x in seq(-15,-1,1)) BUILD[1:(length(STRIP)+x),(x+16)]=STRIP[((-1*x)+1):length(STRIP)];
   for(x in seq(1,15,1)) BUILD[(x+1):length(STRIP),(x+16)]=STRIP[1:(length(STRIP)-x)];
   ## now, compare correlation coefficients to the MEAN.PROFILE;
   BEST=c(0,0);
   for(x in seq(-15,15,1)) {A=cor.test(PROFILE.TARGET[,2],BUILD[,(x+16)])[4]$est; if(A>BEST[1]) {BEST[1]=A;BEST[2]=x}};
   PANEL.MOVES[PANEL.NUMBER,z]=BEST[2]}


#####################################################################################################
## next panels:

for(PANEL.NUMBER in 2:13)
 {PANEL=REG.DATA2[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],,];
  PROFILE.LIST=matrix(,dim(REG.DATA)[2],dim(REG.DATA)[3]);
  PANEL.factor=PANEL[,,1];
  for(x in 1:ncol(PANEL.factor)) PANEL.factor[,x]=x;
  PANEL.factor=as.factor(PANEL.factor);
  for(z in 1:dim(PANEL)[3]) PROFILE.LIST[,z]=as.numeric(tapply(PANEL[,,z],PANEL.factor,MEANx));
  ## note:this is reading from the bottom to the top
  PROFILE.TARGET=PROFILE.LIST[,1:2];
  PROFILE.TARGET[,]<-NA;
  PROFILE.TARGET[,2]=as.numeric(tapply(TARGET[PANEL.BORDERS[PANEL.NUMBER,1]:PANEL.BORDERS[PANEL.NUMBER,2],],PANEL.factor,MEANx));
  EMPTY=matrix(,dim(REG.DATA2)[2],31);
  for(z in 1:ncol(PANEL.MOVES))
    {BUILD=EMPTY;
     STRIP=PROFILE.LIST[,z];
     BUILD[,16]=STRIP;
     for(x in seq(-15,-1,1)) BUILD[1:(length(STRIP)+x),(x+16)]=STRIP[((-1*x)+1):length(STRIP)];
     for(x in seq(1,15,1)) BUILD[(x+1):length(STRIP),(x+16)]=STRIP[1:(length(STRIP)-x)];
     ## now, compare correlation coefficients to the MEAN.PROFILE;
     BEST=c(0,0);
     for(x in seq(-15,15,1)) {A=cor.test(PROFILE.TARGET[,2],BUILD[,(x+16)])[4]$est; if(A>BEST[1]) {BEST[1]=A;BEST[2]=x}};
     PANEL.MOVES[PANEL.NUMBER,z]=BEST[2]}}


#####################################################################################################
## now, use splines to estimate how much to move each "column" of pixels:


## image by image in the stack

REG.DATA3=REG.DATA2
REG.DATA3[,,]<-NA

for(z in 1:dim(REG.DATA3)[3])
 {make.spline=cbind(PANEL.BORDERS[,3],PANEL.MOVES[,z]);
  SPLINE=smooth.spline(make.spline[,1],make.spline[,2], df=7);
  plot(make.spline[,1],make.spline[,2],ylim=c(-10,10));
  predict.range=seq(1,(dim(REG.DATA3)[1]),1);
  SPLINE.PREDICTED=cbind(predict.range,predict(SPLINE,predict.range)$y);
  matlines(SPLINE.PREDICTED[,1],SPLINE.PREDICTED[,2]);
  ## once we round these, they'll be ready to use!
  SPLINE.PREDICTED[,2]=round(SPLINE.PREDICTED[,2])
  matlines(SPLINE.PREDICTED[,1],SPLINE.PREDICTED[,2],col="red");
  for(x in 1:nrow(SPLINE.PREDICTED))
   {move=SPLINE.PREDICTED[x,2];
    if(move==0) REG.DATA3[x,,z]=REG.DATA2[x,,z];
    if(move>0) REG.DATA3[x,(move+1):(dim(REG.DATA3)[2]),z]=REG.DATA2[x,1:(dim(REG.DATA3)[2]-move),z];
    if(move<0) REG.DATA3[x,1:(dim(REG.DATA3)[2]+move),z]=REG.DATA2[x,((-1*move)+1):(dim(REG.DATA3)[2]),z]}}



############
############  and we're done!
############


REG.DATAdisplay=REG.DATA3
REG.DATAdisplay[which(is.na(REG.DATAdisplay))]=0
f.write.analyze(REG.DATAdisplay[,,],paste(FILENAME,"_REGISTERED",sep=""),size="float")
rm(REG.DATAdisplay)

TARGET=REG.DATA3[,,1]
to.average.across.z=matrix(,dim(REG.DATA3)[2],dim(REG.DATA3)[3])
for(x in 1:nrow(to.average.across.z)) to.average.across.z[x,]=x;
to.average.across.z=as.factor(to.average.across.z)
for(x in 1:nrow(TARGET)) TARGET[x,]=as.numeric(tapply(REG.DATA3[x,,],to.average.across.z,MEANx));

REG.DATAdisplay=TARGET
REG.DATAdisplay[which(is.na(REG.DATAdisplay))]=0
f.write.analyze(REG.DATAdisplay[,],paste(FILENAME,"_MEAN",sep=""),size="float")
rm(REG.DATAdisplay)

