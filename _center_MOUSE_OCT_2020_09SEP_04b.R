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

FILENAME="Stack"



#
#
#
#####
# 2 #
#####
###
###
### state your resolution in the x (towards/from optic nerve within the retina) and y (following the path of light, traversing layers)
### in microns per pixel
### so, in the original dataset, it was 1.4 um / pixel for X (horizontal) and 1.53 um /pixel for Y (vertical).
### in that case, the next lines should say
#x.RESOLUTION=1.4
#y.RESOLUTION=1.53

x.RESOLUTION=1.4
y.RESOLUTION=1.53

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


#
#
#
#####
# OLD, OLD, OLD #
#####
###
### You need to describe the order of light and dark conditions.
### LIGHT=1, DARK=0
### so, if your images went DARK, LIGHT, DARK, LIGHT, DARK, LIGHT
### the next line should be
#ORDER=c(0,1,0,1,0,1)
### on the other hand, if your images went DARK, then LIGHT, then DARK, DARK, LIGHT, LIGHT
### the next line should be
#ORDER=c(0,1,0,0,1,1)
###
### IMPORTANT!
###  the length of this "ORDER" variable needs to match the number of images in the ANALYZE file
###  one of the first outputs is registered images sorted by light and dark.
###  if you don't see this, it's probably because of a mismatch.
###  that next line, "length(ORDER)" can be checked before you initialize the whole program.
###
### IMPORTANT!
###  as of the OCT2019 (and beyond) versions, the script has been heavily modified to just run one image at a time
###  the code still has "light" and "dark" distinctions, but this is residual from old versions.
ORDER=c(0,1)
length(ORDER)


DATA<-f.read.analyze.volume(paste(FILENAME,".hdr",sep=""))
if(length(dim(DATA))==4) DATA=DATA[,,,1]

# The first image will be used to register the others.
FIRST=DATA[,,1]


REG.DATA=DATA
for(y in 2:dim(DATA)[3])
 {REG.DATA[,,y]=niftyreg.linear(DATA[,,y], FIRST, scope = "rigid", init = NULL, 
                sourceMask = NULL, targetMask = NULL, symmetric = FALSE, 
                nLevels = 3, maxIterations = 5, estimateOnly = FALSE)$image}
REG.DATA[which(as.character(REG.DATA)=="NaN")]=0

f.write.analyze(REG.DATA[,,],paste(FILENAME,"_REG",sep=""),size="float")

LIGHT.DATA=REG.DATA[,,]
#f.write.analyze(LIGHT.DATA[,,],paste(FILENAME,"-ALL",sep=""),size="float")

DARK.DATA=REG.DATA[,,]
#f.write.analyze(DARK.DATA[,,],paste(FILENAME,"-DARK",sep=""),size="float")



#################################################
#################################################
REF=read.table(paste(FILENAME,".txt",sep=""), header=FALSE)
REF=as.matrix(REF)
#### o.k., get REF in the same orientation as DATA  ####
R=t(REF)
R=R[,ncol(R):1]
REF=R

# now, find the center of the optic nerve
# NA-out unmarked space...
R=REF
Xs=R
insert=seq(1,nrow(R),1)
for(x in 1:ncol(Xs)) Xs[,x]=insert
# NA-out unmarked space...
R[which(R<243)]=NA
R[which(R>243)]=NA
# set marked space to 1
R[which(R==243)]=1
#
Xcoords.opticN=Xs*R
opticN.X.position=round(mean(na.rm=T,Xcoords.opticN))


#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################


# now, find the RPE
R=REF
Xs=R
Ys=R
insert=seq(1,nrow(R),1)
for(x in 1:ncol(Xs)) Xs[,x]=insert
insert=seq(1,ncol(R),1)
for(x in 1:nrow(Ys)) Ys[x,]=insert
#
# zero-out unmarked space...
R[which(R<255)]=0
# set marked space to 1
R[which(R>=255)]=1

Xcoords=Xs*R
Ycoords=Ys*R

coordlist=cbind(Xcoords[which(Xcoords>0)],Ycoords[which(Ycoords>0)])
plot(coordlist[,1],coordlist[,2])


MEAN.x=function(x) mean(na.rm=TRUE,x)

## and do a quick check for the retina/vitreous border
R=REF
Xs=R
Ys=R
insert=seq(1,nrow(R),1)
for(x in 1:ncol(Xs)) Xs[,x]=insert
insert=seq(1,ncol(R),1)
for(x in 1:nrow(Ys)) Ys[x,]=insert
#
# zero-out unmarked space...
R[which(R<249)]=0
R[which(R>249)]=0
# set marked space to 1
R[which(R==249)]=1

ret.vit.Xcoords=Xs*R
ret.vit.Ycoords=Ys*R

ret.vit.coordlist=cbind(ret.vit.Xcoords[which(ret.vit.Xcoords>0)],ret.vit.Ycoords[which(ret.vit.Ycoords>0)])
plot(ret.vit.coordlist[,1],ret.vit.coordlist[,2])



## to make the arrays more manageable, note that we won't need to go far below the lowest marker for the RPE,
## nor the highest marker for the retina/vitreous broder
## let's go out no more that 40 microns


LOWER.BOUND=min(coordlist[,2])-ceiling(40/y.RESOLUTION)
UPPER.BOUND=max(ret.vit.coordlist[,2])+ceiling(40/y.RESOLUTION)
  ### new for 2019_08AUG_14:
  ### if the retina is linearized beforehand, this
  ### can make the retina "too good" and its markers +/- 40 microns
  ### fit within a 200 micron span presumed much later in the program.
  ### just modify mupper bound to 
if((UPPER.BOUND-LOWER.BOUND)<200) 
   {ADD=ceiling((200-(UPPER.BOUND-LOWER.BOUND))/2);
    UPPER.BOUND=UPPER.BOUND+ADD+1;
    LOWER.BOUND=LOWER.BOUND-ADD;
    rm(ADD)}

LIGHT=LIGHT.DATA[,LOWER.BOUND:UPPER.BOUND,]
DARK=DARK.DATA[,LOWER.BOUND:UPPER.BOUND,]
REF.DATA=REF
REF=REF[,LOWER.BOUND:UPPER.BOUND]


coordlist[,2]=coordlist[,2]-LOWER.BOUND



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############################
####                       ####
####   Skip linearization  ####
####                       ####
###############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rm(list = ls())
#############################################################################
#############################################################################
#############################################################################
#############################################################################


## next up, need a list of locations along the RPE.

Retina.Points=seq(min(na.rm=T,coordlist[,1]),max(na.rm=T,coordlist[,1]),1)
Retina.Points=Retina.Points-opticN.X.position
Retina.Points=cbind(Retina.Points,Retina.Points)
for(x in 1:nrow(Retina.Points))
 {find=cbind(abs(coordlist[,1]-Retina.Points[x,1]),coordlist[,2])
  find=find[order(find[,1]),];
  Retina.Points[x,2]=mean(na.rm=T,find[1:5,2])}
Retina.Points=cbind(seq(min(na.rm=T,coordlist[,1]),max(na.rm=T,coordlist[,1]),1),Retina.Points,Retina.Points[,1]*x.RESOLUTION)
colnames(Retina.Points)<-c("x.pos.original","x.pos.from.ON.pixel","y.pos.pixel","x.pos.from.ON.micron")


## so we can crop anything >641 microns from the optic nerve
## (that's 625 each way, then another just-above-15 is for the moving window

CROP.START=which(Retina.Points[,4]<(-641))[length(which(Retina.Points[,4]<(-641)))]
CROP.END=which(Retina.Points[,4]>(641))[1]

LIGHT.cropped=LIGHT[CROP.START:CROP.END,,]
DARK.cropped=DARK[CROP.START:CROP.END,,]
REF.cropped=REF[CROP.START:CROP.END,]
Retina.Points.cropped=Retina.Points[CROP.START:CROP.END,]

#image(FLATTENED.LIGHT.RETINA.cropped[,,1])

#############################################################################
#############################################################################
#############################################################################
#############################################################################



#############################################################################
#############################################################################
#############################################################################
#############################################################################

## now, moving-window identification of RPE

## the moving window will be a bit more than (rounded up just beyond) 30 um wide 
window.width.in.pixels=ceiling(30/x.RESOLUTION)

start.move=ceiling(window.width.in.pixels/2)
end.move=(nrow(Retina.Points.cropped)-start.move)-1

MEAN.x=function(x) mean(na.rm=TRUE,x)
MAX.x=function(x) max(na.rm=TRUE,x)
window.factor=matrix(,dim(LIGHT.cropped)[2],window.width.in.pixels);
for(y in 1:nrow(window.factor))
 {window.factor[y,]=y}
window.factor=t(window.factor)
window.factor=as.factor(window.factor)

window.factor.2=matrix(,dim(LIGHT.cropped)[2],3);
for(y in 1:nrow(window.factor.2))
 {window.factor.2[y,]=y}
window.factor.2=t(window.factor.2)
window.factor.2=as.factor(window.factor.2)


image(DARK.cropped[,,1]+REF.cropped)

  ## we want OLM position [peak signal near 254], 
  ##  RPE position [peak signal near 255]
  ##  RPE signal
  ##  BM signal
  ## and typical signal (from ~10 to 30 microns interior to the OLM) for the ONL
  ##  
  ## and where availalable:
  ##  retina/vitreous border
  ##  RNFL/GCL border
  ##  ONL/OPL border
  ##  INL/IPL border
RPE.POSITION.LIGHT=matrix(,dim(LIGHT.cropped)[1],dim(LIGHT.cropped)[3])
RPE.SIGNAL.LIGHT=RPE.POSITION.LIGHT
OLM.POSITION.LIGHT=RPE.POSITION.LIGHT
TYPICAL.ONL.SIGNAL.LIGHT=RPE.POSITION.LIGHT
INL.IPL.POSITION.LIGHT=RPE.POSITION.LIGHT
ONL.OPL.POSITION.LIGHT=RPE.POSITION.LIGHT
RNFL.GCL.POSITION.LIGHT=RPE.POSITION.LIGHT
RNFL.PEAK.SIGNAL.LIGHT=RPE.POSITION.LIGHT
VITREOUS.RETINA.POSITION.LIGHT=RPE.POSITION.LIGHT
RPE.POSITION.DARK=matrix(,dim(DARK.cropped)[1],dim(DARK.cropped)[3])
RPE.SIGNAL.DARK=RPE.POSITION.DARK
TYPICAL.ONL.SIGNAL.DARK=RPE.POSITION.DARK
INL.IPL.POSITION.DARK=RPE.POSITION.DARK
ONL.OPL.POSITION.DARK=RPE.POSITION.DARK
RNFL.GCL.POSITION.DARK=RPE.POSITION.DARK
VITREOUS.RETINA.POSITION.DARK=RPE.POSITION.DARK
OLM.POSITION.DARK=RPE.POSITION.DARK
RNFL.PEAK.SIGNAL.DARK=RPE.POSITION.DARK

LIGHT.INDEX.2=cbind(seq(1,dim(LIGHT)[3],1),seq(1,dim(LIGHT)[3],1))
DARK.INDEX.2=cbind(seq(1,dim(DARK)[3],1),seq(1,dim(DARK)[3],1))

## down about 40 microns from the RPE, and up 300
startMARK=floor(40/y.RESOLUTION)*y.RESOLUTION*-1
endMARK=ceiling(300/y.RESOLUTION)*y.RESOLUTION


## do left side, then right side...
break.left=which(Retina.Points.cropped[,4]>(-125))[1]
break.right=which(Retina.Points.cropped[,4]<125)[length(which(Retina.Points.cropped[,4]<125))]


#### 
#### 
#### 
#### 
#### 
#### 
#### 
#### 
#### 
#### BIG NOTE:
#### 
####    For mice, it looks like the basement membrane is a much better landmark than peak RPE. 
####     so, I'm using the variable "RPE" to refer to this (I know, it's bound to confuse, but I don't want to #
###        painstakingly update the variable names. If it helps you think about it, "RPE" means "the exterior border of the RPE")
#### 
#### 
#### 
#### 
#### 
#### 
#### 

for(x in start.move:break.left)
 {light.window=LIGHT.cropped[(x-floor(window.width.in.pixels/2)+1):(x+floor(window.width.in.pixels/2)),,];
  light.profile=matrix(,dim(light.window)[2],nrow(LIGHT.INDEX.2));
  for(y in 1:dim(light.window)[3])
   {light.is=light.window[,,y];
    light.profile[,y]=as.vector(tapply(light.is,window.factor,MEAN.x))};
  light.profile=light.profile[2:nrow(light.profile),];
  #plot(light.profile[,2])
  dark.window=DARK.cropped[(x-floor(window.width.in.pixels/2)+1):(x+floor(window.width.in.pixels/2)),,];
  dark.profile=matrix(,dim(dark.window)[2],nrow(DARK.INDEX.2));
  for(y in 1:dim(dark.window)[3])
   {dark.is=dark.window[,,y];
    dark.profile[,y]=as.vector(tapply(dark.is,window.factor,MEAN.x))};
  dark.profile=dark.profile[2:nrow(dark.profile),];
  plot(dark.profile[,2])
  #
  #profile is organized from out to in. 
  # 
  # get my markings...
  markings.window=REF.cropped[(x-1):(x+1),];
  markings.profile=matrix(,(nrow(light.profile)+1),2);
  markings.profile[,2]=as.vector(tapply(markings.window,window.factor.2,MAX.x));
  markings.profile[which(markings.profile[,2]<249),2]=0;
  #matlines(markings.profile[,2],col="red");
  markings.profile[,1]=seq(1,nrow(markings.profile),1)*y.RESOLUTION;
  markings.profile=markings.profile[2:nrow(markings.profile),];
  #
  ## we want OLM position [peak signal near 254], 
  ##  RPE position [peak signal near 255]
  ##  RPE signal
  ## and typical signal (from ~10 to 30 microns interior to the OLM) for the ONL
  ##  
  ## and where availalable:
  ##  retina/vitreous border
  ##  RNFL/GCL border
  ##  ONL/OPL border
  ##  INL/IPL border
  #
  plot(markings.profile[,1],dark.profile[,1]);
  matlines(markings.profile[,1],markings.profile[,2]);
  RPE.C=round(mean(which(markings.profile[,2]==255)));
  ## look back ~8 microns and forward ~10, and find the local (based non 2rd order polnomial)
  ## 
           ## FOR MICE:
           ## two step: find the local peak, and then go back and pretend that *that* is where the operator drew a line
  ## look back ~8 microns and forward ~10, and find the peak (based non 2rd order polnomial)
  RPE.CHECK=markings.profile[seq(RPE.C-6,RPE.C+6,1),1];
  RPE.dark.CHECK=dark.profile[seq(RPE.C-6,RPE.C+6,1),];
  RPE.CHECK.2=markings.profile[seq(RPE.C-12,RPE.C+6,1),1];
  RPE.dark.CHECK.2=dark.profile[seq(RPE.C-12,RPE.C+6,1),];
  matlines(RPE.CHECK,RPE.dark.CHECK[,1],col="red")
  #matlines(RPE.CHECK,RPE.dark.CHECK[,1],col="red")
  ## 
           ## FOR MICE:
           ## the way the line is drawn, a local minimum exterior to the RPE line
           ## will consitently be basement membrane, and I can cut out anything exterior to that. 
  for(y in 1:ncol(RPE.dark.CHECK))
    {TEST=cbind(RPE.CHECK,RPE.dark.CHECK[,y]);
     ## I know the line is drawn at whatever's at the 7th row of TEST
     cut=which(TEST[1:7,2]==min(na.rm=T,TEST[1:7,2]));
     ## assuming the min value in that range is NOT the one at position, discard stuff below the basement membrane
     if((cut>1)&(cut<7)) TEST=TEST[(cut-1):nrow(TEST),];
     ##
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2))$coef;
     peak=MODEL[2]/(MODEL[3]*(-2)); ##gives vertex x position;
     ## if the parabola is in the wrong direction, just use the location of the highest value.
     if(MODEL[3]>0) peak=TEST[which(TEST[,2]==max(na.rm=T,TEST[,2])),1];
     ## if peak is more than 2 microns from the swath identified, just use upper or lower bound
     if(peak>max(na.rm=T,TEST[,1]+2)) peak=max(na.rm=T,TEST[,1]+(y.RESOLUTION/2));
     if(peak<min(na.rm=T,TEST[,1]-2)) peak=min(na.rm=T,TEST[,1]-(y.RESOLUTION/2));
     #if(y==1) {matlines(TEST[,1],(MODEL[3]*TEST[,1]^2+MODEL[2]*TEST[,1]+MODEL[1]),lwd=3);abline(v=peak,col="red",lwd=2)};
     peak.value=MODEL[3]*peak^2+MODEL[2]*peak+MODEL[1];
     #RPE.POSITION.DARK[x,y]=peak;
     RPE.SIGNAL.DARK[x,y]=peak.value[1];
     TEST=cbind(RPE.CHECK.2,RPE.dark.CHECK.2[,y]);
     # grab everything behind the peak
     testROW=which(abs(TEST[,1]-peak)==min(abs(TEST[,1]-peak)));
     testROW=testROW[length(testROW)]; 
     if(testROW<4) testROW=4;
     TEST=TEST[1:testROW,];
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2)+I(TEST[,1]^3))$coef;
     DERIV=seq(TEST[1,1],TEST[nrow(TEST),1],(TEST[nrow(TEST),1]-TEST[1,1])/30);
     DERIV=cbind(DERIV, (MODEL[4]*3*DERIV^2+MODEL[3]*2*DERIV+MODEL[2]));
     # 
     DERIV=cbind(DERIV,sign(DERIV[,2]));
     DERIV=cbind(DERIV,c(DERIV[2:nrow(DERIV),3],DERIV[nrow(DERIV),3]));
     DERIV[,4]<-ifelse(DERIV[,3]<DERIV[,4],1,1000);
     DERIV[,4]=DERIV[,4]+c(seq(30,0,-1));
     near=which(DERIV[,4]==min(na.rm=T,DERIV[,4]))[1]; 
     if(near==nrow(DERIV)) near=nrow(DERIV)-1;
     DERIV2=as.matrix(DERIV[c(near,near+1),c(1,2,3)]);
     if(sum(DERIV2[,3])==0) troughISH=summary(lm(DERIV2[,1]~DERIV2[,2]))$coef[1,1] else troughISH=DERIV2[1,1];
     ## look within 1.9 pixels of presumed peak
     finalize=TEST[which( (TEST[,1]>(troughISH-(y.RESOLUTION*1.9)) & (TEST[,1]<(troughISH+(y.RESOLUTION*1.9)) ))),];
     if(length(finalize)>2) trough=mean(finalize[which(finalize[,2]==min(finalize[,2])),1]) else trough=finalize[1];
     RPE.POSITION.DARK[x,y]=trough;
     if(y==1) abline(v=trough,col="red",lwd=2)};
  ##
  OLM.C=round(mean(which(markings.profile[,2]==254)));
           ## FOR MICE:
           ## the operator should be within a couple pixels (really, within 1, but the registration step may have been imperfect)
           ## less confusion if we grab a bit extra interior... data are sharp-enough
           ## that the usual polynomial function to seach for the peak 
           ## then just puts us within 1 pixel of the actual peak, and we use the local max. 
           ## 
  ## look back ~8 microns and forward ~8, and find the position of lowest slope (usually a small peak) (based non 3rd order polnomial)
  OLM.CHECK=markings.profile[seq(OLM.C-5,OLM.C+5,1),1];
  OLM.dark.CHECK=dark.profile[seq(OLM.C-5,OLM.C+5,1),];
  matlines(OLM.CHECK,OLM.dark.CHECK[,1],col="red");
  for(y in 1:ncol(OLM.dark.CHECK))
    {TEST=cbind(OLM.CHECK,OLM.dark.CHECK[,y]);
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2)+I(TEST[,1]^3))$coef;
     DERIV=seq(TEST[1,1],TEST[11,1],(TEST[11,1]-TEST[1,1])/30);
     DERIV=cbind(DERIV, (MODEL[4]*3*DERIV^2+MODEL[3]*2*DERIV+MODEL[2]));
     # 
     DERIV=cbind(DERIV,sign(DERIV[,2]));
     DERIV=cbind(DERIV,c(DERIV[2:nrow(DERIV),3],DERIV[nrow(DERIV),3]));
     DERIV[,4]<-ifelse(DERIV[,3]>DERIV[,4],1,1000);
     DERIV[,4]=DERIV[,4]+c(seq(15,0,-1),seq(1,15,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
     near=which(DERIV[,4]==min(na.rm=T,DERIV[,4]))[1]; 
     if(near==nrow(DERIV)) near=nrow(DERIV)-1;
     DERIV2=as.matrix(DERIV[c(near,near+1),c(1,2,3)]);
     if(sum(DERIV2[,3])==0) peakISH=summary(lm(DERIV2[,1]~DERIV2[,2]))$coef[1,1] else peakISH=DERIV2[1,1];
     ## look within 1.9 pixels of presumed peak
     finalize=TEST[which( (TEST[,1]>(peakISH-(y.RESOLUTION*1.9)) & (TEST[,1]<(peakISH+(y.RESOLUTION*1.9)) ))),];
     if(length(finalize)>2) peak=mean(finalize[which(finalize[,2]==max(finalize[,2])),1]) else peak=finalize[1];
     if(y==1) {matlines(TEST[,1],(MODEL[4]*TEST[,1]^3+MODEL[3]*TEST[,1]^2+MODEL[2]*TEST[,1]+MODEL[1]),lwd=3);abline(v=peak,col="red",lwd=2)};
     OLM.POSITION.DARK[x,y]=peak;
     ##
     ## now look ~12 to ~48 microns interior (six to 25 pix interior)
     ##
     startONL=which(abs(markings.profile[,1]-peak)==min(abs(markings.profile[,1]-peak)))[1];
     TYPICAL.ONL.SIGNAL.DARK[x,y]=mean(na.rm=TRUE, dark.profile[(startONL-6):(startONL-25),y])};
  ##  
     ## retype that for light
  ## look back ~8 microns and forward ~10, and find the local (based non 2rd order polnomial)
  RPE.light.CHECK=light.profile[seq(RPE.C-6,RPE.C+6,1),];
  RPE.light.CHECK.2=light.profile[seq(RPE.C-12,RPE.C+6,1),];
  #matlines(RPE.CHECK,RPE.light.CHECK[,1],col="red")
  ## 
           ## FOR MICE:
           ## the way the line is drawn, a local minimum exterior to the RPE line
           ## will consitently be basement membrane, and I can cut out anything exterior to that. 
  for(y in 1:ncol(RPE.light.CHECK))
    {TEST=cbind(RPE.CHECK,RPE.light.CHECK[,y]);
     ## I know the line is drawn at whatever's at the 7th row of TEST
     cut=which(TEST[1:7,2]==min(na.rm=T,TEST[1:7,2]));
     ## assuming the min value in that range is NOT the one at position, discard stuff below the basement membrane
     if((cut>1)&(cut<7)) TEST=TEST[(cut-1):nrow(TEST),];
     ##
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2))$coef;
     peak=MODEL[2]/(MODEL[3]*(-2)); ##gives vertex x position;
     ## if the parabola is in the wrong direction, just use the location of the highest value.
     if(MODEL[3]>0) peak=TEST[which(TEST[,2]==max(na.rm=T,TEST[,2])),1];
     ## if peak is more than 2 microns from the swath identified, just use upper or lower bound
     if(peak>max(na.rm=T,TEST[,1]+2)) peak=max(na.rm=T,TEST[,1]+(y.RESOLUTION/2));
     if(peak<min(na.rm=T,TEST[,1]-2)) peak=min(na.rm=T,TEST[,1]-(y.RESOLUTION/2));
     #if(y==1) {matlines(TEST[,1],(MODEL[3]*TEST[,1]^2+MODEL[2]*TEST[,1]+MODEL[1]),lwd=3);abline(v=peak,col="red",lwd=2)};
     peak.value=MODEL[3]*peak^2+MODEL[2]*peak+MODEL[1];
     #RPE.POSITION.LIGHT[x,y]=peak;
     RPE.SIGNAL.LIGHT[x,y]=peak.value[1];
     TEST=cbind(RPE.CHECK.2,RPE.light.CHECK.2[,y]);
     # grab everything behind the peak
     testROW=which(abs(TEST[,1]-peak)==min(abs(TEST[,1]-peak))); 
     testROW=testROW[length(testROW)]; 
     if(testROW<4) testROW=4;
     TEST=TEST[1:testROW,];
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2)+I(TEST[,1]^3))$coef;
     DERIV=seq(TEST[1,1],TEST[nrow(TEST),1],(TEST[nrow(TEST),1]-TEST[1,1])/30);
     DERIV=cbind(DERIV, (MODEL[4]*3*DERIV^2+MODEL[3]*2*DERIV+MODEL[2]));
     # 
     DERIV=cbind(DERIV,sign(DERIV[,2]));
     DERIV=cbind(DERIV,c(DERIV[2:nrow(DERIV),3],DERIV[nrow(DERIV),3]));
     DERIV[,4]<-ifelse(DERIV[,3]<DERIV[,4],1,1000);
     DERIV[,4]=DERIV[,4]+c(seq(30,0,-1));
     near=which(DERIV[,4]==min(na.rm=T,DERIV[,4]))[1]; 
     if(near==nrow(DERIV)) near=nrow(DERIV)-1;
     DERIV2=as.matrix(DERIV[c(near,near+1),c(1,2,3)]);
     if(sum(DERIV2[,3])==0) troughISH=summary(lm(DERIV2[,1]~DERIV2[,2]))$coef[1,1] else troughISH=DERIV2[1,1];
     ## look within 1.9 pixels of presumed peak
     finalize=TEST[which( (TEST[,1]>(troughISH-(y.RESOLUTION*1.9)) & (TEST[,1]<(troughISH+(y.RESOLUTION*1.9)) ))),];
     if(length(finalize)>2) trough=mean(finalize[which(finalize[,2]==min(finalize[,2])),1]) else trough=finalize[1];
     RPE.POSITION.LIGHT[x,y]=trough};
  ##
  OLM.light.CHECK=light.profile[seq(OLM.C-5,OLM.C+5,1),];
  #matlines(OLM.CHECK,OLM.dark.CHECK[,1],col="red");
  for(y in 1:ncol(OLM.dark.CHECK))
    {TEST=cbind(OLM.CHECK,OLM.dark.CHECK[,y]);
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2)+I(TEST[,1]^3))$coef;
     DERIV=seq(TEST[1,1],TEST[11,1],(TEST[11,1]-TEST[1,1])/30);
     DERIV=cbind(DERIV, (MODEL[4]*3*DERIV^2+MODEL[3]*2*DERIV+MODEL[2]));
     # 
     DERIV=cbind(DERIV,sign(DERIV[,2]));
     DERIV=cbind(DERIV,c(DERIV[2:nrow(DERIV),3],DERIV[nrow(DERIV),3]));
     DERIV[,4]<-ifelse(DERIV[,3]>DERIV[,4],1,1000);
     DERIV[,4]=DERIV[,4]+c(seq(15,0,-1),seq(1,15,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
     near=which(DERIV[,4]==min(na.rm=T,DERIV[,4]))[1]; 
     if(near==nrow(DERIV)) near=nrow(DERIV)-1;
     DERIV2=as.matrix(DERIV[c(near,near+1),c(1,2,3)]);
     if(sum(DERIV2[,3])==0) peakISH=summary(lm(DERIV2[,1]~DERIV2[,2]))$coef[1,1] else peakISH=DERIV2[1,1];
     ## look within 1.9 pixels of presumed peak
     finalize=TEST[which( (TEST[,1]>(peakISH-(y.RESOLUTION*1.9)) & (TEST[,1]<(peakISH+(y.RESOLUTION*1.9)) ))),];
     if(length(finalize)>2) peak=mean(finalize[which(finalize[,2]==max(finalize[,2])),1]) else peak=finalize[1];
     #if(y==1) {matlines(TEST[,1],(MODEL[4]*TEST[,1]^3+MODEL[3]*TEST[,1]^2+MODEL[2]*TEST[,1]+MODEL[1]),lwd=3);abline(v=peak,col="red",lwd=2)};
     OLM.POSITION.LIGHT[x,y]=peak;
     ##
     ## now look ~12 to ~48 microns interior (six to 25 pix interior)
     ##
     startONL=which(abs(markings.profile[,1]-peak)==min(abs(markings.profile[,1]-peak)))[1];
     TYPICAL.ONL.SIGNAL.LIGHT[x,y]=mean(na.rm=TRUE, light.profile[(startONL-6):(startONL-25),y])};
  ## 
  ## 
  ## CONTINUE ONLY IF WE HAVE ALL THE OTHER BORDERS PRESENT
  CONTINUE=1
  if(length(which(markings.profile[,2]==249))==0) CONTINUE=0;
  if(length(which(markings.profile[,2]==250))==0) CONTINUE=0;
  if(length(which(markings.profile[,2]==252))==0) CONTINUE=0;
  if(length(which(markings.profile[,2]==253))==0) CONTINUE=0;
  if(CONTINUE==1)
   {ONL.OPL.C=round(mean(which(markings.profile[,2]==253)));
    ## look back ~24 microns and forward ~24, and use the nearest-to-marking half-height to get the right position
    ONL.OPL.CHECK=markings.profile[seq(ONL.OPL.C-12,ONL.OPL.C+12,1),1];
    ONL.OPL.dark.CHECK=dark.profile[seq(ONL.OPL.C-12,ONL.OPL.C+12,1),];
    for(y in 1:ncol(ONL.OPL.dark.CHECK))
      {TEST=cbind(ONL.OPL.CHECK,ONL.OPL.dark.CHECK[,y]);
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       TEST[,5]=TEST[,5]+c(seq(12,0,-1),seq(1,12,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) {abline(v=border,col="red",lwd=2)};
       ONL.OPL.POSITION.DARK[x,y]=border};
    INL.IPL.C=round(mean(which(markings.profile[,2]==252)));
    ## look back ~24 microns and forward ~24, and use the nearest-to-marking half-height to get the right position
    INL.IPL.CHECK=markings.profile[seq(INL.IPL.C-12,INL.IPL.C+12,1),1];
    INL.IPL.dark.CHECK=dark.profile[seq(INL.IPL.C-12,INL.IPL.C+12,1),];
    for(y in 1:ncol(INL.IPL.dark.CHECK))
      {TEST=cbind(INL.IPL.CHECK,INL.IPL.dark.CHECK[,y]);
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       TEST[,5]=TEST[,5]+c(seq(12,0,-1),seq(1,12,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) {abline(v=border,col="red",lwd=2)};
       INL.IPL.POSITION.DARK[x,y]=border};
    RNFL.GCL.C=round(mean(which(markings.profile[,2]==250)));
    ## look back ~24 microns and forward ~16, and use the nearest-to-marking half-height to get the right position
    RNFL.GCL.CHECK=markings.profile[seq(RNFL.GCL.C-12,RNFL.GCL.C+8,1),1];
    RNFL.GCL.dark.CHECK=dark.profile[seq(RNFL.GCL.C-12,RNFL.GCL.C+8,1),];
    for(y in 1:ncol(RNFL.GCL.dark.CHECK))
      {TEST=cbind(RNFL.GCL.CHECK,RNFL.GCL.dark.CHECK[,y]);
       ## for GCL, toss everything out past the max
       TEST=TEST[1:which(TEST[,2]==max(TEST[3:nrow(TEST),2]))[1],];
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       ## unlike the last two broders, this is biased against more-exterior crossings (since we aleady removed anything past the peak. 
       TEST[,5]=TEST[,5]+seq(nrow(TEST),1,-1);
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) abline(v=border,col="red",lwd=2);
       RNFL.GCL.POSITION.DARK[x,y]=border};
    VITREOUS.RETINA.C=round(mean(which(markings.profile[,2]==249)));
    VITREOUS.RETINA.CHECK=markings.profile[seq(VITREOUS.RETINA.C-12,VITREOUS.RETINA.C+12,1),1];
    VITREOUS.RETINA.dark.CHECK=dark.profile[seq(VITREOUS.RETINA.C-12,VITREOUS.RETINA.C+12,1),];
    for(y in 1:ncol(VITREOUS.RETINA.dark.CHECK))
      {TEST=cbind(VITREOUS.RETINA.CHECK,VITREOUS.RETINA.dark.CHECK[,y]);
       ## for vitreous, toss everything out before the max
       TEST=TEST[which(TEST[,2]==max(TEST[1:(nrow(TEST)-2),2]))[1]:nrow(TEST),];
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       ## unlike the last two broders, this is biased against more-exterior (from the perspective of the eye center) crossings
       TEST[,5]=TEST[,5]+seq(nrow(TEST),1,-1);
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) abline(v=border,col="red",lwd=2);
       VITREOUS.RETINA.POSITION.DARK[x,y]=border};
    ### duplicate for light. 
    ## look back ~24 microns and forward ~12, and find the position of lowest slope (usually a small peak) (based non 3rd order polnomial)
    ONL.OPL.light.CHECK=light.profile[seq(ONL.OPL.C-12,ONL.OPL.C+12,1),];
    for(y in 1:ncol(ONL.OPL.light.CHECK))
      {TEST=cbind(ONL.OPL.CHECK,ONL.OPL.light.CHECK[,y]);
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       TEST[,5]=TEST[,5]+c(seq(12,0,-1),seq(1,12,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) {abline(v=border,col="red",lwd=2)};
       ONL.OPL.POSITION.LIGHT[x,y]=border};
    INL.IPL.light.CHECK=light.profile[seq(INL.IPL.C-12,INL.IPL.C+12,1),];
    for(y in 1:ncol(INL.IPL.light.CHECK))
      {TEST=cbind(INL.IPL.CHECK,INL.IPL.light.CHECK[,y]);
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       TEST[,5]=TEST[,5]+c(seq(12,0,-1),seq(1,12,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) {abline(v=border,col="red",lwd=2)};
       INL.IPL.POSITION.LIGHT[x,y]=border};
    RNFL.GCL.light.CHECK=light.profile[seq(RNFL.GCL.C-12,RNFL.GCL.C+8,1),];
    for(y in 1:ncol(RNFL.GCL.light.CHECK))
      {TEST=cbind(RNFL.GCL.CHECK,RNFL.GCL.light.CHECK[,y]);
       ## for GCL, toss everything out past the max
       TEST=TEST[1:which(TEST[,2]==max(TEST[3:nrow(TEST),2]))[1],];
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       ## unlike the last two broders, this is biased against more-exterior crossings (since we aleady removed anything past the peak. 
       TEST[,5]=TEST[,5]+seq(nrow(TEST),1,-1);
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) abline(v=border,col="red",lwd=2);
       RNFL.GCL.POSITION.LIGHT[x,y]=border};
    VITREOUS.RETINA.light.CHECK=light.profile[seq(VITREOUS.RETINA.C-12,VITREOUS.RETINA.C+12,1),];
    for(y in 1:ncol(VITREOUS.RETINA.light.CHECK))
      {TEST=cbind(VITREOUS.RETINA.CHECK,VITREOUS.RETINA.light.CHECK[,y]);
       ## for vitreous, toss everything out before the max
       TEST=TEST[which(TEST[,2]==max(TEST[1:(nrow(TEST)-2),2]))[1]:nrow(TEST),];
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       ## unlike the last two broders, this is biased against more-exterior (from the perspective of the eye center) crossings
       TEST[,5]=TEST[,5]+seq(nrow(TEST),1,-1);
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) abline(v=border,col="red",lwd=2);
       VITREOUS.RETINA.POSITION.LIGHT[x,y]=border}}}

### now, the right half.

for(x in break.right:end.move)
 {light.window=LIGHT.cropped[(x-floor(window.width.in.pixels/2)+1):(x+floor(window.width.in.pixels/2)),,];
  light.profile=matrix(,dim(light.window)[2],nrow(LIGHT.INDEX.2));
  for(y in 1:dim(light.window)[3])
   {light.is=light.window[,,y];
    light.profile[,y]=as.vector(tapply(light.is,window.factor,MEAN.x))};
  light.profile=light.profile[2:nrow(light.profile),];
  #plot(light.profile[,2])
  dark.window=DARK.cropped[(x-floor(window.width.in.pixels/2)+1):(x+floor(window.width.in.pixels/2)),,];
  dark.profile=matrix(,dim(dark.window)[2],nrow(DARK.INDEX.2));
  for(y in 1:dim(dark.window)[3])
   {dark.is=dark.window[,,y];
    dark.profile[,y]=as.vector(tapply(dark.is,window.factor,MEAN.x))};
  dark.profile=dark.profile[2:nrow(dark.profile),];
  plot(dark.profile[,2])
  #
  #profile is organized from out to in. 
  # 
  # get my markings...
  markings.window=REF.cropped[(x-1):(x+1),];
  markings.profile=matrix(,(nrow(light.profile)+1),2);
  markings.profile[,2]=as.vector(tapply(markings.window,window.factor.2,MAX.x));
  markings.profile[which(markings.profile[,2]<249),2]=0;
  #matlines(markings.profile[,2],col="red");
  markings.profile[,1]=seq(1,nrow(markings.profile),1)*y.RESOLUTION;
  markings.profile=markings.profile[2:nrow(markings.profile),];
  #
  ## we want OLM position [peak signal near 254], 
  ##  RPE position [peak signal near 255]
  ##  RPE signal
  ## and typical signal (from ~10 to 30 microns interior to the OLM) for the ONL
  ##  
  ## and where availalable:
  ##  retina/vitreous border
  ##  RNFL/GCL border
  ##  ONL/OPL border
  ##  INL/IPL border
  #
  plot(markings.profile[,1],dark.profile[,1]);
  matlines(markings.profile[,1],markings.profile[,2]);
  RPE.C=round(mean(which(markings.profile[,2]==255)));
  ## look back ~8 microns and forward ~10, and find the local (based non 2rd order polnomial)
  ## 
           ## FOR MICE:
           ## two step: find the local peak, and then go back and pretend that *that* is where the operator drew a line
  ## look back ~8 microns and forward ~10, and find the peak (based non 2rd order polnomial)
  RPE.CHECK=markings.profile[seq(RPE.C-6,RPE.C+6,1),1];
  RPE.dark.CHECK=dark.profile[seq(RPE.C-6,RPE.C+6,1),];
  RPE.CHECK.2=markings.profile[seq(RPE.C-12,RPE.C+6,1),1];
  RPE.dark.CHECK.2=dark.profile[seq(RPE.C-12,RPE.C+6,1),];
  matlines(RPE.CHECK,RPE.dark.CHECK[,1],col="red")
  #matlines(RPE.CHECK,RPE.dark.CHECK[,1],col="red")
  ## 
           ## FOR MICE:
           ## the way the line is drawn, a local minimum exterior to the RPE line
           ## will consitently be basement membrane, and I can cut out anything exterior to that. 
  for(y in 1:ncol(RPE.dark.CHECK))
    {TEST=cbind(RPE.CHECK,RPE.dark.CHECK[,y]);
     ## I know the line is drawn at whatever's at the 7th row of TEST
     cut=which(TEST[1:7,2]==min(na.rm=T,TEST[1:7,2]));
     ## assuming the min value in that range is NOT the one at position, discard stuff below the basement membrane
     if((cut>1)&(cut<7)) TEST=TEST[(cut-1):nrow(TEST),];
     ##
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2))$coef;
     peak=MODEL[2]/(MODEL[3]*(-2)); ##gives vertex x position;
     ## if the parabola is in the wrong direction, just use the location of the highest value.
     if(MODEL[3]>0) peak=TEST[which(TEST[,2]==max(na.rm=T,TEST[,2])),1];
     ## if peak is more than 2 microns from the swath identified, just use upper or lower bound
     if(peak>max(na.rm=T,TEST[,1]+2)) peak=max(na.rm=T,TEST[,1]+(y.RESOLUTION/2));
     if(peak<min(na.rm=T,TEST[,1]-2)) peak=min(na.rm=T,TEST[,1]-(y.RESOLUTION/2));
     #if(y==1) {matlines(TEST[,1],(MODEL[3]*TEST[,1]^2+MODEL[2]*TEST[,1]+MODEL[1]),lwd=3);abline(v=peak,col="red",lwd=2)};
     peak.value=MODEL[3]*peak^2+MODEL[2]*peak+MODEL[1];
     #RPE.POSITION.DARK[x,y]=peak;
     RPE.SIGNAL.DARK[x,y]=peak.value[1];
     TEST=cbind(RPE.CHECK.2,RPE.dark.CHECK.2[,y]);
     # grab everything behind the peak
     testROW=which(abs(TEST[,1]-peak)==min(abs(TEST[,1]-peak)));
     testROW=testROW[length(testROW)]; 
     if(testROW<4) testROW=4;
     TEST=TEST[1:testROW,];
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2)+I(TEST[,1]^3))$coef;
     DERIV=seq(TEST[1,1],TEST[nrow(TEST),1],(TEST[nrow(TEST),1]-TEST[1,1])/30);
     DERIV=cbind(DERIV, (MODEL[4]*3*DERIV^2+MODEL[3]*2*DERIV+MODEL[2]));
     # 
     DERIV=cbind(DERIV,sign(DERIV[,2]));
     DERIV=cbind(DERIV,c(DERIV[2:nrow(DERIV),3],DERIV[nrow(DERIV),3]));
     DERIV[,4]<-ifelse(DERIV[,3]<DERIV[,4],1,1000);
     DERIV[,4]=DERIV[,4]+c(seq(30,0,-1));
     near=which(DERIV[,4]==min(na.rm=T,DERIV[,4]))[1]; 
     if(near==nrow(DERIV)) near=nrow(DERIV)-1;
     DERIV2=as.matrix(DERIV[c(near,near+1),c(1,2,3)]);
     if(sum(DERIV2[,3])==0) troughISH=summary(lm(DERIV2[,1]~DERIV2[,2]))$coef[1,1] else troughISH=DERIV2[1,1];
     ## look within 1.9 pixels of presumed peak
     finalize=TEST[which( (TEST[,1]>(troughISH-(y.RESOLUTION*1.9)) & (TEST[,1]<(troughISH+(y.RESOLUTION*1.9)) ))),];
     if(length(finalize)>2) trough=mean(finalize[which(finalize[,2]==min(finalize[,2])),1]) else trough=finalize[1];
     RPE.POSITION.DARK[x,y]=trough;
     if(y==1) abline(v=trough,col="red",lwd=2)};
  ##
  OLM.C=round(mean(which(markings.profile[,2]==254)));
           ## FOR MICE:
           ## the operator should be within a couple pixels (really, within 1, but the registration step may have been imperfect)
           ## less confusion if we grab a bit extra interior... data are sharp-enough
           ## that the usual polynomial function to seach for the peak 
           ## then just puts us within 1 pixel of the actual peak, and we use the local max. 
           ## 
  ## look back ~8 microns and forward ~8, and find the position of lowest slope (usually a small peak) (based non 3rd order polnomial)
  OLM.CHECK=markings.profile[seq(OLM.C-5,OLM.C+5,1),1];
  OLM.dark.CHECK=dark.profile[seq(OLM.C-5,OLM.C+5,1),];
  matlines(OLM.CHECK,OLM.dark.CHECK[,1],col="red");
  for(y in 1:ncol(OLM.dark.CHECK))
    {TEST=cbind(OLM.CHECK,OLM.dark.CHECK[,y]);
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2)+I(TEST[,1]^3))$coef;
     DERIV=seq(TEST[1,1],TEST[11,1],(TEST[11,1]-TEST[1,1])/30);
     DERIV=cbind(DERIV, (MODEL[4]*3*DERIV^2+MODEL[3]*2*DERIV+MODEL[2]));
     # 
     DERIV=cbind(DERIV,sign(DERIV[,2]));
     DERIV=cbind(DERIV,c(DERIV[2:nrow(DERIV),3],DERIV[nrow(DERIV),3]));
     DERIV[,4]<-ifelse(DERIV[,3]>DERIV[,4],1,1000);
     DERIV[,4]=DERIV[,4]+c(seq(15,0,-1),seq(1,15,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
     near=which(DERIV[,4]==min(na.rm=T,DERIV[,4]))[1]; 
     if(near==nrow(DERIV)) near=nrow(DERIV)-1;
     DERIV2=as.matrix(DERIV[c(near,near+1),c(1,2,3)]);
     if(sum(DERIV2[,3])==0) peakISH=summary(lm(DERIV2[,1]~DERIV2[,2]))$coef[1,1] else peakISH=DERIV2[1,1];
     ## look within 1.9 pixels of presumed peak
     finalize=TEST[which( (TEST[,1]>(peakISH-(y.RESOLUTION*1.9)) & (TEST[,1]<(peakISH+(y.RESOLUTION*1.9)) ))),];
     if(length(finalize)>2) peak=mean(finalize[which(finalize[,2]==max(finalize[,2])),1]) else peak=finalize[1];
     if(y==1) {matlines(TEST[,1],(MODEL[4]*TEST[,1]^3+MODEL[3]*TEST[,1]^2+MODEL[2]*TEST[,1]+MODEL[1]),lwd=3);abline(v=peak,col="red",lwd=2)};
     OLM.POSITION.DARK[x,y]=peak;
     ##
     ## now look ~12 to ~48 microns interior (six to 25 pix interior)
     ##
     startONL=which(abs(markings.profile[,1]-peak)==min(abs(markings.profile[,1]-peak)))[1];
     TYPICAL.ONL.SIGNAL.DARK[x,y]=mean(na.rm=TRUE, dark.profile[(startONL-6):(startONL-25),y])};
  ##  
     ## retype that for light
  ## look back ~8 microns and forward ~10, and find the local (based non 2rd order polnomial)
  RPE.light.CHECK=light.profile[seq(RPE.C-6,RPE.C+6,1),];
  RPE.light.CHECK.2=light.profile[seq(RPE.C-12,RPE.C+6,1),];
  #matlines(RPE.CHECK,RPE.light.CHECK[,1],col="red")
  ## 
           ## FOR MICE:
           ## the way the line is drawn, a local minimum exterior to the RPE line
           ## will consitently be basement membrane, and I can cut out anything exterior to that. 
  for(y in 1:ncol(RPE.light.CHECK))
    {TEST=cbind(RPE.CHECK,RPE.light.CHECK[,y]);
     ## I know the line is drawn at whatever's at the 7th row of TEST
     cut=which(TEST[1:7,2]==min(na.rm=T,TEST[1:7,2]));
     ## assuming the min value in that range is NOT the one at position, discard stuff below the basement membrane
     if((cut>1)&(cut<7)) TEST=TEST[(cut-1):nrow(TEST),];
     ##
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2))$coef;
     peak=MODEL[2]/(MODEL[3]*(-2)); ##gives vertex x position;
     ## if the parabola is in the wrong direction, just use the location of the highest value.
     if(MODEL[3]>0) peak=TEST[which(TEST[,2]==max(na.rm=T,TEST[,2])),1];
     ## if peak is more than 2 microns from the swath identified, just use upper or lower bound
     if(peak>max(na.rm=T,TEST[,1]+2)) peak=max(na.rm=T,TEST[,1]+(y.RESOLUTION/2));
     if(peak<min(na.rm=T,TEST[,1]-2)) peak=min(na.rm=T,TEST[,1]-(y.RESOLUTION/2));
     #if(y==1) {matlines(TEST[,1],(MODEL[3]*TEST[,1]^2+MODEL[2]*TEST[,1]+MODEL[1]),lwd=3);abline(v=peak,col="red",lwd=2)};
     peak.value=MODEL[3]*peak^2+MODEL[2]*peak+MODEL[1];
     #RPE.POSITION.LIGHT[x,y]=peak;
     RPE.SIGNAL.LIGHT[x,y]=peak.value[1];
     TEST=cbind(RPE.CHECK.2,RPE.light.CHECK.2[,y]);
     # grab everything behind the peak
     testROW=which(abs(TEST[,1]-peak)==min(abs(TEST[,1]-peak))); 
     testROW=testROW[length(testROW)]; 
     if(testROW<4) testROW=4;
     TEST=TEST[1:testROW,];
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2)+I(TEST[,1]^3))$coef;
     DERIV=seq(TEST[1,1],TEST[nrow(TEST),1],(TEST[nrow(TEST),1]-TEST[1,1])/30);
     DERIV=cbind(DERIV, (MODEL[4]*3*DERIV^2+MODEL[3]*2*DERIV+MODEL[2]));
     # 
     DERIV=cbind(DERIV,sign(DERIV[,2]));
     DERIV=cbind(DERIV,c(DERIV[2:nrow(DERIV),3],DERIV[nrow(DERIV),3]));
     DERIV[,4]<-ifelse(DERIV[,3]<DERIV[,4],1,1000);
     DERIV[,4]=DERIV[,4]+c(seq(30,0,-1));
     near=which(DERIV[,4]==min(na.rm=T,DERIV[,4]))[1]; 
     if(near==nrow(DERIV)) near=nrow(DERIV)-1;
     DERIV2=as.matrix(DERIV[c(near,near+1),c(1,2,3)]);
     if(sum(DERIV2[,3])==0) troughISH=summary(lm(DERIV2[,1]~DERIV2[,2]))$coef[1,1] else troughISH=DERIV2[1,1];
     ## look within 1.9 pixels of presumed peak
     finalize=TEST[which( (TEST[,1]>(troughISH-(y.RESOLUTION*1.9)) & (TEST[,1]<(troughISH+(y.RESOLUTION*1.9)) ))),];
     if(length(finalize)>2) trough=mean(finalize[which(finalize[,2]==min(finalize[,2])),1]) else trough=finalize[1];
     RPE.POSITION.LIGHT[x,y]=trough};
  ##
  OLM.light.CHECK=light.profile[seq(OLM.C-5,OLM.C+5,1),];
  #matlines(OLM.CHECK,OLM.dark.CHECK[,1],col="red");
  for(y in 1:ncol(OLM.dark.CHECK))
    {TEST=cbind(OLM.CHECK,OLM.dark.CHECK[,y]);
     MODEL=lm(TEST[,2] ~ TEST[,1]+I(TEST[,1]^2)+I(TEST[,1]^3))$coef;
     DERIV=seq(TEST[1,1],TEST[11,1],(TEST[11,1]-TEST[1,1])/30);
     DERIV=cbind(DERIV, (MODEL[4]*3*DERIV^2+MODEL[3]*2*DERIV+MODEL[2]));
     # 
     DERIV=cbind(DERIV,sign(DERIV[,2]));
     DERIV=cbind(DERIV,c(DERIV[2:nrow(DERIV),3],DERIV[nrow(DERIV),3]));
     DERIV[,4]<-ifelse(DERIV[,3]>DERIV[,4],1,1000);
     DERIV[,4]=DERIV[,4]+c(seq(15,0,-1),seq(1,15,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
     near=which(DERIV[,4]==min(na.rm=T,DERIV[,4]))[1]; 
     if(near==nrow(DERIV)) near=nrow(DERIV)-1;
     DERIV2=as.matrix(DERIV[c(near,near+1),c(1,2,3)]);
     if(sum(DERIV2[,3])==0) peakISH=summary(lm(DERIV2[,1]~DERIV2[,2]))$coef[1,1] else peakISH=DERIV2[1,1];
     ## look within 1.9 pixels of presumed peak
     finalize=TEST[which( (TEST[,1]>(peakISH-(y.RESOLUTION*1.9)) & (TEST[,1]<(peakISH+(y.RESOLUTION*1.9)) ))),];
     if(length(finalize)>2) peak=mean(finalize[which(finalize[,2]==max(finalize[,2])),1]) else peak=finalize[1];
     #if(y==1) {matlines(TEST[,1],(MODEL[4]*TEST[,1]^3+MODEL[3]*TEST[,1]^2+MODEL[2]*TEST[,1]+MODEL[1]),lwd=3);abline(v=peak,col="red",lwd=2)};
     OLM.POSITION.LIGHT[x,y]=peak;
     ##
     ## now look ~12 to ~48 microns interior (six to 25 pix interior)
     ##
     startONL=which(abs(markings.profile[,1]-peak)==min(abs(markings.profile[,1]-peak)))[1];
     TYPICAL.ONL.SIGNAL.LIGHT[x,y]=mean(na.rm=TRUE, light.profile[(startONL-6):(startONL-25),y])};
  ## 
  ## 
  ## CONTINUE ONLY IF WE HAVE ALL THE OTHER BORDERS PRESENT
  CONTINUE=1
  if(length(which(markings.profile[,2]==249))==0) CONTINUE=0;
  if(length(which(markings.profile[,2]==250))==0) CONTINUE=0;
  if(length(which(markings.profile[,2]==252))==0) CONTINUE=0;
  if(length(which(markings.profile[,2]==253))==0) CONTINUE=0;
  if(CONTINUE==1)
   {ONL.OPL.C=round(mean(which(markings.profile[,2]==253)));
    ## look back ~24 microns and forward ~24, and use the nearest-to-marking half-height to get the right position
    ONL.OPL.CHECK=markings.profile[seq(ONL.OPL.C-12,ONL.OPL.C+12,1),1];
    ONL.OPL.dark.CHECK=dark.profile[seq(ONL.OPL.C-12,ONL.OPL.C+12,1),];
    for(y in 1:ncol(ONL.OPL.dark.CHECK))
      {TEST=cbind(ONL.OPL.CHECK,ONL.OPL.dark.CHECK[,y]);
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       TEST[,5]=TEST[,5]+c(seq(12,0,-1),seq(1,12,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) {abline(v=border,col="red",lwd=2)};
       ONL.OPL.POSITION.DARK[x,y]=border};
    INL.IPL.C=round(mean(which(markings.profile[,2]==252)));
    ## look back ~24 microns and forward ~24, and use the nearest-to-marking half-height to get the right position
    INL.IPL.CHECK=markings.profile[seq(INL.IPL.C-12,INL.IPL.C+12,1),1];
    INL.IPL.dark.CHECK=dark.profile[seq(INL.IPL.C-12,INL.IPL.C+12,1),];
    for(y in 1:ncol(INL.IPL.dark.CHECK))
      {TEST=cbind(INL.IPL.CHECK,INL.IPL.dark.CHECK[,y]);
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       TEST[,5]=TEST[,5]+c(seq(12,0,-1),seq(1,12,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) {abline(v=border,col="red",lwd=2)};
       INL.IPL.POSITION.DARK[x,y]=border};
    RNFL.GCL.C=round(mean(which(markings.profile[,2]==250)));
    ## look back ~24 microns and forward ~16, and use the nearest-to-marking half-height to get the right position
    RNFL.GCL.CHECK=markings.profile[seq(RNFL.GCL.C-12,RNFL.GCL.C+8,1),1];
    RNFL.GCL.dark.CHECK=dark.profile[seq(RNFL.GCL.C-12,RNFL.GCL.C+8,1),];
    for(y in 1:ncol(RNFL.GCL.dark.CHECK))
      {TEST=cbind(RNFL.GCL.CHECK,RNFL.GCL.dark.CHECK[,y]);
       ## for GCL, toss everything out past the max
       TEST=TEST[1:which(TEST[,2]==max(TEST[3:nrow(TEST),2]))[1],];
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       ## unlike the last two broders, this is biased against more-exterior crossings (since we aleady removed anything past the peak. 
       TEST[,5]=TEST[,5]+seq(nrow(TEST),1,-1);
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) abline(v=border,col="red",lwd=2);
       RNFL.GCL.POSITION.DARK[x,y]=border};
    VITREOUS.RETINA.C=round(mean(which(markings.profile[,2]==249)));
    VITREOUS.RETINA.CHECK=markings.profile[seq(VITREOUS.RETINA.C-12,VITREOUS.RETINA.C+12,1),1];
    VITREOUS.RETINA.dark.CHECK=dark.profile[seq(VITREOUS.RETINA.C-12,VITREOUS.RETINA.C+12,1),];
    for(y in 1:ncol(VITREOUS.RETINA.dark.CHECK))
      {TEST=cbind(VITREOUS.RETINA.CHECK,VITREOUS.RETINA.dark.CHECK[,y]);
       ## for vitreous, toss everything out before the max
       TEST=TEST[which(TEST[,2]==max(TEST[1:(nrow(TEST)-2),2]))[1]:nrow(TEST),];
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       ## unlike the last two broders, this is biased against more-exterior (from the perspective of the eye center) crossings
       TEST[,5]=TEST[,5]+seq(nrow(TEST),1,-1);
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) abline(v=border,col="red",lwd=2);
       VITREOUS.RETINA.POSITION.DARK[x,y]=border};
    ### duplicate for light. 
    ## look back ~24 microns and forward ~12, and find the position of lowest slope (usually a small peak) (based non 3rd order polnomial)
    ONL.OPL.light.CHECK=light.profile[seq(ONL.OPL.C-12,ONL.OPL.C+12,1),];
    for(y in 1:ncol(ONL.OPL.light.CHECK))
      {TEST=cbind(ONL.OPL.CHECK,ONL.OPL.light.CHECK[,y]);
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       TEST[,5]=TEST[,5]+c(seq(12,0,-1),seq(1,12,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) {abline(v=border,col="red",lwd=2)};
       ONL.OPL.POSITION.LIGHT[x,y]=border};
    INL.IPL.light.CHECK=light.profile[seq(INL.IPL.C-12,INL.IPL.C+12,1),];
    for(y in 1:ncol(INL.IPL.light.CHECK))
      {TEST=cbind(INL.IPL.CHECK,INL.IPL.light.CHECK[,y]);
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       TEST[,5]=TEST[,5]+c(seq(12,0,-1),seq(1,12,1));
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) {abline(v=border,col="red",lwd=2)};
       INL.IPL.POSITION.LIGHT[x,y]=border};
    RNFL.GCL.light.CHECK=light.profile[seq(RNFL.GCL.C-12,RNFL.GCL.C+8,1),];
    for(y in 1:ncol(RNFL.GCL.light.CHECK))
      {TEST=cbind(RNFL.GCL.CHECK,RNFL.GCL.light.CHECK[,y]);
       ## for GCL, toss everything out past the max
       TEST=TEST[1:which(TEST[,2]==max(TEST[3:nrow(TEST),2]))[1],];
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       ## unlike the last two broders, this is biased against more-exterior crossings (since we aleady removed anything past the peak. 
       TEST[,5]=TEST[,5]+seq(nrow(TEST),1,-1);
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) abline(v=border,col="red",lwd=2);
       RNFL.GCL.POSITION.LIGHT[x,y]=border};
    VITREOUS.RETINA.light.CHECK=light.profile[seq(VITREOUS.RETINA.C-12,VITREOUS.RETINA.C+12,1),];
    for(y in 1:ncol(VITREOUS.RETINA.light.CHECK))
      {TEST=cbind(VITREOUS.RETINA.CHECK,VITREOUS.RETINA.light.CHECK[,y]);
       ## for vitreous, toss everything out before the max
       TEST=TEST[which(TEST[,2]==max(TEST[1:(nrow(TEST)-2),2]))[1]:nrow(TEST),];
       mid=(max(na.rm=T,TEST[,2])+min(na.rm=T,TEST[,2]))/2;
       TEST=cbind(TEST,TEST[,2]-mid);
       ## it's at the closest-to-center transition from negative to positive
       TEST=cbind(TEST,sign(TEST[,3]));
       TEST=cbind(TEST,c(TEST[2:nrow(TEST),4],TEST[nrow(TEST),4]));
       TEST[,5]<-ifelse(TEST[,5]==TEST[,4],1000,1);
       ## unlike the last two broders, this is biased against more-exterior (from the perspective of the eye center) crossings
       TEST[,5]=TEST[,5]+seq(nrow(TEST),1,-1);
       ## if no answer, the default is where the marking is placed.
       ## once we know where the transition is, pick where it would be 0.
       near=which(TEST[,5]==min(na.rm=T,TEST[,5]))[1]; 
       if(near==nrow(TEST)) near=nrow(TEST)-1;
       TEST2=as.matrix(TEST[c(near,near+1),c(1,3)]);
       border=summary(lm(TEST2[,1]~TEST2[,2]))$coef[1,1];
       if(y==1) abline(v=border,col="red",lwd=2);
       VITREOUS.RETINA.POSITION.LIGHT[x,y]=border}}}


## post.Rdata


x=1
plot(RPE.SIGNAL.LIGHT[,x])
plot(TYPICAL.ONL.SIGNAL.LIGHT[,x])
plot(VITREOUS.RETINA.POSITION.LIGHT[,x],ylim=c(0,max(na.rm=T,VITREOUS.RETINA.POSITION.LIGHT[,1])))
matlines(RNFL.GCL.POSITION.LIGHT[,x],col="red")
matlines(INL.IPL.POSITION.LIGHT[,x],col="blue")
matlines(ONL.OPL.POSITION.LIGHT[,x],col="red")
matlines(OLM.POSITION.LIGHT[,x],col="blue")
matlines(RPE.POSITION.LIGHT[,x],col="red")

x=1
plot(RPE.SIGNAL.DARK[,x])
plot(TYPICAL.ONL.SIGNAL.DARK[,x])
plot(VITREOUS.RETINA.POSITION.DARK[,x],ylim=c(0,max(na.rm=T,VITREOUS.RETINA.POSITION.DARK[,1])))
matlines(RNFL.GCL.POSITION.DARK[,x],col="red")
matlines(INL.IPL.POSITION.DARK[,x],col="blue")
matlines(ONL.OPL.POSITION.DARK[,x],col="red")
matlines(OLM.POSITION.DARK[,x],col="blue")
matlines(OLM.POSITION.LIGHT[,x],col="blue",lwd=2)
matlines(RPE.POSITION.DARK[,x],col="red")
matlines(RPE.POSITION.LIGHT[,x],col="red",lwd=2)



#########

### B. for before adjustment
B.VITREOUS.RETINA.POSITION.LIGHT=VITREOUS.RETINA.POSITION.LIGHT
B.RNFL.GCL.POSITION.LIGHT=RNFL.GCL.POSITION.LIGHT
B.INL.IPL.POSITION.LIGHT=INL.IPL.POSITION.LIGHT
B.ONL.OPL.POSITION.LIGHT=ONL.OPL.POSITION.LIGHT
B.OLM.POSITION.LIGHT=OLM.POSITION.LIGHT
B.RPE.POSITION.LIGHT=RPE.POSITION.LIGHT

B.VITREOUS.RETINA.POSITION.DARK=VITREOUS.RETINA.POSITION.DARK
B.RNFL.GCL.POSITION.DARK=RNFL.GCL.POSITION.DARK
B.INL.IPL.POSITION.DARK=INL.IPL.POSITION.DARK
B.ONL.OPL.POSITION.DARK=ONL.OPL.POSITION.DARK
B.OLM.POSITION.DARK=OLM.POSITION.DARK
B.RPE.POSITION.DARK=RPE.POSITION.DARK

x=1
plot(VITREOUS.RETINA.POSITION.DARK[,x],ylim=c(0,max(na.rm=T,VITREOUS.RETINA.POSITION.DARK[,1])))
matlines(RNFL.GCL.POSITION.DARK[,x],col="red")
matlines(INL.IPL.POSITION.DARK[,x],col="blue")
matlines(ONL.OPL.POSITION.DARK[,x],col="red")
matlines(OLM.POSITION.DARK[,x],col="blue")
matlines(OLM.POSITION.LIGHT[,x],col="blue",lwd=2)
matlines(RPE.POSITION.DARK[,x],col="red")
matlines(RPE.POSITION.LIGHT[,x],col="red",lwd=2)
matlines(OLM.POSITION.DARK[,x]-RPE.POSITION.DARK[,x])
matlines(VITREOUS.RETINA.POSITION.DARK[,x]-RPE.POSITION.DARK[,x])

## how do we know the error is in one layer ID and not another?
## well, probably won't mis-specify multiple in the same way
## first, look for errors in the RPE identification by comparison with the 
## VITREOUS.RETINA.POSITION, the ONL.OPL.POSITION, the OLM.POSITION, and the INL.IPL.POSITION (all are fairly reliable.

z.threshold=1
split=round(dim(VITREOUS.RETINA.POSITION.DARK)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.DARK)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.DARK)[2])
 {LEFT=cbind(VITREOUS.RETINA.POSITION.DARK[1:split,x],
             ONL.OPL.POSITION.DARK[1:split,x],
             INL.IPL.POSITION.DARK[1:split,x],
             OLM.POSITION.DARK[1:split,x],
             RPE.POSITION.DARK[1:split,x]);
  RIGHT=cbind(VITREOUS.RETINA.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              ONL.OPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              INL.IPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              OLM.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              RPE.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,4]<-LEFT[,4]-LEFT[,5];
  LEFT[,3]<-LEFT[,3]-LEFT[,5];
  LEFT[,2]<-LEFT[,2]-LEFT[,5];
  LEFT[,1]<-LEFT[,1]-LEFT[,5];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,5];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,5];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,5];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,5];
  plot(seq(1,split,1),LEFT[,4],ylim=c(0,max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=4);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
## double weight of OLM.POSITION when calculating z score and the estimate for the replacement
## overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5]+overall.z[,5])/5;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
## quadruple weight of OLM.
##  RPE.POSITION.DARK[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4;
  RPE.POSITION.DARK[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4]+GIVEBACK.4[,4])/5;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(0,max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
## double weight of OLM.POSITION when calculating z score and the estimate for the replacement
## overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5]+overall.z[,5])/5;
   ##
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
## quadruple weight of OLM.
##  RPE.POSITION.DARK[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4
  RPE.POSITION.DARK[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4]+GIVEBACK.4[,4])/5}



###
###
### for the remaining steps, I'll raise the z threshold to 2
### now, repeat that process for the other borders

z.threshold=2
split=round(dim(VITREOUS.RETINA.POSITION.DARK)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.DARK)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.DARK)[2])
 {LEFT=cbind(RPE.POSITION.DARK[1:split,x],
             OLM.POSITION.DARK[1:split,x],
             VITREOUS.RETINA.POSITION.DARK[1:split,x],
             ONL.OPL.POSITION.DARK[1:split,x],
             INL.IPL.POSITION.DARK[1:split,x]);
  RIGHT=cbind(RPE.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              OLM.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              VITREOUS.RETINA.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              ONL.OPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              INL.IPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,4]<-LEFT[,4]-LEFT[,5];
  LEFT[,3]<-LEFT[,3]-LEFT[,5];
  LEFT[,2]<-LEFT[,2]-LEFT[,5];
  LEFT[,1]<-LEFT[,1]-LEFT[,5];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,5];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,5];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,5];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,5];
  plot(seq(1,split,1),LEFT[,4],ylim=c(0,max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  INL.IPL.POSITION.DARK[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(0,max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  INL.IPL.POSITION.DARK[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4}

### now, repeat that process for the other borders



###
###
### for the remaining steps, I'll raise the z threshold to 2
### now, repeat that process for the other borders

z.threshold=2
split=round(dim(VITREOUS.RETINA.POSITION.DARK)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.DARK)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.DARK)[2])
 {LEFT=cbind(INL.IPL.POSITION.DARK[1:split,x],
             RPE.POSITION.DARK[1:split,x],
             OLM.POSITION.DARK[1:split,x],
             VITREOUS.RETINA.POSITION.DARK[1:split,x],
             ONL.OPL.POSITION.DARK[1:split,x]);
  RIGHT=cbind(INL.IPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              RPE.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              OLM.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              VITREOUS.RETINA.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              ONL.OPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,4]<-LEFT[,4]-LEFT[,5];
  LEFT[,3]<-LEFT[,3]-LEFT[,5];
  LEFT[,2]<-LEFT[,2]-LEFT[,5];
  LEFT[,1]<-LEFT[,1]-LEFT[,5];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,5];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,5];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,5];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,5];
  plot(seq(1,split,1),LEFT[,4],ylim=c(0,max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  ONL.OPL.POSITION.DARK[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(0,max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  ONL.OPL.POSITION.DARK[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4}



###
###
### for the remaining steps, I'll raise the z threshold to 2
### now, repeat that process for the other borders

z.threshold=2
split=round(dim(VITREOUS.RETINA.POSITION.DARK)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.DARK)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.DARK)[2])
 {LEFT=cbind(ONL.OPL.POSITION.DARK[1:split,x],
             INL.IPL.POSITION.DARK[1:split,x],
             RPE.POSITION.DARK[1:split,x],
             OLM.POSITION.DARK[1:split,x],
             VITREOUS.RETINA.POSITION.DARK[1:split,x]);
  RIGHT=cbind(ONL.OPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              INL.IPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              RPE.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              OLM.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              VITREOUS.RETINA.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,4]<-LEFT[,4]-LEFT[,5];
  LEFT[,3]<-LEFT[,3]-LEFT[,5];
  LEFT[,2]<-LEFT[,2]-LEFT[,5];
  LEFT[,1]<-LEFT[,1]-LEFT[,5];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,5];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,5];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,5];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,5];
  plot(seq(1,split,1),LEFT[,4],ylim=c(0,max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  VITREOUS.RETINA.POSITION.DARK[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(0,max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  VITREOUS.RETINA.POSITION.DARK[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4}


### for the OLM, set z-threshold to 3
z.threshold=3
split=round(dim(VITREOUS.RETINA.POSITION.DARK)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.DARK)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.DARK)[2])
 {LEFT=cbind(VITREOUS.RETINA.POSITION.DARK[1:split,x],
             ONL.OPL.POSITION.DARK[1:split,x],
             INL.IPL.POSITION.DARK[1:split,x],
             RPE.POSITION.DARK[1:split,x],
             OLM.POSITION.DARK[1:split,x]);
  RIGHT=cbind(VITREOUS.RETINA.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              ONL.OPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              INL.IPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              RPE.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              OLM.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,4]<-LEFT[,4]-LEFT[,5];
  LEFT[,3]<-LEFT[,3]-LEFT[,5];
  LEFT[,2]<-LEFT[,2]-LEFT[,5];
  LEFT[,1]<-LEFT[,1]-LEFT[,5];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,5];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,5];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,5];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,5];
  plot(seq(1,split,1),LEFT[,4],ylim=c(0,max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  OLM.POSITION.DARK[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(0,max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  OLM.POSITION.DARK[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4}




### great!  now I can try and work on the least-stable layer, RNFL.GCL.POSITION


### z.threshold back down to 1.5
z.threshold=1.5
split=round(dim(VITREOUS.RETINA.POSITION.DARK)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.DARK)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.DARK)[2])
 {LEFT=cbind(VITREOUS.RETINA.POSITION.DARK[1:split,x],
             ONL.OPL.POSITION.DARK[1:split,x],
             INL.IPL.POSITION.DARK[1:split,x],
             RPE.POSITION.DARK[1:split,x],
             OLM.POSITION.DARK[1:split,x],
             RNFL.GCL.POSITION.DARK[1:split,x]);
  RIGHT=cbind(VITREOUS.RETINA.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              ONL.OPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              INL.IPL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              RPE.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              OLM.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x],
              RNFL.GCL.POSITION.DARK[split:dim(VITREOUS.RETINA.POSITION.DARK)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,5]<-LEFT[,5]-LEFT[,6];
  LEFT[,4]<-LEFT[,4]-LEFT[,6];
  LEFT[,3]<-LEFT[,3]-LEFT[,6];
  LEFT[,2]<-LEFT[,2]-LEFT[,6];
  LEFT[,1]<-LEFT[,1]-LEFT[,6];
  RIGHT[,5]<-RIGHT[,5]-RIGHT[,6];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,6];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,6];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,6];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,6];
  plot(seq(1,split,1),LEFT[,4],ylim=c(min(na.rm=T,LEFT[,4]),max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,5]);
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  #
  FIVE=cbind(seq(1,split,1), LEFT[,5]);
  FIVE=FIVE[which(!(is.na(FIVE[,2]))),];
  MODEL.5=smooth.spline(FIVE[,1],FIVE[,2],df=5);
  RESULT.5=as.numeric(predict(MODEL.5,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.5,col="red");
  DEVIATION.5=LEFT[,5]-RESULT.5;
  z.DEVIATION.5=DEVIATION.5;
  z.DEVIATION.5=(z.DEVIATION.5-mean(na.rm=T,DEVIATION.5))/sd(na.rm=T,DEVIATION.5);
  #
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.5,z.DEVIATION.1);
  overall.z[,7]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5]+overall.z[,6])/5;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  FIVE=cbind(seq(1,split,1), LEFT[,5]);
  #FIVE=FIVE[which(!(is.na(FIVE[,2]))),];
  MODEL.5=smooth.spline(FIVE[which(overall.z[,ncol(overall.z)]==0),1],FIVE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.5=as.numeric(predict(MODEL.5,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,6]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,6]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,6]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,6]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  GIVEBACK.5=cbind(LEFT[,5],RESULT.5);
  GIVEBACK.5[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.5[,2],GIVEBACK.5[,1]);
  GIVEBACK.5=cbind(GIVEBACK.5,LEFT.ORIG[,5],LEFT[,6]);
  GIVEBACK.5[,4]<-(GIVEBACK.5[,3]-GIVEBACK.5[,1]);
  ###
  ### and give it back!
  RNFL.GCL.POSITION.DARK[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4]+GIVEBACK.5[,4])/5;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(min(na.rm=T,RIGHT[,4]),max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,5]);
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  #
  FIVE=cbind(seq(split,end,1), RIGHT[,5]);
  FIVE=FIVE[which(!(is.na(FIVE[,2]))),];
  MODEL.5=smooth.spline(FIVE[,1],FIVE[,2],df=5);
  RESULT.5=as.numeric(predict(MODEL.5,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.5,col="red");
  DEVIATION.5=RIGHT[,5]-RESULT.5;
  z.DEVIATION.5=DEVIATION.5;
  z.DEVIATION.5=(z.DEVIATION.5-mean(na.rm=T,DEVIATION.5))/sd(na.rm=T,DEVIATION.5);
  #
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.5,z.DEVIATION.1);
  overall.z[,7]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5]+overall.z[,6])/5;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  FIVE=cbind(seq(split,end,1), RIGHT[,5]);
  #FIVE=FIVE[which(!(is.na(FIVE[,2]))),];
  MODEL.5=smooth.spline(FIVE[which(overall.z[,ncol(overall.z)]==0),1],FIVE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.5=as.numeric(predict(MODEL.5,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,6]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.DARK
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,6]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,6]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,6]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  GIVEBACK.5=cbind(RIGHT[,5],RESULT.5);
  GIVEBACK.5[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.5[,2],GIVEBACK.5[,1]);
  GIVEBACK.5=cbind(GIVEBACK.5,RIGHT.ORIG[,5],RIGHT[,6]);
  GIVEBACK.5[,4]<-(GIVEBACK.5[,3]-GIVEBACK.5[,1]);
  ###
  ### and give it back!
  RNFL.GCL.POSITION.DARK[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4]+GIVEBACK.5[,4])/5}




###



### now, repeat all that for LIGHT


## how do we know the error is in one layer ID and not another?
## well, probably won't mis-specify multiple in the same way
## first, look for errors in the RPE identification by comparison with the 
## VITREOUS.RETINA.POSITION, the ONL.OPL.POSITION, the OLM.POSITION, and the INL.IPL.POSITION (all are fairly reliable.

z.threshold=1
split=round(dim(VITREOUS.RETINA.POSITION.LIGHT)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.LIGHT)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.LIGHT)[2])
 {LEFT=cbind(VITREOUS.RETINA.POSITION.LIGHT[1:split,x],
             ONL.OPL.POSITION.LIGHT[1:split,x],
             INL.IPL.POSITION.LIGHT[1:split,x],
             OLM.POSITION.LIGHT[1:split,x],
             RPE.POSITION.LIGHT[1:split,x]);
  RIGHT=cbind(VITREOUS.RETINA.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              ONL.OPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              INL.IPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              OLM.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              RPE.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,4]<-LEFT[,4]-LEFT[,5];
  LEFT[,3]<-LEFT[,3]-LEFT[,5];
  LEFT[,2]<-LEFT[,2]-LEFT[,5];
  LEFT[,1]<-LEFT[,1]-LEFT[,5];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,5];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,5];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,5];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,5];
  plot(seq(1,split,1),LEFT[,4],ylim=c(0,max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=4);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
## double weight of OLM.POSITION when calculating z score and the estimate for the replacement
## overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5]+overall.z[,5])/5;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
## quadruple weight of OLM.
##  RPE.POSITION.LIGHT[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4;
  RPE.POSITION.LIGHT[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4]+GIVEBACK.4[,4])/5;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(0,max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
## double weight of OLM.POSITION when calculating z score and the estimate for the replacement
## overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5]+overall.z[,5])/5;
   ##
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
## quadruple weight of OLM.
##  RPE.POSITION.LIGHT[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4
  RPE.POSITION.LIGHT[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4]+GIVEBACK.4[,4])/5}



###
###
### for the remaining steps, I'll raise the z threshold to 2
### now, repeat that process for the other borders

z.threshold=2
split=round(dim(VITREOUS.RETINA.POSITION.LIGHT)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.LIGHT)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.LIGHT)[2])
 {LEFT=cbind(RPE.POSITION.LIGHT[1:split,x],
             OLM.POSITION.LIGHT[1:split,x],
             VITREOUS.RETINA.POSITION.LIGHT[1:split,x],
             ONL.OPL.POSITION.LIGHT[1:split,x],
             INL.IPL.POSITION.LIGHT[1:split,x]);
  RIGHT=cbind(RPE.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              OLM.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              VITREOUS.RETINA.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              ONL.OPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              INL.IPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,4]<-LEFT[,4]-LEFT[,5];
  LEFT[,3]<-LEFT[,3]-LEFT[,5];
  LEFT[,2]<-LEFT[,2]-LEFT[,5];
  LEFT[,1]<-LEFT[,1]-LEFT[,5];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,5];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,5];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,5];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,5];
  plot(seq(1,split,1),LEFT[,4],ylim=c(0,max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  INL.IPL.POSITION.LIGHT[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(0,max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  INL.IPL.POSITION.LIGHT[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4}

### now, repeat that process for the other borders



###
###
### for the remaining steps, I'll raise the z threshold to 2
### now, repeat that process for the other borders

z.threshold=2
split=round(dim(VITREOUS.RETINA.POSITION.LIGHT)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.LIGHT)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.LIGHT)[2])
 {LEFT=cbind(INL.IPL.POSITION.LIGHT[1:split,x],
             RPE.POSITION.LIGHT[1:split,x],
             OLM.POSITION.LIGHT[1:split,x],
             VITREOUS.RETINA.POSITION.LIGHT[1:split,x],
             ONL.OPL.POSITION.LIGHT[1:split,x]);
  RIGHT=cbind(INL.IPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              RPE.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              OLM.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              VITREOUS.RETINA.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              ONL.OPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,4]<-LEFT[,4]-LEFT[,5];
  LEFT[,3]<-LEFT[,3]-LEFT[,5];
  LEFT[,2]<-LEFT[,2]-LEFT[,5];
  LEFT[,1]<-LEFT[,1]-LEFT[,5];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,5];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,5];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,5];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,5];
  plot(seq(1,split,1),LEFT[,4],ylim=c(0,max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  ONL.OPL.POSITION.LIGHT[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(0,max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  ONL.OPL.POSITION.LIGHT[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4}



###
###
### for the remaining steps, I'll raise the z threshold to 2
### now, repeat that process for the other borders

z.threshold=2
split=round(dim(VITREOUS.RETINA.POSITION.LIGHT)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.LIGHT)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.LIGHT)[2])
 {LEFT=cbind(ONL.OPL.POSITION.LIGHT[1:split,x],
             INL.IPL.POSITION.LIGHT[1:split,x],
             RPE.POSITION.LIGHT[1:split,x],
             OLM.POSITION.LIGHT[1:split,x],
             VITREOUS.RETINA.POSITION.LIGHT[1:split,x]);
  RIGHT=cbind(ONL.OPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              INL.IPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              RPE.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              OLM.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              VITREOUS.RETINA.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,4]<-LEFT[,4]-LEFT[,5];
  LEFT[,3]<-LEFT[,3]-LEFT[,5];
  LEFT[,2]<-LEFT[,2]-LEFT[,5];
  LEFT[,1]<-LEFT[,1]-LEFT[,5];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,5];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,5];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,5];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,5];
  plot(seq(1,split,1),LEFT[,4],ylim=c(0,max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  VITREOUS.RETINA.POSITION.LIGHT[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(0,max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  VITREOUS.RETINA.POSITION.LIGHT[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4}


### for the OLM, set z-threshold to 3
z.threshold=3
split=round(dim(VITREOUS.RETINA.POSITION.LIGHT)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.LIGHT)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.LIGHT)[2])
 {LEFT=cbind(VITREOUS.RETINA.POSITION.LIGHT[1:split,x],
             ONL.OPL.POSITION.LIGHT[1:split,x],
             INL.IPL.POSITION.LIGHT[1:split,x],
             RPE.POSITION.LIGHT[1:split,x],
             OLM.POSITION.LIGHT[1:split,x]);
  RIGHT=cbind(VITREOUS.RETINA.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              ONL.OPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              INL.IPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              RPE.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              OLM.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,4]<-LEFT[,4]-LEFT[,5];
  LEFT[,3]<-LEFT[,3]-LEFT[,5];
  LEFT[,2]<-LEFT[,2]-LEFT[,5];
  LEFT[,1]<-LEFT[,1]-LEFT[,5];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,5];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,5];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,5];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,5];
  plot(seq(1,split,1),LEFT[,4],ylim=c(0,max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  OLM.POSITION.LIGHT[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(0,max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.1);
  overall.z[,6]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5])/4;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,5]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,5]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,5]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,5]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  ###
  ### and give it back!
  OLM.POSITION.LIGHT[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4])/4}




### great!  now I can try and work on the least-stable layer, RNFL.GCL.POSITION


### z.threshold back down to 1.5
z.threshold=1.5
split=round(dim(VITREOUS.RETINA.POSITION.LIGHT)[1]/2)
end=dim(VITREOUS.RETINA.POSITION.LIGHT)[1]
for(x in 1:dim(VITREOUS.RETINA.POSITION.LIGHT)[2])
 {LEFT=cbind(VITREOUS.RETINA.POSITION.LIGHT[1:split,x],
             ONL.OPL.POSITION.LIGHT[1:split,x],
             INL.IPL.POSITION.LIGHT[1:split,x],
             RPE.POSITION.LIGHT[1:split,x],
             OLM.POSITION.LIGHT[1:split,x],
             RNFL.GCL.POSITION.LIGHT[1:split,x]);
  RIGHT=cbind(VITREOUS.RETINA.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              ONL.OPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              INL.IPL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              RPE.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              OLM.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x],
              RNFL.GCL.POSITION.LIGHT[split:dim(VITREOUS.RETINA.POSITION.LIGHT)[1],x]);
  LEFT.ORIG=LEFT;
  RIGHT.ORIG=RIGHT;
  LEFT[,5]<-LEFT[,5]-LEFT[,6];
  LEFT[,4]<-LEFT[,4]-LEFT[,6];
  LEFT[,3]<-LEFT[,3]-LEFT[,6];
  LEFT[,2]<-LEFT[,2]-LEFT[,6];
  LEFT[,1]<-LEFT[,1]-LEFT[,6];
  RIGHT[,5]<-RIGHT[,5]-RIGHT[,6];
  RIGHT[,4]<-RIGHT[,4]-RIGHT[,6];
  RIGHT[,3]<-RIGHT[,3]-RIGHT[,6];
  RIGHT[,2]<-RIGHT[,2]-RIGHT[,6];
  RIGHT[,1]<-RIGHT[,1]-RIGHT[,6];
  plot(seq(1,split,1),LEFT[,4],ylim=c(min(na.rm=T,LEFT[,4]),max(na.rm=T,LEFT[,1])));
  matlines(seq(1,split,1),LEFT[,5]);
  matlines(seq(1,split,1),LEFT[,4]);
  matlines(seq(1,split,1),LEFT[,3]);
  matlines(seq(1,split,1),LEFT[,2]);
  matlines(seq(1,split,1),LEFT[,1]);
  #
  FIVE=cbind(seq(1,split,1), LEFT[,5]);
  FIVE=FIVE[which(!(is.na(FIVE[,2]))),];
  MODEL.5=smooth.spline(FIVE[,1],FIVE[,2],df=5);
  RESULT.5=as.numeric(predict(MODEL.5,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.5,col="red");
  DEVIATION.5=LEFT[,5]-RESULT.5;
  z.DEVIATION.5=DEVIATION.5;
  z.DEVIATION.5=(z.DEVIATION.5-mean(na.rm=T,DEVIATION.5))/sd(na.rm=T,DEVIATION.5);
  #
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.4,col="red");
  DEVIATION.4=LEFT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.3,col="red");
  DEVIATION.3=LEFT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.2,col="red");
  DEVIATION.2=LEFT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  matlines(seq(1,split,1),RESULT.1,col="red");
  DEVIATION.1=LEFT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(1,split,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.5,z.DEVIATION.1);
  overall.z[,7]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5]+overall.z[,6])/5;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(1,split,1), LEFT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(1,split,1))$y);
  TWO=cbind(seq(1,split,1), LEFT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(1,split,1))$y);
  THREE=cbind(seq(1,split,1), LEFT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(1,split,1))$y);
  FOUR=cbind(seq(1,split,1), LEFT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(1,split,1))$y);
  FIVE=cbind(seq(1,split,1), LEFT[,5]);
  #FIVE=FIVE[which(!(is.na(FIVE[,2]))),];
  MODEL.5=smooth.spline(FIVE[which(overall.z[,ncol(overall.z)]==0),1],FIVE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.5=as.numeric(predict(MODEL.5,seq(1,split,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(LEFT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,LEFT.ORIG[,1],LEFT[,6]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(LEFT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,LEFT.ORIG[,2],LEFT[,6]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(LEFT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,LEFT.ORIG[,3],LEFT[,6]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(LEFT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,LEFT.ORIG[,4],LEFT[,6]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  GIVEBACK.5=cbind(LEFT[,5],RESULT.5);
  GIVEBACK.5[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.5[,2],GIVEBACK.5[,1]);
  GIVEBACK.5=cbind(GIVEBACK.5,LEFT.ORIG[,5],LEFT[,6]);
  GIVEBACK.5[,4]<-(GIVEBACK.5[,3]-GIVEBACK.5[,1]);
  ###
  ### and give it back!
  RNFL.GCL.POSITION.LIGHT[1:split,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4]+GIVEBACK.5[,4])/5;
  ###
  ### now the right side:
  plot(seq(split,end,1),RIGHT[,4],ylim=c(min(na.rm=T,RIGHT[,4]),max(na.rm=T,RIGHT[,1])));
  matlines(seq(split,end,1),RIGHT[,5]);
  matlines(seq(split,end,1),RIGHT[,4]);
  matlines(seq(split,end,1),RIGHT[,3]);
  matlines(seq(split,end,1),RIGHT[,2]);
  matlines(seq(split,end,1),RIGHT[,1]);
  #
  FIVE=cbind(seq(split,end,1), RIGHT[,5]);
  FIVE=FIVE[which(!(is.na(FIVE[,2]))),];
  MODEL.5=smooth.spline(FIVE[,1],FIVE[,2],df=5);
  RESULT.5=as.numeric(predict(MODEL.5,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.5,col="red");
  DEVIATION.5=RIGHT[,5]-RESULT.5;
  z.DEVIATION.5=DEVIATION.5;
  z.DEVIATION.5=(z.DEVIATION.5-mean(na.rm=T,DEVIATION.5))/sd(na.rm=T,DEVIATION.5);
  #
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.4,col="red");
  DEVIATION.4=RIGHT[,4]-RESULT.4
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  #
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[,1],THREE[,2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.3,col="red");
  DEVIATION.3=RIGHT[,3]-RESULT.3
  z.DEVIATION.3=DEVIATION.3;
  z.DEVIATION.3=(z.DEVIATION.3-mean(na.rm=T,DEVIATION.3))/sd(na.rm=T,DEVIATION.3);
  #
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[,1],TWO[,2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.2,col="red");
  DEVIATION.2=RIGHT[,2]-RESULT.2
  z.DEVIATION.2=DEVIATION.2;
  z.DEVIATION.2=(z.DEVIATION.2-mean(na.rm=T,DEVIATION.2))/sd(na.rm=T,DEVIATION.2);
  #
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[,1],ONE[,2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  matlines(seq(split,end,1),RESULT.1,col="red");
  DEVIATION.1=RIGHT[,1]-RESULT.1
  z.DEVIATION.1=DEVIATION.1;
  z.DEVIATION.1=(z.DEVIATION.1-mean(na.rm=T,DEVIATION.1))/sd(na.rm=T,DEVIATION.1);
  # 
  # 
  overall.z=cbind(seq(split,end,1),z.DEVIATION.1,z.DEVIATION.2,z.DEVIATION.3,z.DEVIATION.4,z.DEVIATION.5,z.DEVIATION.1);
  overall.z[,7]<-abs(overall.z[,2]+overall.z[,3]+overall.z[,4]+overall.z[,5]+overall.z[,6])/5;
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,ncol(overall.z)]<-ifelse(overall.z[,(ncol(overall.z)-1)]>z.threshold,1,0);
  ##
  ##
  ## ### now, re-calculate RESULT without that above-threshold data to see what we should replace it with:
  ##
  ONE=cbind(seq(split,end,1), RIGHT[,1]);
  #ONE=ONE[which(!(is.na(ONE[,2]))),];
  MODEL.1=smooth.spline(ONE[which(overall.z[,ncol(overall.z)]==0),1],ONE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.1=as.numeric(predict(MODEL.1,seq(split,end,1))$y);
  TWO=cbind(seq(split,end,1), RIGHT[,2]);
  #TWO=TWO[which(!(is.na(TWO[,2]))),];
  MODEL.2=smooth.spline(TWO[which(overall.z[,ncol(overall.z)]==0),1],TWO[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.2=as.numeric(predict(MODEL.2,seq(split,end,1))$y);
  THREE=cbind(seq(split,end,1), RIGHT[,3]);
  #THREE=THREE[which(!(is.na(THREE[,2]))),];
  MODEL.3=smooth.spline(THREE[which(overall.z[,ncol(overall.z)]==0),1],THREE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.3=as.numeric(predict(MODEL.3,seq(split,end,1))$y);
  FOUR=cbind(seq(split,end,1), RIGHT[,4]);
  #FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[which(overall.z[,ncol(overall.z)]==0),1],FOUR[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,seq(split,end,1))$y);
  FIVE=cbind(seq(split,end,1), RIGHT[,5]);
  #FIVE=FIVE[which(!(is.na(FIVE[,2]))),];
  MODEL.5=smooth.spline(FIVE[which(overall.z[,ncol(overall.z)]==0),1],FIVE[which(overall.z[,ncol(overall.z)]==0),2],df=5);
  RESULT.5=as.numeric(predict(MODEL.5,seq(split,end,1))$y);
  ##
  ## ### 
  ##
  ##
  GIVEBACK.1=cbind(RIGHT[,1],RESULT.1);
  GIVEBACK.1[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.1[,2],GIVEBACK.1[,1]);
  GIVEBACK.1=cbind(GIVEBACK.1,RIGHT.ORIG[,1],RIGHT[,6]);
  GIVEBACK.1[,4]<-(GIVEBACK.1[,3]-GIVEBACK.1[,1]);
  ## now, GIVEBACK.1[,4] is one of four estimates we can generate for replacing RPE.POSITION.LIGHT
  ## let's make the other three
  GIVEBACK.2=cbind(RIGHT[,2],RESULT.2);
  GIVEBACK.2[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.2[,2],GIVEBACK.2[,1]);
  GIVEBACK.2=cbind(GIVEBACK.2,RIGHT.ORIG[,2],RIGHT[,6]);
  GIVEBACK.2[,4]<-(GIVEBACK.2[,3]-GIVEBACK.2[,1]);
  GIVEBACK.3=cbind(RIGHT[,3],RESULT.3);
  GIVEBACK.3[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.3[,2],GIVEBACK.3[,1]);
  GIVEBACK.3=cbind(GIVEBACK.3,RIGHT.ORIG[,3],RIGHT[,6]);
  GIVEBACK.3[,4]<-(GIVEBACK.3[,3]-GIVEBACK.3[,1]);
  GIVEBACK.4=cbind(RIGHT[,4],RESULT.4);
  GIVEBACK.4[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.4[,2],GIVEBACK.4[,1]);
  GIVEBACK.4=cbind(GIVEBACK.4,RIGHT.ORIG[,4],RIGHT[,6]);
  GIVEBACK.4[,4]<-(GIVEBACK.4[,3]-GIVEBACK.4[,1]);
  GIVEBACK.5=cbind(RIGHT[,5],RESULT.5);
  GIVEBACK.5[,1]<-ifelse(overall.z[,ncol(overall.z)]==1,GIVEBACK.5[,2],GIVEBACK.5[,1]);
  GIVEBACK.5=cbind(GIVEBACK.5,RIGHT.ORIG[,5],RIGHT[,6]);
  GIVEBACK.5[,4]<-(GIVEBACK.5[,3]-GIVEBACK.5[,1]);
  ###
  ### and give it back!
  RNFL.GCL.POSITION.LIGHT[split:end,x]=(GIVEBACK.1[,4]+GIVEBACK.2[,4]+GIVEBACK.3[,4]+GIVEBACK.4[,4]+GIVEBACK.5[,4])/5}







## DEC2019





#########  In this next step,
#########  we use a high-degree-of-freedom spline to smooth out any jitter in the localization.

### R. for revised
R.RPE.SIGNAL.LIGHT=RPE.SIGNAL.LIGHT
R.TYPICAL.ONL.SIGNAL.LIGHT=TYPICAL.ONL.SIGNAL.LIGHT
R.VITREOUS.RETINA.POSITION.LIGHT=VITREOUS.RETINA.POSITION.LIGHT
R.RNFL.GCL.POSITION.LIGHT=RNFL.GCL.POSITION.LIGHT
R.INL.IPL.POSITION.LIGHT=INL.IPL.POSITION.LIGHT
R.ONL.OPL.POSITION.LIGHT=ONL.OPL.POSITION.LIGHT
R.OLM.POSITION.LIGHT=OLM.POSITION.LIGHT
R.RPE.POSITION.LIGHT=RPE.POSITION.LIGHT

R.RPE.SIGNAL.DARK=RPE.SIGNAL.DARK
R.TYPICAL.ONL.SIGNAL.DARK=TYPICAL.ONL.SIGNAL.DARK
R.VITREOUS.RETINA.POSITION.DARK=VITREOUS.RETINA.POSITION.DARK
R.RNFL.GCL.POSITION.DARK=RNFL.GCL.POSITION.DARK
R.INL.IPL.POSITION.DARK=INL.IPL.POSITION.DARK
R.ONL.OPL.POSITION.DARK=ONL.OPL.POSITION.DARK
R.OLM.POSITION.DARK=OLM.POSITION.DARK
R.RPE.POSITION.DARK=RPE.POSITION.DARK

##

# start.move:break.left
# break.right:end.move
for(x in 1:dim(VITREOUS.RETINA.POSITION.LIGHT)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],VITREOUS.RETINA.POSITION.LIGHT[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.VITREOUS.RETINA.POSITION.LIGHT[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],VITREOUS.RETINA.POSITION.LIGHT[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.VITREOUS.RETINA.POSITION.LIGHT[break.right:end.move,x]=RESULT[,2]}

for(x in 1:dim(VITREOUS.RETINA.POSITION.DARK)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],VITREOUS.RETINA.POSITION.DARK[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.VITREOUS.RETINA.POSITION.DARK[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],VITREOUS.RETINA.POSITION.DARK[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.VITREOUS.RETINA.POSITION.DARK[break.right:end.move,x]=RESULT[,2]}




#############################
#############################


# start.move:break.left
# break.right:end.move
for(x in 1:dim(RNFL.GCL.POSITION.LIGHT)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],RNFL.GCL.POSITION.LIGHT[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.RNFL.GCL.POSITION.LIGHT[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],RNFL.GCL.POSITION.LIGHT[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.RNFL.GCL.POSITION.LIGHT[break.right:end.move,x]=RESULT[,2]}

for(x in 1:dim(RNFL.GCL.POSITION.DARK)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],RNFL.GCL.POSITION.DARK[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.RNFL.GCL.POSITION.DARK[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],RNFL.GCL.POSITION.DARK[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.RNFL.GCL.POSITION.DARK[break.right:end.move,x]=RESULT[,2]}

#############################
#############################

# start.move:break.left
# break.right:end.move
for(x in 1:dim(INL.IPL.POSITION.LIGHT)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],INL.IPL.POSITION.LIGHT[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.INL.IPL.POSITION.LIGHT[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],INL.IPL.POSITION.LIGHT[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.INL.IPL.POSITION.LIGHT[break.right:end.move,x]=RESULT[,2]}

for(x in 1:dim(INL.IPL.POSITION.DARK)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],INL.IPL.POSITION.DARK[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.INL.IPL.POSITION.DARK[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],INL.IPL.POSITION.DARK[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.INL.IPL.POSITION.DARK[break.right:end.move,x]=RESULT[,2]}



#############################
#############################


# start.move:break.left
# break.right:end.move
for(x in 1:dim(ONL.OPL.POSITION.LIGHT)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],ONL.OPL.POSITION.LIGHT[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.ONL.OPL.POSITION.LIGHT[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],ONL.OPL.POSITION.LIGHT[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.ONL.OPL.POSITION.LIGHT[break.right:end.move,x]=RESULT[,2]}

for(x in 1:dim(ONL.OPL.POSITION.DARK)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],ONL.OPL.POSITION.DARK[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.ONL.OPL.POSITION.DARK[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],ONL.OPL.POSITION.DARK[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.ONL.OPL.POSITION.DARK[break.right:end.move,x]=RESULT[,2]}


#############################
#############################


# start.move:break.left
# break.right:end.move
for(x in 1:dim(OLM.POSITION.LIGHT)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],OLM.POSITION.LIGHT[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.OLM.POSITION.LIGHT[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],OLM.POSITION.LIGHT[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.OLM.POSITION.LIGHT[break.right:end.move,x]=RESULT[,2]}

for(x in 1:dim(OLM.POSITION.DARK)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],OLM.POSITION.DARK[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.OLM.POSITION.DARK[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],OLM.POSITION.DARK[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.OLM.POSITION.DARK[break.right:end.move,x]=RESULT[,2]}


#############################
#############################

# start.move:break.left
# break.right:end.move
for(x in 1:dim(RPE.POSITION.LIGHT)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],RPE.POSITION.LIGHT[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.RPE.POSITION.LIGHT[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],RPE.POSITION.LIGHT[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.RPE.POSITION.LIGHT[break.right:end.move,x]=RESULT[,2]}

for(x in 1:dim(RPE.POSITION.DARK)[2])
 {DISPLAY=cbind(Retina.Points.cropped[start.move:break.left,4],RPE.POSITION.DARK[start.move:break.left,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[start.move:break.left,4],as.numeric(predict(MODEL,Retina.Points.cropped[start.move:break.left,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.RPE.POSITION.DARK[start.move:break.left,x]=RESULT[,2];
  ##
  DISPLAY=cbind(Retina.Points.cropped[break.right:end.move,4],RPE.POSITION.DARK[break.right:end.move,x]);
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  M=mean(na.rm=T,DISPLAY[,2]);
  SD=sd(na.rm=T,DISPLAY[,2]);
  DISPLAY[which(DISPLAY[,2]<(M-3*SD)),2]=NA;
  DISPLAY[which(DISPLAY[,2]>(M+3*SD)),2]=NA;
  ## anchor a bit...
  SAVE=which(!(is.na(DISPLAY[,2])))[1]-1;
  plot(DISPLAY[,1],DISPLAY[,2]);
  DISPLAY=DISPLAY[which(!(is.na(DISPLAY[,2]))),]
  MODEL=smooth.spline(DISPLAY[,1], DISPLAY[,2],df=11);
  RESULT=cbind(Retina.Points.cropped[break.right:end.move,4],as.numeric(predict(MODEL,Retina.Points.cropped[break.right:end.move,4])$y));
  #RESULT[1:SAVE,2]=NA;
  matlines(RESULT[,1],RESULT[,2],col="red");
  R.RPE.POSITION.DARK[break.right:end.move,x]=RESULT[,2]}


#############################
############################# zz.Rdata DEC2019



##### make a file that holds the info for the thicknesses for just the first image (which, in current use, is usually the only image)
## distance between the two is 
SPECIAL.OUTPUT=as.data.frame(matrix(,(5+length(start.move:break.left)+length(break.right:end.move)),10))
SPECIAL.OUTPUT[6:(5+length(start.move:break.left)),1]="LEFT"
SPECIAL.OUTPUT[(6+length(start.move:break.left)):nrow(SPECIAL.OUTPUT),1]="RIGHT"
SPECIAL.OUTPUT[6:(5+length(start.move:break.left)),2]=Retina.Points.cropped[start.move:break.left,4]
SPECIAL.OUTPUT[(6+length(start.move:break.left)):nrow(SPECIAL.OUTPUT),2]=Retina.Points.cropped[break.right:end.move,4]

SPECIAL.OUTPUT[5,1]="SIDE"
SPECIAL.OUTPUT[5,2]="POSITION_microns"
SPECIAL.OUTPUT[6:(5+length(start.move:break.left)),3]=round((R.VITREOUS.RETINA.POSITION.DARK[start.move:break.left,1]-R.RPE.POSITION.DARK[start.move:break.left,1])*100)/100
SPECIAL.OUTPUT[(6+length(start.move:break.left)):nrow(SPECIAL.OUTPUT),3]=round((R.VITREOUS.RETINA.POSITION.DARK[break.right:end.move,1]-R.RPE.POSITION.DARK[break.right:end.move,1])*100)/100
SPECIAL.OUTPUT[5,3]="WHOLE_RETINA"
#
SPECIAL.OUTPUT[6:(5+length(start.move:break.left)),4]=round((R.VITREOUS.RETINA.POSITION.DARK[start.move:break.left,1]-R.RNFL.GCL.POSITION.DARK[start.move:break.left,1])*100)/100
SPECIAL.OUTPUT[(6+length(start.move:break.left)):nrow(SPECIAL.OUTPUT),4]=round((R.VITREOUS.RETINA.POSITION.DARK[break.right:end.move,1]-R.RNFL.GCL.POSITION.DARK[break.right:end.move,1])*100)/100
SPECIAL.OUTPUT[5,4]="RNFL_THICK"
#
SPECIAL.OUTPUT[6:(5+length(start.move:break.left)),5]=round((R.RNFL.GCL.POSITION.DARK[start.move:break.left,1]-R.INL.IPL.POSITION.DARK[start.move:break.left,1])*100)/100
SPECIAL.OUTPUT[(6+length(start.move:break.left)):nrow(SPECIAL.OUTPUT),5]=round((R.RNFL.GCL.POSITION.DARK[break.right:end.move,1]-R.INL.IPL.POSITION.DARK[break.right:end.move,1])*100)/100
SPECIAL.OUTPUT[5,5]="INL_IPL_border_to_GCL_RNFL_border"
#
SPECIAL.OUTPUT[6:(5+length(start.move:break.left)),6]=round((R.INL.IPL.POSITION.DARK[start.move:break.left,1]-R.ONL.OPL.POSITION.DARK[start.move:break.left,1])*100)/100
SPECIAL.OUTPUT[(6+length(start.move:break.left)):nrow(SPECIAL.OUTPUT),6]=round((R.INL.IPL.POSITION.DARK[break.right:end.move,1]-R.ONL.OPL.POSITION.DARK[break.right:end.move,1])*100)/100
SPECIAL.OUTPUT[5,6]="ONL_OPL_border_to_INL_IPL_border"
#
SPECIAL.OUTPUT[6:(5+length(start.move:break.left)),7]=round((R.ONL.OPL.POSITION.DARK[start.move:break.left,1]-R.OLM.POSITION.DARK[start.move:break.left,1])*100)/100
SPECIAL.OUTPUT[(6+length(start.move:break.left)):nrow(SPECIAL.OUTPUT),7]=round((R.ONL.OPL.POSITION.DARK[break.right:end.move,1]-R.OLM.POSITION.DARK[break.right:end.move,1])*100)/100
SPECIAL.OUTPUT[5,7]="OLM_to_ONL_OPL_border"
#
SPECIAL.OUTPUT[6:(5+length(start.move:break.left)),8]=round((R.OLM.POSITION.DARK[start.move:break.left,1]-R.RPE.POSITION.DARK[start.move:break.left,1])*100)/100
SPECIAL.OUTPUT[(6+length(start.move:break.left)):nrow(SPECIAL.OUTPUT),8]=round((R.OLM.POSITION.DARK[break.right:end.move,1]-R.RPE.POSITION.DARK[break.right:end.move,1])*100)/100
SPECIAL.OUTPUT[5,8]="base_of_RPE_to_OLM"
#
SPECIAL.OUTPUT=SPECIAL.OUTPUT[5:nrow(SPECIAL.OUTPUT),1:8]
write(t(SPECIAL.OUTPUT),ncol=ncol(SPECIAL.OUTPUT),file=paste(FILENAME,"_first-image-only_thickness_details_c.txt",sep=""))




IMAGE.INDEX.DARK=seq(1,dim(DARK.DATA)[3],1)
IMAGE.INDEX.LIGHT=seq(1,dim(LIGHT.DATA)[3],1)

Zz.MEAN.RPE.OML.DIST.DARK=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)
Zz.MEAN.RPE.OML.DIST.LIGHT=cbind(IMAGE.INDEX.LIGHT,IMAGE.INDEX.LIGHT)

for(x in 1:nrow(Zz.MEAN.RPE.OML.DIST.DARK))
  {Zz.MEAN.RPE.OML.DIST.DARK[x,2]=mean(na.rm=TRUE,(R.OLM.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]-R.RPE.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.RPE.OML.DIST.LIGHT))
  {Zz.MEAN.RPE.OML.DIST.LIGHT[x,2]=mean(na.rm=TRUE,(R.OLM.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]-R.RPE.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]))}

mean(Zz.MEAN.RPE.OML.DIST.DARK[,2])
mean(Zz.MEAN.RPE.OML.DIST.LIGHT[,2])
#t.test(Zz.MEAN.RPE.OML.DIST.DARK[,2],Zz.MEAN.RPE.OML.DIST.LIGHT[,2])



Zz.MEAN.RPE.VITREOUS.DIST.DARK=Zz.MEAN.RPE.OML.DIST.DARK
Zz.MEAN.RPE.VITREOUS.DIST.LIGHT=Zz.MEAN.RPE.OML.DIST.LIGHT
for(x in 1:nrow(Zz.MEAN.RPE.VITREOUS.DIST.DARK))
  {Zz.MEAN.RPE.VITREOUS.DIST.DARK[x,2]=mean(na.rm=TRUE,(R.VITREOUS.RETINA.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]-R.RPE.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.RPE.VITREOUS.DIST.LIGHT))
  {Zz.MEAN.RPE.VITREOUS.DIST.LIGHT[x,2]=mean(na.rm=TRUE,(R.VITREOUS.RETINA.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]-R.RPE.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]))}

mean(Zz.MEAN.RPE.VITREOUS.DIST.DARK[,2])
mean(Zz.MEAN.RPE.VITREOUS.DIST.LIGHT[,2])
#t.test(Zz.MEAN.RPE.VITREOUS.DIST.DARK[,2],Zz.MEAN.RPE.VITREOUS.DIST.LIGHT[,2])
###############################
###############################

###############################
###############################

###############################
###############################


###############################
###############################



## assign 6 points from -30microns to 0 [inclusive of lower value and upper, all others are inclusive of upper value only]
## assign 19 points to the RPE-OLM span
## assign 18 points to the OLM-ONL/OPL span
## assign 18 points to the ONL/OPL - INL/IPL span
## assign 25 points to the INL/IPL - RNFL/GCL span
## assign 9 points to the RNFL/GCL - retina/vitreous span
## assign 5 points from retina/vitreous to 30.41 into the vitreous.
## so, that's 100 points, slightly over-representing the outer retina span because that's where our hypothesis is focused.

depthstrip=c(y.RESOLUTION,markings.profile[,1])
Zz.distances.from.leftedge.pixels=seq(1,nrow(Retina.Points.cropped),1)


### now, the much more challenging step of making a spatially normalized image

## "RAS" for revised,adjusted,spatially normalized 
## (which will have RPE at sample spot at the 6th pixel above "bottom", and the retina/vitreous border )
RAS.FLATTENED.LIGHT.RETINA.cropped=LIGHT.cropped[,1:200,]
RAS.FLATTENED.DARK.RETINA.cropped=DARK.cropped[,1:200,]
RAS.FLATTENED.LIGHT.RETINA.cropped[,,]<-NA
RAS.FLATTENED.DARK.RETINA.cropped[,,]<-NA


WHICH.INDEX=function(x) which(abs(depthstrip-x)==min(na.rm=T,abs(depthstrip-x)))[1]
position.of.RPE.light=as.factor(seq(1,nrow(R.RPE.POSITION.LIGHT),1))
position.of.RPE.dark=as.factor(seq(1,nrow(R.RPE.POSITION.DARK),1))

##starts=as.vector(tapply(R.RPE.POSITION.LIGHT[,z],position.of.RPE.light,WHICH.INDEX))


### we'll want some "flattened" cropped images, while we're at it...
### (in reality, we just shift RPE peak to a standard position....)
LIGHT.flat.cropped=LIGHT.cropped[,1:180,]
LIGHT.flat.cropped[,,]<-NA
DARK.flat.cropped=DARK.cropped[,1:180,]
DARK.flat.cropped[,,]<-NA
END=dim(DARK.flat.cropped)[1]

for(z in 1:dim(LIGHT.flat.cropped)[3])
 {RPE.START.LIGHT=as.vector(tapply(R.RPE.POSITION.LIGHT[,z],position.of.RPE.light,WHICH.INDEX))[1:END];
  for(x in 1:length(RPE.START.LIGHT))
   {if(!(is.na(RPE.START.LIGHT[x])))
     {A=seq(RPE.START.LIGHT[x]-19,RPE.START.LIGHT[x]+160,1);
      A[which(A<1)]=1;
      A[which(A>dim(LIGHT.cropped)[2])]=dim(LIGHT.cropped)[2];
      LIGHT.flat.cropped[x,,z]=LIGHT.cropped[x,A,z]}}}
  ## ignore the warnings()
for(z in 1:dim(DARK.flat.cropped)[3])
 {RPE.START.DARK=as.vector(tapply(R.RPE.POSITION.DARK[,z],position.of.RPE.light,WHICH.INDEX))[1:END];
  for(x in 1:length(RPE.START.DARK))
   {if(!(is.na(RPE.START.DARK[x])))
     {A=seq(RPE.START.DARK[x]-19,RPE.START.DARK[x]+160,1);
      A[which(A<1)]=1;
      A[which(A>dim(LIGHT.cropped)[2])]=dim(LIGHT.cropped)[2];
      DARK.flat.cropped[x,,z]=DARK.cropped[x,A,z]}}}
  ## ignore the warnings()

LIGHT.flat.cropped[which(is.na(LIGHT.flat.cropped))]=0
DARK.flat.cropped[which(is.na(DARK.flat.cropped))]=0

###
###
###
### NEW FOR 2019_06_JUN13f:
### redo the above using revised estimates of where everything is at.
### in essence, until now, we've been working to get good-enough guesses of where everything is at. 

















#f.write.analyze(LIGHT.flat.cropped[,,],paste(FILENAME,"_AS--LIGHT",sep=""),size="float")
f.write.analyze(DARK.flat.cropped[,,],paste(FILENAME,"_AS--ALL",sep=""),size="float")







Zz.distances.from.opticN.microns=Retina.Points.cropped[,4]



### NEW AS OF DEC-2019
#image(DARK.flat.cropped[,,1])
image(Zz.distances.from.opticN.microns,(seq(-19,160,1)*y.RESOLUTION),DARK.flat.cropped[,,1],col=gray.colors(33),main="red = basement membrane",xlab="Distance from Optic Nerve (microns)",ylab="Distance from Basement Membrane (microns)")
abline(h=0,col=rgb(1,0,0,0.5),lwd=1,lty=2)
matlines(Zz.distances.from.opticN.microns,(R.INL.IPL.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.OLM.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.RNFL.GCL.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.VITREOUS.RETINA.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.ONL.OPL.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)


lock.to.BM.mean=t(DARK.flat.cropped[1:3,,1])
lock.to.BM.mean[,1]=seq(-19,160,1)*y.RESOLUTION
lock.to.BM.mean[,3]=seq(1,nrow(lock.to.BM.mean),1)
colnames(lock.to.BM.mean)<-c("Dist.micron.from.BM","signal.mean","column.of.flat.cropped")
for(x in 1:nrow(lock.to.BM.mean)) lock.to.BM.mean[x,2]=mean(na.rm=T, DARK.flat.cropped[,x,1])


dev.new()
plot(lock.to.BM.mean[,1],lock.to.BM.mean[,2],xlab="Distance from BM (microns)",ylab="Reflectivity",sub="red=BM")
abline(v=0,col=rgb(1,0,0,0.5),lwd=1,lty=2)
#savePlot(filename = paste("_lock_to_BM_plot_",FILENAME,"__",floor(runif(1000, 1,101))[1],".png",sep=""),
#         type = "png")

## we're going to do a moving window analysis to find features between the OLM and the RPE.
## the location of the first peak (dervative <= 0) as we move internal from the basement membrane will be used as the peak for the RPE signal.
## (+/- 1 pixel)

NEW.LANDMARKS.by.col.position.flat.cropped=DARK.flat.cropped[,1:10,1]
NEW.LANDMARKS.by.col.position.flat.cropped[,]<-NA
colnames(NEW.LANDMARKS.by.col.position.flat.cropped)<-c("RPE.peak","a","a","a","sub.IS.OS.trough.est.1","sub.IS.OS.trough.est.2","IS.OS.est.1","IS.OS.est.2","TROUGH.external.to.OLM","OLM")

lock.to.BM.mean=cbind(lock.to.BM.mean,lock.to.BM.mean[,2]-lock.to.BM.mean[c(2:nrow(lock.to.BM.mean),nrow(lock.to.BM.mean)),2])
calc=lock.to.BM.mean[20:nrow(lock.to.BM.mean),]
NEW.LANDMARKS.by.col.position.flat.cropped[,1]<-calc[which(calc[,4]>0)[1],3]

NEW.LANDMARKS.by.col.position.flat.cropped[,10]<-20+round((R.OLM.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1])/y.RESOLUTION)


## so, take the most vitread point in 
thickest.OLM.part.microns=mean(na.rm=T,(R.OLM.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]))
thickest.OLM.part.pix=ceiling(thickest.OLM.part.microns/y.RESOLUTION)
firstSLICE.OLM.to.RPE=DARK.flat.cropped[,20:(20+thickest.OLM.part.pix),1]

image(Zz.distances.from.opticN.microns,(seq(0,thickest.OLM.part.pix,1)*y.RESOLUTION),firstSLICE.OLM.to.RPE[,],col=gray.colors(33),main="red = basement membrane",xlab="Distance from Optic Nerve (microns)",ylab="Distance from Basement Membrane (microns)")

## we'll actually want things locked to the OLM for a moment....

lock.to.OLM=DARK.flat.cropped[,,1]
lock.to.OLM[,]<-NA
 OLM.START.DARK=as.vector(tapply(R.OLM.POSITION.DARK[,1],position.of.RPE.light,WHICH.INDEX))[1:END];
  ## ignore the warnings()
 for(x in 1:length(OLM.START.DARK))
  {if(!(is.na(OLM.START.DARK[x])))
    {A=seq((OLM.START.DARK[x]-19)-thickest.OLM.part.pix,(OLM.START.DARK[x]+160)-thickest.OLM.part.pix,1);
     A[which(A<1)]=1;
     A[which(A>dim(LIGHT.cropped)[2])]=dim(LIGHT.cropped)[2];
     lock.to.OLM[x,]=DARK.cropped[x,A,1]}}
image(Zz.distances.from.opticN.microns,(seq(-19,160,1)*y.RESOLUTION),lock.to.OLM[,],col=gray.colors(33),main="red = basement membrane",xlab="Distance from Optic Nerve (microns)",ylab="Distance from Basement Membrane (microns)")

## and OLM in that matrix is locked to "thickest.OLM.part.pix"+20 (so, the 55th column, if equal to thickest.OLM.part.pix=35)

lock.to.OLM.mean=t(lock.to.OLM[1:2,])
lock.to.OLM.mean[,1]=seq((-19-thickest.OLM.part.pix),160-thickest.OLM.part.pix,1)*y.RESOLUTION
colnames(lock.to.OLM.mean)<-c("Dist.micron.from.ELM","signal.mean")
for(x in 1:nrow(lock.to.OLM.mean)) lock.to.OLM.mean[x,2]=mean(na.rm=T, lock.to.OLM[,x])

image(Zz.distances.from.opticN.microns,(seq((-19)-thickest.OLM.part.pix,160-thickest.OLM.part.pix,1)*y.RESOLUTION),lock.to.OLM[,],col=gray.colors(33),main="red = ELM",xlab="Distance from Optic Nerve (microns)",ylab="Distance from ELM (microns)")
abline(h=0,col=rgb(1,0,0,0.5),lwd=1,lty=2)
savePlot(filename = paste("_lock_to_ELM_",FILENAME,"__",floor(runif(1000, 1,101))[1],".png",sep=""),
         type = "png")
dev.new()
plot(lock.to.OLM.mean[,1],lock.to.OLM.mean[,2],xlab="Distance from ELM (microns)",ylab="Reflectivity",sub="green=basement membrane, red=ELM")
abline(v=0,col=rgb(1,0,0,0.5),lwd=1,lty=2)
abline(v=(-1*thickest.OLM.part.microns),col=rgb(0,1,0,0.5),lwd=1,lty=2)
#savePlot(filename = paste("_lock_to_ELM_plot_",FILENAME,"__",floor(runif(1000, 1,101))[1],".png",sep=""),
#         type = "png")

lock.to.OLM.mean=cbind(lock.to.OLM.mean,lock.to.OLM.mean[,2]-lock.to.OLM.mean[c(2:nrow(lock.to.OLM.mean),nrow(lock.to.OLM.mean)),2])
## starting one pixel exterior to OLM, and moving exterior, find the first place with a positive difference (col 3) and move one interior to that. 
## that's the local trough. 
calc=lock.to.OLM.mean[which(lock.to.OLM.mean[,1]<0.01),]
calc[(nrow(calc)-1):nrow(calc),3]=NA
## how many pixels exterior?
trough=nrow(calc)-(which(calc[,3]>0)[length(which(calc[,3]>0))]+1)

NEW.LANDMARKS.by.col.position.flat.cropped[,9]<-(NEW.LANDMARKS.by.col.position.flat.cropped[,10]-trough)

## production rules: 
# in the OLM-locked profile (1) find the closest instance of RPE peak, and toss everything exterior to it.
toss.ext=min(na.rm=T, (NEW.LANDMARKS.by.col.position.flat.cropped[,10]-NEW.LANDMARKS.by.col.position.flat.cropped[,1]))
calc=calc[(nrow(calc)-toss.ext):nrow(calc),]
# and NA-out through to the trough adjacent to the OLM (+one pixel exterior)
calc[(which(calc[,3]>0)[length(which(calc[,3]>0))]):nrow(calc),3]=NA
# and due to how things are calculated, there's no guarantee we're at the expected RPE peak signal...
# so look for a peak within the first ~four values
peak=max(na.rm=T,calc[1:4,2])
calc[1:which(calc[1:4,2]==peak)[length(which(calc[1:4,2]==peak))],3]=NA
## now, NA-out a couple adjacent pixles
whatsleft=which(!(is.na(calc[,3])))
calc[whatsleft[1],3]=NA
calc[whatsleft[length(whatsleft)],3]=NA
calc[whatsleft[2],3]=NA
calc[whatsleft[(length(whatsleft)-1)],3]=NA

# in parallel, look for surface features within the BM-aligned profile, and shave off as much as possible...
calc2=lock.to.BM.mean[NEW.LANDMARKS.by.col.position.flat.cropped[1,1]:nrow(lock.to.BM.mean),]
calc2=calc2[1:toss.ext,]
# and NA-out through to the trough adjacent to the OLM (+one pixel exterior)
calc2[(nrow(calc2)-trough):nrow(calc2),4]=NA
## and the way things are calculated, we might not be quite at the true trough...
calc2=calc2[nrow(calc2):1,]
peak.inverse=min(na.rm=T,calc2[1:6,2])
calc2[1:which(calc2[1:6,2]==peak.inverse)[length(which(calc2[1:6,2]==peak.inverse))],4]=NA
calc2=calc2[nrow(calc2):1,]
## now, NA-out a couple adjacent pixles
whatsleft=which(!(is.na(calc2[,4])))
calc2[whatsleft[1],4]=NA
calc2[whatsleft[2],4]=NA
calc2[whatsleft[length(whatsleft)],4]=NA
calc2[whatsleft[(length(whatsleft)-1)],4]=NA

## production rules:
##  calc and calc2 currently hold the slopes of the mean contour (walking from inside to outsid the eye). 
## if there's at least one peak (walking in the same direction, we go from positive to negative)
## then (1) the first time this happens is where we should look for the IS/OS peak.
##      (2) the next place there's a trough is where we should look for the sub-IS/OS trough, 
##      (3) only between the sub-IS/OS trough and the RPE peak should we look for a third dip in light-exposed eyes
## but if there's not at least one peak (walking in the same direction, we go from positive to negative)
## then (1) the place with the gentlist slope is the best estimate for the IS/OS peak. 
##      (2) no further guesses about trough position, though as a crude estimate, if there is space between the RPE and the 
##          photoreceptors, then we shouldn't find it in the interior half of the span between IS/OS and RPE peak; it should be closer to the RPE.

calc=cbind(calc,calc[,3])
calc2=cbind(calc2,calc2[,4])

subISOS=NA
if(length(which(calc[,3]<0))>0)
  {ISOS=which(calc[,3]<0)[length(which(calc[,3]<0))]+1;
   calc[ISOS:nrow(calc),4]=NA;
   subISOS=which(calc[,4]>0)[length(which(calc[,4]>0))]+1;
   NEW.LANDMARKS.by.col.position.flat.cropped[,8]<-(NEW.LANDMARKS.by.col.position.flat.cropped[,10]+round(calc[ISOS,1]/y.RESOLUTION));
   NEW.LANDMARKS.by.col.position.flat.cropped[,6]<-(NEW.LANDMARKS.by.col.position.flat.cropped[,10]+round(calc[subISOS,1]/y.RESOLUTION))}
if(length(which(calc[,3]<0))==0)
  {ISOS=which(calc[,3]==min(na.rm=T,calc[,3]))[length(which(calc[,3]==min(na.rm=T,calc[,3])))]+1;
   NEW.LANDMARKS.by.col.position.flat.cropped[,8]<-(NEW.LANDMARKS.by.col.position.flat.cropped[,10]+round(calc[ISOS,1]/y.RESOLUTION))}

if(length(which(calc2[,4]<0))>0)
  {ISOS=which(calc2[,4]<0)[length(which(calc2[,4]<0))]+1;
   calc2[ISOS:nrow(calc2),5]=NA;
   subISOS=which(calc2[,5]>0)[length(which(calc2[,5]>0))]+1;
   NEW.LANDMARKS.by.col.position.flat.cropped[,7]<-calc2[ISOS,3];
   NEW.LANDMARKS.by.col.position.flat.cropped[,5]<-calc2[subISOS,3]}
if(length(which(calc2[,4]<0))==0)
  {ISOS=which(calc2[,4]==min(na.rm=T,calc2[,4]))[length(which(calc2[,4]==min(na.rm=T,calc2[,4])))]+1;
   NEW.LANDMARKS.by.col.position.flat.cropped[,7]<-calc2[ISOS,3]}

## OK! Now to loop through and look for local minimums and maximums:



## the moving window will be a bit more than (rounded up just beyond) 30 um wide 
# ...just use one from before:
#window.width.in.pixels=ceiling(30/x.RESOLUTION)
#MEAN.x=function(x) mean(na.rm=TRUE,x)
window.factor=matrix(,dim(DARK.flat.cropped)[2],window.width.in.pixels);
for(y in 1:nrow(window.factor))
 {window.factor[y,]=y}
window.factor=t(window.factor)
window.factor=as.factor(window.factor)

TROUGH.EXTERIOR.TO.OLM=matrix(,dim(DARK.flat.cropped)[1],5)
TROUGH.EXTERIOR.TO.OLM[,1]=Zz.distances.from.opticN.microns
colnames(TROUGH.EXTERIOR.TO.OLM)<-c("Distance.from.Optic.N","Column.in.flat.cropped","Trough.distance.from.BM.microns","Trough.width.microns","Trough.signal.not.normed")
ISOS=TROUGH.EXTERIOR.TO.OLM
colnames(ISOS)<-c("Distance.from.Optic.N","Column.in.flat.cropped","Peak.distance.from.BM.microns","Peak.width.microns","Peak.signal.not.normed")
GAP.JUST.INTERIOR.TO.RPE=cbind(TROUGH.EXTERIOR.TO.OLM,ISOS[,2:5])
colnames(GAP.JUST.INTERIOR.TO.RPE)<-c("Distance.from.Optic.N","Column.in.flat.cropped","Peak.distance.from.BM.microns","Trough.width.microns","Trough.signal.not.normed","RPE.signal.local.peak.not.normed","ROS.tip.signal.local.peak.not.normed","interior.dip.boder.microns","exterior.dip.border.microns")


DFC=DARK.flat.cropped[,,1]
DFC[which(DFC==0)]=NA




EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION=GAP.JUST.INTERIOR.TO.RPE[,1:3]
EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[,2]<-NA
EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[,3]<-NA
colnames(EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION)<-c("Distance.from.Optic.N","BM.to.RPE.peak.microns","BM.to.peak.for.OS.tip.microns")


### 2020_08AUG_12 note:


for(x in start.move:break.left)
 {window=DFC[(x-floor(window.width.in.pixels/2)+1):(x+floor(window.width.in.pixels/2)),];
  profile=cbind((seq(-19,160,1)*y.RESOLUTION),seq(1,180,1),as.vector(tapply(window,window.factor,MEAN.x)));
  profile=cbind(profile,profile[,3]-profile[c(2:nrow(profile),nrow(profile)),3]);
  plot(profile[,2],profile[,3])
  abline(v=(NEW.LANDMARKS.by.col.position.flat.cropped[x,]),lwd=2)
  ## focus down:
  ## we know that none of the landmarks we're logging will be beyond the OLM nor the RPE PEAK
  profile2=profile[(NEW.LANDMARKS.by.col.position.flat.cropped[x,1]):(NEW.LANDMARKS.by.col.position.flat.cropped[x,10]),];
  matlines(profile2[,2],profile2[,3],col="red");
  profile2.spline=smooth.spline(profile2[,1],profile2[,3],df=nrow(profile2));
  profile2.spline.detail=seq(min(profile2[,1]),max(profile2[,1]),y.RESOLUTION/10);
  profile2.spline.detail=cbind(profile2.spline.detail,as.numeric(predict(profile2.spline,profile2.spline.detail)$y));
  matlines((profile2.spline.detail[,1]/y.RESOLUTION)+20,profile2.spline.detail[,2],col="red",lwd=2);
  ##
  ## first, get trough exterior to OLM
  if(!(is.na(NEW.LANDMARKS.by.col.position.flat.cropped[x,9])))
    {search=profile[(NEW.LANDMARKS.by.col.position.flat.cropped[x,9]-1):(NEW.LANDMARKS.by.col.position.flat.cropped[x,9]+1),];
     TROUGH.EXTERIOR.TO.OLM[x,c(2,3,5)]=search[which(search[,3]==min(search[,3]))[1],c(2,1,3)];
     ## now, get ISOS (if present)(here, count as maximum in the proscribed area, even if not a true peak
     if( (!(is.na(NEW.LANDMARKS.by.col.position.flat.cropped[x,7]))) | (!(is.na(NEW.LANDMARKS.by.col.position.flat.cropped[x,8]))))
       {minny=min(na.rm=T,NEW.LANDMARKS.by.col.position.flat.cropped[x,7:8]);
        maxxy=max(na.rm=T,NEW.LANDMARKS.by.col.position.flat.cropped[x,7:8]);
        search=profile[(minny-1):(maxxy+1),];
        ISOS[x,c(2,3,5)]=search[which(search[,3]==max(search[,3]))[1],c(2,1,3)];
        ## as long as I have a value for ISOS peak and position, I can now use half-height method to extract width of 
        ## trough exterior to OLM
        search.exterior=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(ISOS[x,3],2)):which(round(profile2.spline.detail[,1],2)==round(TROUGH.EXTERIOR.TO.OLM[x,3],2)),];
        search.interior=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(TROUGH.EXTERIOR.TO.OLM[x,3],2)):nrow(profile2.spline.detail),];
        #
        mid.ext=(max(na.rm=T,search.exterior[,2])+min(na.rm=T,search.exterior[,2]))/2;
        search.exterior[,2]=abs(search.exterior[,2]-mid.ext)
        ext.border=search.exterior[which(search.exterior[,2]==min(na.rm=T,search.exterior[,2]))[1],1];
        #
        mid.int=(max(na.rm=T,search.interior[,2])+min(na.rm=T,search.interior[,2]))/2;
        search.interior[,2]=abs(search.interior[,2]-mid.int)
        int.border=search.interior[which(search.interior[,2]==min(na.rm=T,search.interior[,2]))[1],1];
        ## 
        TROUGH.EXTERIOR.TO.OLM[x,4]=int.border-ext.border;
        ## to proceed, there needs to be an actual peak in the vicinity of the IS/OS 
        ## (expand search view, and look for at least one negative slope (walking from inside to outside the eye)
        ## (and at least one positive slope interior to that).
        search=profile[(minny-6):(maxxy+4),];
        find=which(search[,4]<0)
        if(length(find)>0)
          {seek=search[find[1]:nrow(search),];
           if(find[1]==nrow(search)) {seek=rbind(search[nrow(search),],search[nrow(search),]); seek[,1]<-NA};
           if(length(which(seek[,4]>0))>0)
             {## find the local minimum
              if( (!(is.na(NEW.LANDMARKS.by.col.position.flat.cropped[x,5]))) | (!(is.na(NEW.LANDMARKS.by.col.position.flat.cropped[x,6]))))
                {minny=min(na.rm=T,NEW.LANDMARKS.by.col.position.flat.cropped[x,5:6]);
                 maxxy=max(na.rm=T,NEW.LANDMARKS.by.col.position.flat.cropped[x,5:6]);
                 search=profile[(minny-1):(maxxy+1),];
                 THE.TROUGH=search[which(search[,3]==min(search[,3]))[1],];
                  #*
                 if(round(ISOS[x,3],2)==THE.TROUGH[1]) THE.TROUGH=search[which(search[,3]==min(search[,3]))[1]-1,];
                  #*
                 ## as long as I have a these values for the trough position and value, I can now use half-height method to extract width of 
                 ## ISOS
                 search.exteriorISOS=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(THE.TROUGH[1],2)):which(round(profile2.spline.detail[,1],2)==round(ISOS[x,3],2)),];
                 mid.extISOS=(max(na.rm=T,search.exteriorISOS[,2])+min(na.rm=T,search.exteriorISOS[,2]))/2;
                 search.exteriorISOS[,2]=abs(search.exteriorISOS[,2]-mid.extISOS)
                 extISOS.border=search.exteriorISOS[which(search.exteriorISOS[,2]==min(na.rm=T,search.exteriorISOS[,2]))[1],1];
                 ISOS[x,4]=ext.border-extISOS.border;
                 ## NEXT UP, seach for a dip in the contour between THE.TROUGH and the position of the RPE peak
                 ## it must appear in the original profile
                 search=profile[(NEW.LANDMARKS.by.col.position.flat.cropped[x,1]):THE.TROUGH[2],];
                 ## and if by chance we're not starting at the peak signal within the first few points, fix that:
                 local.max=which(search[,3]==max(na.rm=T,search[1:4,3]))[1];
                 search=search[local.max:nrow(search),];
                 ## and when we look for a negative value for the slope (indicating a dip on an otherwise upward trend, walking from the inside to outside the eye)
                 ## toss out the slope within 1 of the peak and within 2 of the trough (asymmetric because of how it's calculated)
                 search[(nrow(search)-1):nrow(search),4]=NA;
                 search[1,4]=NA;
                 ## now, is there a negative slope?
                 dip=which(search[,4]<0);
                 ## if so, go to the spline
                 if(length(dip)>0)
                   {the.exterior.peak=search[1,];
                    ## if there's more than one dip, prioritize the one with > 1 contiguous negative slope.
                    ## if there's more than one dip, and they both have equal number of contiguous negative slopes,
                    ## the pick the dip with the greatest negative slope.
                    ## assign the interior-most point of a dip to be "THEDIP"
                    if(length(dip)==1) THEDIP=dip
                    if(length(dip)>1)
                      {discontinuous=cbind(search[dip,],search[dip[c(2:length(dip),length(dip))],2],search[dip[c(2:length(dip),length(dip))],2]);
                       discontinuous[,6]<-discontinuous[,5]-discontinuous[,2];
                       if(max(discontinuous[,6])==1) THEDIP=dip[length(dip)]
                       if(max(discontinuous[,6])>1)
                         {## for coding purposes, assume no more than three seperate dips
                          dip1=discontinuous[1:2,];
                          dip1[,]<-NA;
                          dip2=dip1;
                          dip3=dip1;
                          ## work in reverse
                          counter=1
                          for(y in nrow(discontinuous):1)
                            {if(discontinuous[y,6]>1) counter=counter+1;
                             if(counter==1) dip1=rbind(dip1,discontinuous[y,]);
                             if(counter==2) dip2=rbind(dip2,discontinuous[y,]);
                             if(counter==3) dip3=rbind(dip3,discontinuous[y,])};
                          dip1=cbind(dip1,matrix(1,nrow(dip1),1));
                          dip2=cbind(dip2,matrix(2,nrow(dip2),1));
                          dip3=cbind(dip3,matrix(3,nrow(dip3),1));
                          dipANALYZE=rbind(dip1,dip2,dip3);
                          freqTable=table(dipANALYZE[,7]);
                          toExclude=as.numeric(which(freqTable==min(freqTable)));
                          ## if the length of toExclude is < 3, then it means one of the three possible dips has fewer contiguous negative slopes
                          ## if the length of toExclude is == 3, then there are (at least) three dips of equal congtiguity... 
                          if(length(toExclude)<3) for(z in 1:length(toExclude)) dipANALYZE=dipANALYZE[which(!(dipANALYZE[,7]==toExclude[z])),];
                          ## and from whatever is left, find the most negative slope
                          dipANALYZE=dipANALYZE[which(!(is.na(dipANALYZE[,1]))),];
                          WHICHDIP=dipANALYZE[which(dipANALYZE[,4]==min(dipANALYZE[,4]))[1],7];
                          for(z in 1:nrow(dipANALYZE)) if(!(dipANALYZE[z,7]==WHICHDIP)) dipANALYZE[z,]<-NA;
                          THEDIP=which(search[,1]==max(na.rm=T,dipANALYZE[,1]))}};
                    the.interior.peak=search[THEDIP+1,];
                    dip.search=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(the.exterior.peak[1],2)):which(round(profile2.spline.detail[,1],2)==round(the.interior.peak[1],2)),];
                    ## find the min, and the interior and exterior half-heights
                    the.dip=dip.search[which(dip.search[,2]==min(dip.search[,2])[1]),]
                    GAP.JUST.INTERIOR.TO.RPE[x,c(3,5)]=the.dip[c(1,2)];
                    GAP.JUST.INTERIOR.TO.RPE[x,6]=the.exterior.peak[3];
                    GAP.JUST.INTERIOR.TO.RPE[x,7]=the.interior.peak[3];
                    find.column=profile;
                    find.column[,1]<-abs(find.column[,1]-the.dip[1]);
                    GAP.JUST.INTERIOR.TO.RPE[x,2]=find.column[which(find.column[,1]==min(na.rm=T,find.column[,1]))[1],2];
                    interior.to.dip=dip.search[which(dip.search[,1]==the.dip[1]):nrow(dip.search),];
                    midway=(min(interior.to.dip[,2])+max(interior.to.dip[,2]))/2;
                    interior.to.dip[,2]<-abs(interior.to.dip[,2]-midway);
                    inner.dip.border=interior.to.dip[which(interior.to.dip[,2]==min(interior.to.dip[,2]))[1],1];
                    exterior.to.dip=dip.search[1:which(dip.search[,1]==the.dip[1]),];
			  ### <new for 2020_08AUG_12>
                    if(length(exterior.to.dip)==2) exterior.to.dip=rbind(exterior.to.dip,exterior.to.dip);
			  ### </new for 2020_08AUG_12>
                    midway=(min(exterior.to.dip[,2])+max(exterior.to.dip[,2]))/2;
                    exterior.to.dip[,2]<-abs(exterior.to.dip[,2]-midway);
                    outer.dip.border=exterior.to.dip[which(exterior.to.dip[,2]==min(exterior.to.dip[,2]))[1],1];
                    GAP.JUST.INTERIOR.TO.RPE[x,8]=inner.dip.border;
                    GAP.JUST.INTERIOR.TO.RPE[x,9]=outer.dip.border;
                    GAP.JUST.INTERIOR.TO.RPE[x,4]=inner.dip.border-outer.dip.border;
                    EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[x,3]=search[1,1];
                    EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[x,2]=search[THEDIP+1,1]}}}}}}}

for(x in break.right:end.move)
 {window=DFC[(x-floor(window.width.in.pixels/2)+1):(x+floor(window.width.in.pixels/2)),];
  profile=cbind((seq(-19,160,1)*y.RESOLUTION),seq(1,180,1),as.vector(tapply(window,window.factor,MEAN.x)));
  profile=cbind(profile,profile[,3]-profile[c(2:nrow(profile),nrow(profile)),3]);
  plot(profile[,2],profile[,3])
  abline(v=(NEW.LANDMARKS.by.col.position.flat.cropped[x,]),lwd=2)
  ## focus down:
  ## we know that none of the landmarks we're logging will be beyond the OLM nor the RPE PEAK
  profile2=profile[(NEW.LANDMARKS.by.col.position.flat.cropped[x,1]):(NEW.LANDMARKS.by.col.position.flat.cropped[x,10]),];
  matlines(profile2[,2],profile2[,3],col="red");
  profile2.spline=smooth.spline(profile2[,1],profile2[,3],df=nrow(profile2));
  profile2.spline.detail=seq(min(profile2[,1]),max(profile2[,1]),y.RESOLUTION/10);
  profile2.spline.detail=cbind(profile2.spline.detail,as.numeric(predict(profile2.spline,profile2.spline.detail)$y));
  matlines((profile2.spline.detail[,1]/y.RESOLUTION)+20,profile2.spline.detail[,2],col="red",lwd=2);
  ##
  ## first, get trough exterior to OLM
  if(!(is.na(NEW.LANDMARKS.by.col.position.flat.cropped[x,9])))
    {search=profile[(NEW.LANDMARKS.by.col.position.flat.cropped[x,9]-1):(NEW.LANDMARKS.by.col.position.flat.cropped[x,9]+1),];
     TROUGH.EXTERIOR.TO.OLM[x,c(2,3,5)]=search[which(search[,3]==min(search[,3]))[1],c(2,1,3)];
     ## now, get ISOS (if present)(here, count as maximum in the proscribed area, even if not a true peak
     if( (!(is.na(NEW.LANDMARKS.by.col.position.flat.cropped[x,7]))) | (!(is.na(NEW.LANDMARKS.by.col.position.flat.cropped[x,8]))))
       {minny=min(na.rm=T,NEW.LANDMARKS.by.col.position.flat.cropped[x,7:8]);
        maxxy=max(na.rm=T,NEW.LANDMARKS.by.col.position.flat.cropped[x,7:8]);
        search=profile[(minny-1):(maxxy+1),];
        ISOS[x,c(2,3,5)]=search[which(search[,3]==max(search[,3]))[1],c(2,1,3)];
        ## as long as I have a value for ISOS peak and position, I can now use half-height method to extract width of 
        ## trough exterior to OLM
        search.exterior=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(ISOS[x,3],2)):which(round(profile2.spline.detail[,1],2)==round(TROUGH.EXTERIOR.TO.OLM[x,3],2)),];
        search.interior=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(TROUGH.EXTERIOR.TO.OLM[x,3],2)):nrow(profile2.spline.detail),];
        #
        mid.ext=(max(na.rm=T,search.exterior[,2])+min(na.rm=T,search.exterior[,2]))/2;
        search.exterior[,2]=abs(search.exterior[,2]-mid.ext)
        ext.border=search.exterior[which(search.exterior[,2]==min(na.rm=T,search.exterior[,2]))[1],1];
        #
        mid.int=(max(na.rm=T,search.interior[,2])+min(na.rm=T,search.interior[,2]))/2;
        search.interior[,2]=abs(search.interior[,2]-mid.int)
        int.border=search.interior[which(search.interior[,2]==min(na.rm=T,search.interior[,2]))[1],1];
        ## 
        TROUGH.EXTERIOR.TO.OLM[x,4]=int.border-ext.border;
        ## to proceed, there needs to be an actual peak in the vicinity of the IS/OS 
        ## (expand search view, and look for at least one negative slope (walking from inside to outside the eye)
        ## (and at least one positive slope interior to that).
        search=profile[(minny-6):(maxxy+4),];
        find=which(search[,4]<0)
        if(length(find)>0)
          {seek=search[find[1]:nrow(search),];
           if(find[1]==nrow(search)) {seek=rbind(search[nrow(search),],search[nrow(search),]); seek[,1]<-NA};
           if(length(which(seek[,4]>0))>0)
             {## find the local minimum
              if( (!(is.na(NEW.LANDMARKS.by.col.position.flat.cropped[x,5]))) | (!(is.na(NEW.LANDMARKS.by.col.position.flat.cropped[x,6]))))
                {minny=min(na.rm=T,NEW.LANDMARKS.by.col.position.flat.cropped[x,5:6]);
                 maxxy=max(na.rm=T,NEW.LANDMARKS.by.col.position.flat.cropped[x,5:6]);
                 search=profile[(minny-1):(maxxy+1),];
                 THE.TROUGH=search[which(search[,3]==min(search[,3]))[1],];
                  #*
                 if(round(ISOS[x,3],2)==THE.TROUGH[1]) THE.TROUGH=search[which(search[,3]==min(search[,3]))[1]-1,];
                  #*
                 ## as long as I have a these values for the trough position and value, I can now use half-height method to extract width of 
                 ## ISOS
                 search.exteriorISOS=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(THE.TROUGH[1],2)):which(round(profile2.spline.detail[,1],2)==round(ISOS[x,3],2)),];
                 mid.extISOS=(max(na.rm=T,search.exteriorISOS[,2])+min(na.rm=T,search.exteriorISOS[,2]))/2;
                 search.exteriorISOS[,2]=abs(search.exteriorISOS[,2]-mid.extISOS)
                 extISOS.border=search.exteriorISOS[which(search.exteriorISOS[,2]==min(na.rm=T,search.exteriorISOS[,2]))[1],1];
                 ISOS[x,4]=ext.border-extISOS.border;
                 ## NEXT UP, seach for a dip in the contour between THE.TROUGH and the position of the RPE peak
                 ## it must appear in the original profile
                 search=profile[(NEW.LANDMARKS.by.col.position.flat.cropped[x,1]):THE.TROUGH[2],];
                 ## and if by chance we're not starting at the peak signal within the first few points, fix that:
                 local.max=which(search[,3]==max(na.rm=T,search[1:4,3]))[1];
                 search=search[local.max:nrow(search),];
                 ## and when we look for a negative value for the slope (indicating a dip on an otherwise upward trend, walking from the inside to outside the eye)
                 ## toss out the slope within 1 of the peak and within 2 of the trough (asymmetric because of how it's calculated)
                 search[(nrow(search)-1):nrow(search),4]=NA;
                 search[1,4]=NA;
                 ## now, is there a negative slope?
                 dip=which(search[,4]<0);
                 ## if so, go to the spline
                 if(length(dip)>0)
                   {the.exterior.peak=search[1,];
                    ## if there's more than one dip, prioritize the one with > 1 contiguous negative slope.
                    ## if there's more than one dip, and they both have equal number of contiguous negative slopes,
                    ## the pick the dip with the greatest negative slope.
                    ## assign the interior-most point of a dip to be "THEDIP"
                    if(length(dip)==1) THEDIP=dip
                    if(length(dip)>1)
                      {discontinuous=cbind(search[dip,],search[dip[c(2:length(dip),length(dip))],2],search[dip[c(2:length(dip),length(dip))],2]);
                       discontinuous[,6]<-discontinuous[,5]-discontinuous[,2];
                       if(max(discontinuous[,6])==1) THEDIP=dip[length(dip)]
                       if(max(discontinuous[,6])>1)
                         {## for coding purposes, assume no more than three seperate dips
                          dip1=discontinuous[1:2,];
                          dip1[,]<-NA;
                          dip2=dip1;
                          dip3=dip1;
                          ## work in reverse
                          counter=1
                          for(y in nrow(discontinuous):1)
                            {if(discontinuous[y,6]>1) counter=counter+1;
                             if(counter==1) dip1=rbind(dip1,discontinuous[y,]);
                             if(counter==2) dip2=rbind(dip2,discontinuous[y,]);
                             if(counter==3) dip3=rbind(dip3,discontinuous[y,])};
                          dip1=cbind(dip1,matrix(1,nrow(dip1),1));
                          dip2=cbind(dip2,matrix(2,nrow(dip2),1));
                          dip3=cbind(dip3,matrix(3,nrow(dip3),1));
                          dipANALYZE=rbind(dip1,dip2,dip3);
                          freqTable=table(dipANALYZE[,7]);
                          toExclude=as.numeric(which(freqTable==min(freqTable)));
                          ## if the length of toExclude is < 3, then it means one of the three possible dips has fewer contiguous negative slopes
                          ## if the length of toExclude is == 3, then there are (at least) three dips of equal congtiguity... 
                          if(length(toExclude)<3) for(z in 1:length(toExclude)) dipANALYZE=dipANALYZE[which(!(dipANALYZE[,7]==toExclude[z])),];
                          ## and from whatever is left, find the most negative slope
                          dipANALYZE=dipANALYZE[which(!(is.na(dipANALYZE[,1]))),];
                          WHICHDIP=dipANALYZE[which(dipANALYZE[,4]==min(dipANALYZE[,4]))[1],7];
                          for(z in 1:nrow(dipANALYZE)) if(!(dipANALYZE[z,7]==WHICHDIP)) dipANALYZE[z,]<-NA;
                          THEDIP=which(search[,1]==max(na.rm=T,dipANALYZE[,1]))}};
                    the.interior.peak=search[THEDIP+1,];
                    dip.search=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(the.exterior.peak[1],2)):which(round(profile2.spline.detail[,1],2)==round(the.interior.peak[1],2)),];
                    ## find the min, and the interior and exterior half-heights
                    the.dip=dip.search[which(dip.search[,2]==min(dip.search[,2])[1]),]
                    GAP.JUST.INTERIOR.TO.RPE[x,c(3,5)]=the.dip[c(1,2)];
                    GAP.JUST.INTERIOR.TO.RPE[x,6]=the.exterior.peak[3];
                    GAP.JUST.INTERIOR.TO.RPE[x,7]=the.interior.peak[3];
                    find.column=profile;
                    find.column[,1]<-abs(find.column[,1]-the.dip[1]);
                    GAP.JUST.INTERIOR.TO.RPE[x,2]=find.column[which(find.column[,1]==min(na.rm=T,find.column[,1]))[1],2];
                    interior.to.dip=dip.search[which(dip.search[,1]==the.dip[1]):nrow(dip.search),];
                    midway=(min(interior.to.dip[,2])+max(interior.to.dip[,2]))/2;
                    interior.to.dip[,2]<-abs(interior.to.dip[,2]-midway);
                    inner.dip.border=interior.to.dip[which(interior.to.dip[,2]==min(interior.to.dip[,2]))[1],1];
                    exterior.to.dip=dip.search[1:which(dip.search[,1]==the.dip[1]),];
			  ### <new for 2020_08AUG_12>
                    if(length(exterior.to.dip)==2) exterior.to.dip=rbind(exterior.to.dip,exterior.to.dip);
			  ### </new for 2020_08AUG_12>
                    midway=(min(exterior.to.dip[,2])+max(exterior.to.dip[,2]))/2;
                    exterior.to.dip[,2]<-abs(exterior.to.dip[,2]-midway);
                    outer.dip.border=exterior.to.dip[which(exterior.to.dip[,2]==min(exterior.to.dip[,2]))[1],1];
                    GAP.JUST.INTERIOR.TO.RPE[x,8]=inner.dip.border;
                    GAP.JUST.INTERIOR.TO.RPE[x,9]=outer.dip.border;
                    GAP.JUST.INTERIOR.TO.RPE[x,4]=inner.dip.border-outer.dip.border;
                    EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[x,3]=search[1,1];
                    EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[x,2]=search[THEDIP+1,1]}}}}}}}


#image(DARK.flat.cropped[,,1])
image(Zz.distances.from.opticN.microns,(seq(-19,160,1)*y.RESOLUTION),DARK.flat.cropped[,,1],col=gray.colors(33),main="red = basement membrane",xlab="Distance from Optic Nerve (microns)",ylab="Distance from Basement Membrane (microns)")
abline(h=0,col=rgb(1,0,0,0.5),lwd=1,lty=2)
matlines(Zz.distances.from.opticN.microns,(R.INL.IPL.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.OLM.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.RNFL.GCL.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.VITREOUS.RETINA.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.ONL.OPL.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
if(length(which(!(is.na(GAP.JUST.INTERIOR.TO.RPE[,8]))))>0) matlines(Zz.distances.from.opticN.microns,GAP.JUST.INTERIOR.TO.RPE[,8],col=rgb(0,0,1,0.5),lty=2)
if(length(which(!(is.na(GAP.JUST.INTERIOR.TO.RPE[,9]))))>0) matlines(Zz.distances.from.opticN.microns,GAP.JUST.INTERIOR.TO.RPE[,9],col=rgb(0,0,1,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,TROUGH.EXTERIOR.TO.OLM[,3],col=rgb(0,0,1),lty=1)
matlines(Zz.distances.from.opticN.microns,ISOS[,3],col=rgb(0,0,1),lty=2)

## we'll use some similar techniques as before to filter out bad values (probably won't replace them with anything, though). 
## assuming there are enough points...
if(length(which(!(is.na(GAP.JUST.INTERIOR.TO.RPE[,8]))))>10) 
 {z.threshold=2;
  split=round(nrow(GAP.JUST.INTERIOR.TO.RPE)/2);
  end=nrow(GAP.JUST.INTERIOR.TO.RPE);
  x=1;
  LEFT=GAP.JUST.INTERIOR.TO.RPE[1:split,c(1,3)];
  RIGHT=GAP.JUST.INTERIOR.TO.RPE[split:end,c(1,3)];
  FOUR=LEFT;
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,LEFT[,1])$y);
  #plot(LEFT[,1],RESULT.4,col="red");
  DEVIATION.4=LEFT[,2]-RESULT.4;
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  # 
  # 
  overall.z=cbind(seq(1,split,1),abs(z.DEVIATION.4));
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,3]<-ifelse(overall.z[,2]>z.threshold,1,0);
  overall.z.LEFT=overall.z;
  ##
  ###-----
  ## 
  FOUR=RIGHT;
  FOUR=FOUR[which(!(is.na(FOUR[,2]))),];
  MODEL.4=smooth.spline(FOUR[,1],FOUR[,2],df=5);
  RESULT.4=as.numeric(predict(MODEL.4,RIGHT[,1])$y);
  #plot(RIGHT[,1],RESULT.4,col="red");
  DEVIATION.4=RIGHT[,2]-RESULT.4;
  z.DEVIATION.4=DEVIATION.4;
  z.DEVIATION.4=(z.DEVIATION.4-mean(na.rm=T,DEVIATION.4))/sd(na.rm=T,DEVIATION.4);
  # 
  # 
  overall.z=cbind(seq(split,end,1),abs(z.DEVIATION.4));
  ## 
  ## if abs(overall.z) > z.threshold, replace the misbehaving data in RPE.POSITION with whatever a z score of 0 would be. 
  overall.z=cbind(overall.z,overall.z[,1]);
  overall.z[,3]<-ifelse(overall.z[,2]>z.threshold,1,0);
  ## now, I can join overall.z's, and use that to toss errors in finding the values for the gap.
  ZZZ=rbind(overall.z.LEFT,overall.z[2:nrow(overall.z),])
  GAP.JUST.INTERIOR.TO.RPE=cbind(GAP.JUST.INTERIOR.TO.RPE,ZZZ[,3])
  for(x in 1:nrow(GAP.JUST.INTERIOR.TO.RPE)) if(!(is.na(GAP.JUST.INTERIOR.TO.RPE[x,10]))) if(GAP.JUST.INTERIOR.TO.RPE[x,10]==1) GAP.JUST.INTERIOR.TO.RPE[x,2:9]=NA;
  GAP.JUST.INTERIOR.TO.RPE=GAP.JUST.INTERIOR.TO.RPE[,1:9]}


image(Zz.distances.from.opticN.microns,(seq(-19,160,1)*y.RESOLUTION),DARK.flat.cropped[,,1],col=gray.colors(33),main="red = basement membrane",xlab="Distance from Optic Nerve (microns)",ylab="Distance from Basement Membrane (microns)")
abline(h=0,col=rgb(1,0,0,0.5),lwd=1,lty=2)
matlines(Zz.distances.from.opticN.microns,(R.INL.IPL.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.OLM.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.RNFL.GCL.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.VITREOUS.RETINA.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,(R.ONL.OPL.POSITION.LIGHT[,1]-R.RPE.POSITION.LIGHT[,1]),col=rgb(0,1,0,0.5),lty=2)
if(length(which(!(is.na(GAP.JUST.INTERIOR.TO.RPE[,8]))))>0) matlines(Zz.distances.from.opticN.microns,GAP.JUST.INTERIOR.TO.RPE[,8],col=rgb(0,0,1,0.5),lty=2)
if(length(which(!(is.na(GAP.JUST.INTERIOR.TO.RPE[,9]))))>0) matlines(Zz.distances.from.opticN.microns,GAP.JUST.INTERIOR.TO.RPE[,9],col=rgb(0,0,1,0.5),lty=2)
matlines(Zz.distances.from.opticN.microns,TROUGH.EXTERIOR.TO.OLM[,3],col=rgb(0,0,1),lty=1)
matlines(Zz.distances.from.opticN.microns,ISOS[,3],col=rgb(0,0,1),lty=2)

savePlot(filename = paste("_lock_to_ELM_extra_layers",FILENAME,"__",floor(runif(1000, 1,101))[1],".png",sep=""),
         type = "png")


## now, more file outputs....

mean.distance.from.BM.to.ISOS=mean(na.rm=T,ISOS[,3])
mean.width.ISOS=mean(na.rm=T,ISOS[,4])
mean.peak.signal.ISOS=mean(na.rm=T,ISOS[,4])

## total span of retina pictured
span.pictured=length(which(!(is.na(R.OLM.POSITION.LIGHT[,2]))))*x.RESOLUTION
span.with.ISOS.to.RPE.gap=length(which(!(is.na(GAP.JUST.INTERIOR.TO.RPE[,8]))))*x.RESOLUTION
## all the below measures are means only of places where there's a gap.
mean.thickness.ISOS.to.RPE.gap=NA
mean.trough.signal.ISOS.to.RPE.gap=NA
mean.peak.signal.RPE.near.gap=NA
mean.peak.signal.OStip.near.gap=NA
mean.position.RPE.internal.border.microns.when.near.gap=NA
mean.position.OS.external.border.microns.when.near.gap=NA

if(length(which(!(is.na(GAP.JUST.INTERIOR.TO.RPE[,8]))))>0)
  {mean.thickness.ISOS.to.RPE.gap=mean(na.rm=T,GAP.JUST.INTERIOR.TO.RPE[,4]);
   mean.trough.signal.ISOS.to.RPE.gap=mean(na.rm=T,GAP.JUST.INTERIOR.TO.RPE[,5]);
   mean.peak.signal.RPE.near.gap=mean(na.rm=T,GAP.JUST.INTERIOR.TO.RPE[,6]);
   mean.peak.signal.OStip.near.gap=mean(na.rm=T,GAP.JUST.INTERIOR.TO.RPE[,7]);
   mean.position.RPE.internal.border.microns.when.near.gap=mean(na.rm=T,GAP.JUST.INTERIOR.TO.RPE[,9]);
   mean.position.OS.external.border.microns.when.near.gap=mean(na.rm=T,GAP.JUST.INTERIOR.TO.RPE[,8])}

## rmake output
SUBRETINAL.SPACE=as.data.frame(matrix(,(4+ncol(GAP.JUST.INTERIOR.TO.RPE)),(10+nrow(GAP.JUST.INTERIOR.TO.RPE))))
SUBRETINAL.SPACE[1,1]=FILENAME
SUBRETINAL.SPACE[2,1]="all_distances_relative_to_basement_membrane_in_um"
SUBRETINAL.SPACE[3,1]="all_intensities_in_same_units_as_input_images"
SUBRETINAL.SPACE[4,1]="covering_span_from_125_to_625um_from_OpticN"
SUBRETINAL.SPACE[5,1]="assumes_axial_microns_per_pixel_of"
SUBRETINAL.SPACE[8,1]="assumes_x_resolution_microns_per_pixel_of"
SUBRETINAL.SPACE[6,1]=y.RESOLUTION
SUBRETINAL.SPACE[9,1]=x.RESOLUTION

SUBRETINAL.SPACE[1,8]="BELOW"
SUBRETINAL.SPACE[2,8]="HERE"
SUBRETINAL.SPACE[3,8]="CONTINUOUS"
SUBRETINAL.SPACE[4,8]="DATA"
SUBRETINAL.SPACE[5,8]="_"
SUBRETINAL.SPACE[6,8]="_"
SUBRETINAL.SPACE[7,8]="_"
SUBRETINAL.SPACE[8,8]="_"
SUBRETINAL.SPACE[9,8]="_"
SUBRETINAL.SPACE[10,8]="_"
SUBRETINAL.SPACE[11,8]="_"
SUBRETINAL.SPACE[12,8]="_"
SUBRETINAL.SPACE[13,8]="_"

SUBRETINAL.SPACE[1,3]="Total_analyzed_span_microns"
SUBRETINAL.SPACE[2,3]=span.pictured
SUBRETINAL.SPACE[4,3]="span_with_RPE_OStip_gap_microns"
SUBRETINAL.SPACE[5,3]=span.with.ISOS.to.RPE.gap
SUBRETINAL.SPACE[7,3]="PERCENT_span_with_RPE_OStip_gap_microns"
SUBRETINAL.SPACE[8,3]=span.with.ISOS.to.RPE.gap/span.pictured

SUBRETINAL.SPACE[1,5]="-"
SUBRETINAL.SPACE[2,5]="ALL_LOCATIONS"
SUBRETINAL.SPACE[3,5]="ALL_LOCATIONS"
SUBRETINAL.SPACE[4,5]="ALL_LOCATIONS"
SUBRETINAL.SPACE[5,5]="ONLY_WHERE_RPE_OS_GAP_PRESENT"
SUBRETINAL.SPACE[6,5]="ONLY_WHERE_RPE_OS_GAP_PRESENT"
SUBRETINAL.SPACE[7,5]="ONLY_WHERE_RPE_OS_GAP_PRESENT"
SUBRETINAL.SPACE[8,5]="ONLY_WHERE_RPE_OS_GAP_PRESENT"
SUBRETINAL.SPACE[9,5]="ONLY_WHERE_RPE_OS_GAP_PRESENT"
SUBRETINAL.SPACE[10,5]="ONLY_WHERE_RPE_OS_GAP_PRESENT"
SUBRETINAL.SPACE[11,5]="ONLY_WHERE_RPE_OS_GAP_PRESENT"

SUBRETINAL.SPACE[1,10]="Distance_from_optic_nerve_microns"
SUBRETINAL.SPACE[1,11:ncol(SUBRETINAL.SPACE)]=GAP.JUST.INTERIOR.TO.RPE[,1]

SUBRETINAL.SPACE[2,10]="Distance_from_BM_to_ISOS_signal_peak"
SUBRETINAL.SPACE[2,11:ncol(SUBRETINAL.SPACE)]=ISOS[,3]

SUBRETINAL.SPACE[3,10]="Width_ISOS_band_microns"
SUBRETINAL.SPACE[3,11:ncol(SUBRETINAL.SPACE)]=ISOS[,4]

SUBRETINAL.SPACE[4,10]="Signal_at_ISOS_peak"
SUBRETINAL.SPACE[4,11:ncol(SUBRETINAL.SPACE)]=ISOS[,5]

SUBRETINAL.SPACE[5,10]="Distance_from_BM_to_interior_border_RPE_microns"
SUBRETINAL.SPACE[5,11:ncol(SUBRETINAL.SPACE)]=GAP.JUST.INTERIOR.TO.RPE[,9]

SUBRETINAL.SPACE[6,10]="Distance_from_BM_to_exterior_border_OS_tips_microns"
SUBRETINAL.SPACE[6,11:ncol(SUBRETINAL.SPACE)]=GAP.JUST.INTERIOR.TO.RPE[,8]

SUBRETINAL.SPACE[7,10]="width_gap_OStips_and_RPE"
SUBRETINAL.SPACE[7,11:ncol(SUBRETINAL.SPACE)]=GAP.JUST.INTERIOR.TO.RPE[,4]

SUBRETINAL.SPACE[8,10]="local_peak_RPE_signal"
SUBRETINAL.SPACE[8,11:ncol(SUBRETINAL.SPACE)]=GAP.JUST.INTERIOR.TO.RPE[,6]

SUBRETINAL.SPACE[9,10]="local_peak_OS_tip_signal"
SUBRETINAL.SPACE[9,11:ncol(SUBRETINAL.SPACE)]=GAP.JUST.INTERIOR.TO.RPE[,7]

SUBRETINAL.SPACE[10,10]="local_min_signal_in_gap_OStips_and_RPE"
SUBRETINAL.SPACE[10,11:ncol(SUBRETINAL.SPACE)]=GAP.JUST.INTERIOR.TO.RPE[,5]

SUBRETINAL.SPACE[11,10]="distance_from_BM_to_local_min_signal_in_gap"
SUBRETINAL.SPACE[11,11:ncol(SUBRETINAL.SPACE)]=GAP.JUST.INTERIOR.TO.RPE[,3]


EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION.special=EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION
EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION.special[,2]<-ifelse(is.na(GAP.JUST.INTERIOR.TO.RPE[,3]),NA,EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION.special[,2])
EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION.special[,3]<-ifelse(is.na(GAP.JUST.INTERIOR.TO.RPE[,3]),NA,EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION.special[,3])


SUBRETINAL.SPACE[12,10]="distance_from_BM_to_local_max_OStip"
SUBRETINAL.SPACE[12,11:ncol(SUBRETINAL.SPACE)]=EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION.special[,2]

SUBRETINAL.SPACE[13,10]="distance_from_BM_to_local_max_RPE"
SUBRETINAL.SPACE[13,11:ncol(SUBRETINAL.SPACE)]=EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION.special[,3]

rm(EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION.special)

SUBRETINAL.SPACE[2,6]="mean_distance_microns_from_BM_to_ISOS_band"
SUBRETINAL.SPACE[2,7]=mean.distance.from.BM.to.ISOS

SUBRETINAL.SPACE[3,6]="mean_width_ISOS_BAND"
SUBRETINAL.SPACE[3,7]=mean.width.ISOS

SUBRETINAL.SPACE[4,6]="mean_peak_signal_ISOS_BAND"
SUBRETINAL.SPACE[4,7]=mean.peak.signal.ISOS

SUBRETINAL.SPACE[5,6]="mean_distance_from_BM_to_interior_border_RPE_microns"
SUBRETINAL.SPACE[5,7]=mean.position.RPE.internal.border.microns.when.near.gap

SUBRETINAL.SPACE[6,6]="mean_distance_from_BM_to_exterior_border_OS_tips_microns"
SUBRETINAL.SPACE[6,7]=mean.position.OS.external.border.microns.when.near.gap

SUBRETINAL.SPACE[7,6]="width_gap_OStips_and_RPE"
SUBRETINAL.SPACE[7,7]=mean.thickness.ISOS.to.RPE.gap

SUBRETINAL.SPACE[8,6]="local_peak_RPE_signal"
SUBRETINAL.SPACE[8,7]=mean.peak.signal.RPE.near.gap

SUBRETINAL.SPACE[9,6]="local_peak_OS_tip_signal"
SUBRETINAL.SPACE[9,7]=mean.peak.signal.OStip.near.gap

SUBRETINAL.SPACE[10,6]="local_min_signal_in_gap_OStips_and_RPE"
SUBRETINAL.SPACE[10,7]=mean.trough.signal.ISOS.to.RPE.gap




### <END> NEW AS OF DEC-2019
#write(as.matrix(SUBRETINAL.SPACE),ncol=nrow(SUBRETINAL.SPACE),file=paste(FILENAME,"_SubretinalSPACEs.txt",sep=""))


### <START> NEW AS OF AUG-2020
### clean up output to match _first-image-only_thickness_details_c
###
SUBRETINAL.SPACE.OUT=SUBRETINAL.SPACE[,c(1:10,21:380,558:915)]
#write(as.matrix(SUBRETINAL.SPACE.OUT),ncol=nrow(SUBRETINAL.SPACE.OUT),file=paste(FILENAME,"_SubretinalSPACEs.txt",sep=""))

### <END> NEW AS OF AUG-2020


### <START> NEW AS OF SEP-2020

write(as.matrix(SUBRETINAL.SPACE.OUT),ncol=nrow(SUBRETINAL.SPACE.OUT),file=paste(FILENAME,"_tabdelimated_SubretinalSPACEs.txt",sep=""),sep="\t")


### <END> NEW AS OF SEP-2020






## assign 6*2 points from -30 to 0 [inclusive of lower value and upper, all others are inclusive of upper value only]
## assign 19*2 points to the RPE-OLM span
## assign 18*2 points to the OLM-ONL/OPL span
## assign 18*2 points to the ONL/OPL - INL/IPL span
## assign 25*2 points to the INL/IPL - RNFL/GCL span
## assign 9*2 points to the RNFL/GCL - retina/vitreous span
## assign 5*2 points from retina/vitreous to 30 into the vitreous.
## so, that's 200 points, slightly over-representing the outer retina span because that's where our hypothesis is focused.



END=dim(RAS.FLATTENED.LIGHT.RETINA.cropped)[1]
for(z in 1:dim(RAS.FLATTENED.LIGHT.RETINA.cropped)[3])
  {## the value tells me what dim[2] position in R.FLATTENED.LIGHT.RETINA.cropped to get the signal from::       as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-24,position.of.RPE.light,WHICH.INDEX))[258:END]
   AA=as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-30,position.of.RPE.light,WHICH.INDEX))[1:END];
   AB=as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-24,position.of.RPE.light,WHICH.INDEX))[1:END];
   AC=as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-18,position.of.RPE.light,WHICH.INDEX))[1:END];
   AD=as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-12,position.of.RPE.light,WHICH.INDEX))[1:END];
   AE=as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-6,position.of.RPE.light,WHICH.INDEX))[1:END];
   AF=as.vector(tapply(R.RPE.POSITION.LIGHT[,z],position.of.RPE.light,WHICH.INDEX))[1:END];
   BA=as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-37,position.of.RPE.light,WHICH.INDEX))[1:END];
   BB=as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-21,position.of.RPE.light,WHICH.INDEX))[1:END];
   BC=as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-15,position.of.RPE.light,WHICH.INDEX))[1:END];
   BD=as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-9,position.of.RPE.light,WHICH.INDEX))[1:END];
   BE=as.vector(tapply(R.RPE.POSITION.LIGHT[,z]-3,position.of.RPE.light,WHICH.INDEX))[1:END];
   BF=as.vector(tapply(R.RPE.POSITION.LIGHT[,z],position.of.RPE.light,WHICH.INDEX))[1:END];
   AG=as.vector(tapply(R.VITREOUS.RETINA.POSITION.LIGHT[,z]+6,position.of.RPE.light,WHICH.INDEX))[1:END];
   AH=as.vector(tapply(R.VITREOUS.RETINA.POSITION.LIGHT[,z]+12,position.of.RPE.light,WHICH.INDEX))[1:END];
   AI=as.vector(tapply(R.VITREOUS.RETINA.POSITION.LIGHT[,z]+18,position.of.RPE.light,WHICH.INDEX))[1:END];
   AJ=as.vector(tapply(R.VITREOUS.RETINA.POSITION.LIGHT[,z]+24,position.of.RPE.light,WHICH.INDEX))[1:END];
   AK=as.vector(tapply(R.VITREOUS.RETINA.POSITION.LIGHT[,z]+30,position.of.RPE.light,WHICH.INDEX))[1:END];
   BG=as.vector(tapply(R.VITREOUS.RETINA.POSITION.LIGHT[,z]+3,position.of.RPE.light,WHICH.INDEX))[1:END];
   BH=as.vector(tapply(R.VITREOUS.RETINA.POSITION.LIGHT[,z]+9,position.of.RPE.light,WHICH.INDEX))[1:END];
   BI=as.vector(tapply(R.VITREOUS.RETINA.POSITION.LIGHT[,z]+15,position.of.RPE.light,WHICH.INDEX))[1:END];
   BJ=as.vector(tapply(R.VITREOUS.RETINA.POSITION.LIGHT[,z]+21,position.of.RPE.light,WHICH.INDEX))[1:END];
   BK=as.vector(tapply(R.VITREOUS.RETINA.POSITION.LIGHT[,z]+27,position.of.RPE.light,WHICH.INDEX))[1:END];
   HARVEST=LIGHT.cropped[1:END,,z];
   for(x in 1:nrow(HARVEST))
    {RAS.FLATTENED.LIGHT.RETINA.cropped[(x),1,z]=HARVEST[x,AA[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),2,z]=HARVEST[x,BA[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),3,z]=HARVEST[x,AB[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),4,z]=HARVEST[x,BB[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),5,z]=HARVEST[x,AC[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),6,z]=HARVEST[x,BC[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),7,z]=HARVEST[x,AD[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),8,z]=HARVEST[x,BD[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),9,z]=HARVEST[x,AE[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),10,z]=HARVEST[x,BE[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),11,z]=HARVEST[x,AF[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),12,z]=HARVEST[x,BF[x]];
     #
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),191,z]=HARVEST[x,BG[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),192,z]=HARVEST[x,AG[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),193,z]=HARVEST[x,BH[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),194,z]=HARVEST[x,AH[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),195,z]=HARVEST[x,BI[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),196,z]=HARVEST[x,AI[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),197,z]=HARVEST[x,BJ[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),198,z]=HARVEST[x,AJ[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),199,z]=HARVEST[x,BK[x]];
     RAS.FLATTENED.LIGHT.RETINA.cropped[(x),200,z]=HARVEST[x,AK[x]]}
   ##
   ## "position.of.RPE.light" is just a strip of numbers; OK to use as factor here...
   STARTPOINT=R.RPE.POSITION.LIGHT[1:END,z];
   ENDPOINT=R.OLM.POSITION.LIGHT[1:END,z];
   list.index=as.factor(seq(1,38,1));
   for(x in 1:length(STARTPOINT)) 
     {if(!(is.na(STARTPOINT[x]+ENDPOINT[x]))) 
      {list=seq(STARTPOINT[x],ENDPOINT[x],(ENDPOINT[x]-STARTPOINT[x])/38)[2:39];
       RAS.FLATTENED.LIGHT.RETINA.cropped[(x),13:50,z]=HARVEST[x,as.vector(tapply(list,list.index,WHICH.INDEX))]}};
   STARTPOINT=R.OLM.POSITION.LIGHT[1:END,z];
   ENDPOINT=R.ONL.OPL.POSITION.LIGHT[1:END,z];
   list.index=as.factor(seq(1,36,1));
   for(x in 1:length(STARTPOINT)) 
     {if(!(is.na(STARTPOINT[x]+ENDPOINT[x]))) 
      {list=seq(STARTPOINT[x],ENDPOINT[x],(ENDPOINT[x]-STARTPOINT[x])/36)[2:37];
       RAS.FLATTENED.LIGHT.RETINA.cropped[(x),51:86,z]=HARVEST[x,as.vector(tapply(list,list.index,WHICH.INDEX))]}};
   STARTPOINT=R.ONL.OPL.POSITION.LIGHT[1:END,z];
   ENDPOINT=R.INL.IPL.POSITION.LIGHT[1:END,z];
   list.index=as.factor(seq(1,36,1));
   for(x in 1:length(STARTPOINT)) 
     {if(!(is.na(STARTPOINT[x]+ENDPOINT[x]))) 
      {list=seq(STARTPOINT[x],ENDPOINT[x],(ENDPOINT[x]-STARTPOINT[x])/36)[2:37];
       RAS.FLATTENED.LIGHT.RETINA.cropped[(x),87:122,z]=HARVEST[x,as.vector(tapply(list,list.index,WHICH.INDEX))]}};
   STARTPOINT=R.INL.IPL.POSITION.LIGHT[1:END,z];
   ENDPOINT=R.RNFL.GCL.POSITION.LIGHT[1:END,z];
   list.index=as.factor(seq(1,50,1));
   for(x in 1:length(STARTPOINT)) 
     {if(!(is.na(STARTPOINT[x]+ENDPOINT[x]))) 
      {list=seq(STARTPOINT[x],ENDPOINT[x],(ENDPOINT[x]-STARTPOINT[x])/50)[2:51];
       RAS.FLATTENED.LIGHT.RETINA.cropped[(x),123:172,z]=HARVEST[x,as.vector(tapply(list,list.index,WHICH.INDEX))]}};
   STARTPOINT=R.RNFL.GCL.POSITION.LIGHT[1:END,z];
   ENDPOINT=R.VITREOUS.RETINA.POSITION.LIGHT[1:END,z];
   list.index=as.factor(seq(1,18,1));
   for(x in 1:length(STARTPOINT)) 
     {if(!(is.na(STARTPOINT[x]+ENDPOINT[x]))) 
      {list=seq(STARTPOINT[x],ENDPOINT[x],(ENDPOINT[x]-STARTPOINT[x])/18)[2:19];
       RAS.FLATTENED.LIGHT.RETINA.cropped[(x),173:190,z]=HARVEST[x,as.vector(tapply(list,list.index,WHICH.INDEX))]}}}
## ignore the warnings



END=dim(RAS.FLATTENED.DARK.RETINA.cropped)[1]
for(z in 1:dim(RAS.FLATTENED.DARK.RETINA.cropped)[3])
  {## the value tells me what dim[2] position in R.FLATTENED.DARK.RETINA.cropped to get the signal from::       as.vector(tapply(R.RPE.POSITION.DARK[,z]-24,position.of.RPE.light,WHICH.INDEX))[258:END]
   AA=as.vector(tapply(R.RPE.POSITION.DARK[,z]-30,position.of.RPE.light,WHICH.INDEX))[1:END];
   AB=as.vector(tapply(R.RPE.POSITION.DARK[,z]-24,position.of.RPE.light,WHICH.INDEX))[1:END];
   AC=as.vector(tapply(R.RPE.POSITION.DARK[,z]-18,position.of.RPE.light,WHICH.INDEX))[1:END];
   AD=as.vector(tapply(R.RPE.POSITION.DARK[,z]-12,position.of.RPE.light,WHICH.INDEX))[1:END];
   AE=as.vector(tapply(R.RPE.POSITION.DARK[,z]-6,position.of.RPE.light,WHICH.INDEX))[1:END];
   AF=as.vector(tapply(R.RPE.POSITION.DARK[,z],position.of.RPE.light,WHICH.INDEX))[1:END];
   BA=as.vector(tapply(R.RPE.POSITION.DARK[,z]-37,position.of.RPE.light,WHICH.INDEX))[1:END];
   BB=as.vector(tapply(R.RPE.POSITION.DARK[,z]-21,position.of.RPE.light,WHICH.INDEX))[1:END];
   BC=as.vector(tapply(R.RPE.POSITION.DARK[,z]-15,position.of.RPE.light,WHICH.INDEX))[1:END];
   BD=as.vector(tapply(R.RPE.POSITION.DARK[,z]-9,position.of.RPE.light,WHICH.INDEX))[1:END];
   BE=as.vector(tapply(R.RPE.POSITION.DARK[,z]-3,position.of.RPE.light,WHICH.INDEX))[1:END];
   BF=as.vector(tapply(R.RPE.POSITION.DARK[,z],position.of.RPE.light,WHICH.INDEX))[1:END];
   AG=as.vector(tapply(R.VITREOUS.RETINA.POSITION.DARK[,z]+6,position.of.RPE.light,WHICH.INDEX))[1:END];
   AH=as.vector(tapply(R.VITREOUS.RETINA.POSITION.DARK[,z]+12,position.of.RPE.light,WHICH.INDEX))[1:END];
   AI=as.vector(tapply(R.VITREOUS.RETINA.POSITION.DARK[,z]+18,position.of.RPE.light,WHICH.INDEX))[1:END];
   AJ=as.vector(tapply(R.VITREOUS.RETINA.POSITION.DARK[,z]+24,position.of.RPE.light,WHICH.INDEX))[1:END];
   AK=as.vector(tapply(R.VITREOUS.RETINA.POSITION.DARK[,z]+30,position.of.RPE.light,WHICH.INDEX))[1:END];
   BG=as.vector(tapply(R.VITREOUS.RETINA.POSITION.DARK[,z]+3,position.of.RPE.light,WHICH.INDEX))[1:END];
   BH=as.vector(tapply(R.VITREOUS.RETINA.POSITION.DARK[,z]+9,position.of.RPE.light,WHICH.INDEX))[1:END];
   BI=as.vector(tapply(R.VITREOUS.RETINA.POSITION.DARK[,z]+15,position.of.RPE.light,WHICH.INDEX))[1:END];
   BJ=as.vector(tapply(R.VITREOUS.RETINA.POSITION.DARK[,z]+21,position.of.RPE.light,WHICH.INDEX))[1:END];
   BK=as.vector(tapply(R.VITREOUS.RETINA.POSITION.DARK[,z]+27,position.of.RPE.light,WHICH.INDEX))[1:END];
   HARVEST=DARK.cropped[1:END,,z];
   for(x in 1:nrow(HARVEST))
    {RAS.FLATTENED.DARK.RETINA.cropped[(x),1,z]=HARVEST[x,AA[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),2,z]=HARVEST[x,BA[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),3,z]=HARVEST[x,AB[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),4,z]=HARVEST[x,BB[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),5,z]=HARVEST[x,AC[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),6,z]=HARVEST[x,BC[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),7,z]=HARVEST[x,AD[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),8,z]=HARVEST[x,BD[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),9,z]=HARVEST[x,AE[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),10,z]=HARVEST[x,BE[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),11,z]=HARVEST[x,AF[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),12,z]=HARVEST[x,BF[x]];
     #
     RAS.FLATTENED.DARK.RETINA.cropped[(x),191,z]=HARVEST[x,BG[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),192,z]=HARVEST[x,AG[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),193,z]=HARVEST[x,BH[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),194,z]=HARVEST[x,AH[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),195,z]=HARVEST[x,BI[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),196,z]=HARVEST[x,AI[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),197,z]=HARVEST[x,BJ[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),198,z]=HARVEST[x,AJ[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),199,z]=HARVEST[x,BK[x]];
     RAS.FLATTENED.DARK.RETINA.cropped[(x),200,z]=HARVEST[x,AK[x]]}
   ##
   ## "position.of.RPE.light" is just a strip of numbers; OK to use as factor here...
   STARTPOINT=R.RPE.POSITION.DARK[1:END,z];
   ENDPOINT=R.OLM.POSITION.DARK[1:END,z];
   list.index=as.factor(seq(1,38,1));
   for(x in 1:length(STARTPOINT)) 
     {if(!(is.na(STARTPOINT[x]+ENDPOINT[x]))) 
      {list=seq(STARTPOINT[x],ENDPOINT[x],(ENDPOINT[x]-STARTPOINT[x])/38)[2:39];
       RAS.FLATTENED.DARK.RETINA.cropped[(x),13:50,z]=HARVEST[x,as.vector(tapply(list,list.index,WHICH.INDEX))]}};
   STARTPOINT=R.OLM.POSITION.DARK[1:END,z];
   ENDPOINT=R.ONL.OPL.POSITION.DARK[1:END,z];
   list.index=as.factor(seq(1,36,1));
   for(x in 1:length(STARTPOINT)) 
     {if(!(is.na(STARTPOINT[x]+ENDPOINT[x]))) 
      {list=seq(STARTPOINT[x],ENDPOINT[x],(ENDPOINT[x]-STARTPOINT[x])/36)[2:37];
       RAS.FLATTENED.DARK.RETINA.cropped[(x),51:86,z]=HARVEST[x,as.vector(tapply(list,list.index,WHICH.INDEX))]}};
   STARTPOINT=R.ONL.OPL.POSITION.DARK[1:END,z];
   ENDPOINT=R.INL.IPL.POSITION.DARK[1:END,z];
   list.index=as.factor(seq(1,36,1));
   for(x in 1:length(STARTPOINT)) 
     {if(!(is.na(STARTPOINT[x]+ENDPOINT[x]))) 
      {list=seq(STARTPOINT[x],ENDPOINT[x],(ENDPOINT[x]-STARTPOINT[x])/36)[2:37];
       RAS.FLATTENED.DARK.RETINA.cropped[(x),87:122,z]=HARVEST[x,as.vector(tapply(list,list.index,WHICH.INDEX))]}};
   STARTPOINT=R.INL.IPL.POSITION.DARK[1:END,z];
   ENDPOINT=R.RNFL.GCL.POSITION.DARK[1:END,z];
   list.index=as.factor(seq(1,50,1));
   for(x in 1:length(STARTPOINT)) 
     {if(!(is.na(STARTPOINT[x]+ENDPOINT[x]))) 
      {list=seq(STARTPOINT[x],ENDPOINT[x],(ENDPOINT[x]-STARTPOINT[x])/50)[2:51];
       RAS.FLATTENED.DARK.RETINA.cropped[(x),123:172,z]=HARVEST[x,as.vector(tapply(list,list.index,WHICH.INDEX))]}};
   STARTPOINT=R.RNFL.GCL.POSITION.DARK[1:END,z];
   ENDPOINT=R.VITREOUS.RETINA.POSITION.DARK[1:END,z];
   list.index=as.factor(seq(1,18,1));
   for(x in 1:length(STARTPOINT)) 
     {if(!(is.na(STARTPOINT[x]+ENDPOINT[x]))) 
      {list=seq(STARTPOINT[x],ENDPOINT[x],(ENDPOINT[x]-STARTPOINT[x])/18)[2:19];
       RAS.FLATTENED.DARK.RETINA.cropped[(x),173:190,z]=HARVEST[x,as.vector(tapply(list,list.index,WHICH.INDEX))]}}}
## ignore the warnings


RAS.FLATTENED.LIGHT.RETINA.cropped[which(is.na(RAS.FLATTENED.LIGHT.RETINA.cropped))]=0
RAS.FLATTENED.DARK.RETINA.cropped[which(is.na(RAS.FLATTENED.DARK.RETINA.cropped))]=0


#f.write.analyze(RAS.FLATTENED.LIGHT.RETINA.cropped[,,],paste(FILENAME,"_RAS--LIGHT",sep=""),size="float")
f.write.analyze(RAS.FLATTENED.DARK.RETINA.cropped[,,],paste(FILENAME,"_RAS--ALL",sep=""),size="float")


# z.Rdata
### this might be as good as it gets.
### output some of the thickness data:


Zz.MEAN.OLM.ONLopl.DIST.DARK=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)
Zz.MEAN.OLM.ONLopl.DIST.LIGHT=cbind(IMAGE.INDEX.LIGHT,IMAGE.INDEX.LIGHT)

for(x in 1:nrow(Zz.MEAN.OLM.ONLopl.DIST.DARK))
  {Zz.MEAN.OLM.ONLopl.DIST.DARK[x,2]=mean(na.rm=TRUE,(R.ONL.OPL.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]-R.OLM.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.OLM.ONLopl.DIST.LIGHT))
  {Zz.MEAN.OLM.ONLopl.DIST.LIGHT[x,2]=mean(na.rm=TRUE,(R.ONL.OPL.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]-R.OLM.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]))}

mean(Zz.MEAN.OLM.ONLopl.DIST.DARK[,2])
mean(Zz.MEAN.OLM.ONLopl.DIST.LIGHT[,2])
#t.test(Zz.MEAN.OLM.ONLopl.DIST.DARK[,2],Zz.MEAN.OLM.ONLopl.DIST.LIGHT[,2])


###
##

Zz.MEAN.ONLopl.INLipl.DIST.DARK=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)
Zz.MEAN.ONLopl.INLipl.DIST.LIGHT=cbind(IMAGE.INDEX.LIGHT,IMAGE.INDEX.LIGHT)

for(x in 1:nrow(Zz.MEAN.ONLopl.INLipl.DIST.DARK))
  {Zz.MEAN.ONLopl.INLipl.DIST.DARK[x,2]=mean(na.rm=TRUE,(R.INL.IPL.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]-R.ONL.OPL.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.ONLopl.INLipl.DIST.LIGHT))
  {Zz.MEAN.ONLopl.INLipl.DIST.LIGHT[x,2]=mean(na.rm=TRUE,(R.INL.IPL.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]-R.ONL.OPL.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]))}

mean(Zz.MEAN.ONLopl.INLipl.DIST.DARK[,2])
mean(Zz.MEAN.ONLopl.INLipl.DIST.LIGHT[,2])
#t.test(Zz.MEAN.ONLopl.INLipl.DIST.DARK[,2],Zz.MEAN.ONLopl.INLipl.DIST.LIGHT[,2])

###
##

Zz.MEAN.INLipl.GCLrnfl.DIST.DARK=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)
Zz.MEAN.INLipl.GCLrnfl.DIST.LIGHT=cbind(IMAGE.INDEX.LIGHT,IMAGE.INDEX.LIGHT)

for(x in 1:nrow(Zz.MEAN.INLipl.GCLrnfl.DIST.DARK))
  {Zz.MEAN.INLipl.GCLrnfl.DIST.DARK[x,2]=mean(na.rm=TRUE,(R.RNFL.GCL.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]-R.INL.IPL.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.INLipl.GCLrnfl.DIST.LIGHT))
  {Zz.MEAN.INLipl.GCLrnfl.DIST.LIGHT[x,2]=mean(na.rm=TRUE,(R.RNFL.GCL.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]-R.INL.IPL.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]))}

mean(Zz.MEAN.INLipl.GCLrnfl.DIST.DARK[,2])
mean(Zz.MEAN.INLipl.GCLrnfl.DIST.LIGHT[,2])
#t.test(Zz.MEAN.INLipl.GCLrnfl.DIST.DARK[,2],Zz.MEAN.INLipl.GCLrnfl.DIST.LIGHT[,2])


###
##

Zz.MEAN.GCLrnfl.RETvitDIST.DARK=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)
Zz.MEAN.GCLrnfl.RETvitDIST.LIGHT=cbind(IMAGE.INDEX.LIGHT,IMAGE.INDEX.LIGHT)

for(x in 1:nrow(Zz.MEAN.GCLrnfl.RETvitDIST.DARK))
  {Zz.MEAN.GCLrnfl.RETvitDIST.DARK[x,2]=mean(na.rm=TRUE,(R.VITREOUS.RETINA.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]-R.RNFL.GCL.POSITION.DARK[c(start.move:break.left,break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.GCLrnfl.RETvitDIST.LIGHT))
  {Zz.MEAN.GCLrnfl.RETvitDIST.LIGHT[x,2]=mean(na.rm=TRUE,(R.VITREOUS.RETINA.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]-R.RNFL.GCL.POSITION.LIGHT[c(start.move:break.left,break.right:end.move),x]))}

mean(Zz.MEAN.GCLrnfl.RETvitDIST.DARK[,2])
mean(Zz.MEAN.GCLrnfl.RETvitDIST.LIGHT[,2])
#t.test(Zz.MEAN.GCLrnfl.RETvitDIST.DARK[,2],Zz.MEAN.GCLrnfl.RETvitDIST.LIGHT[,2])



DARK.THICKNESSES=as.data.frame(matrix(,8,(2+nrow(Zz.MEAN.GCLrnfl.RETvitDIST.DARK))))
DARK.THICKNESSES[1,1]=FILENAME
DARK.THICKNESSES[1,2]="thicknesses_in_um"
DARK.THICKNESSES[1,3]="from_125_to_625um_from_OpticN"
DARK.THICKNESSES[2,1]="SPAN"
DARK.THICKNESSES[2,2]="MEAN_across_images"
for(x in 1:nrow(Zz.MEAN.GCLrnfl.RETvitDIST.DARK))
 {DARK.THICKNESSES[2,2+x]=paste("image_",Zz.MEAN.GCLrnfl.RETvitDIST.DARK[x,1],sep="")}
DARK.THICKNESSES[3,1]="WHOLE_RETINA"
DARK.THICKNESSES[4,1]="RNFL_THICK"
DARK.THICKNESSES[5,1]="INL_IPL_border_to_GCL_RNFL_border"
DARK.THICKNESSES[6,1]="ONL_OPL_border_to_INL_IPL_border"
DARK.THICKNESSES[7,1]="OLM_to_ONL_OPL_border"
DARK.THICKNESSES[8,1]="base_of_RPE_to_OLM"

DARK.THICKNESSES[3,2]=round(mean(na.rm=T,t(Zz.MEAN.RPE.VITREOUS.DIST.DARK[,2])),2)
DARK.THICKNESSES[3,3:ncol(DARK.THICKNESSES)]=round(t(Zz.MEAN.RPE.VITREOUS.DIST.DARK[,2]),2)

DARK.THICKNESSES[4,2]=round(mean(na.rm=T,t(Zz.MEAN.GCLrnfl.RETvitDIST.DARK[,2])),2)
DARK.THICKNESSES[4,3:ncol(DARK.THICKNESSES)]=round(t(Zz.MEAN.GCLrnfl.RETvitDIST.DARK[,2]),2)

DARK.THICKNESSES[5,2]=round(mean(na.rm=T,t(Zz.MEAN.INLipl.GCLrnfl.DIST.DARK[,2])),2)
DARK.THICKNESSES[5,3:ncol(DARK.THICKNESSES)]=round(t(Zz.MEAN.INLipl.GCLrnfl.DIST.DARK[,2]),2)

DARK.THICKNESSES[6,2]=round(mean(na.rm=T,t(Zz.MEAN.ONLopl.INLipl.DIST.DARK[,2])),2)
DARK.THICKNESSES[6,3:ncol(DARK.THICKNESSES)]=round(t(Zz.MEAN.ONLopl.INLipl.DIST.DARK[,2]),2)

DARK.THICKNESSES[7,2]=round(mean(na.rm=T,t(Zz.MEAN.OLM.ONLopl.DIST.DARK[,2])),2)
DARK.THICKNESSES[7,3:ncol(DARK.THICKNESSES)]=round(t(Zz.MEAN.OLM.ONLopl.DIST.DARK[,2]),2)

DARK.THICKNESSES[8,2]=round(mean(na.rm=T,t(Zz.MEAN.RPE.OML.DIST.DARK[,2])),2)
DARK.THICKNESSES[8,3:ncol(DARK.THICKNESSES)]=round(t(Zz.MEAN.RPE.OML.DIST.DARK[,2]),2)


LIGHT.THICKNESSES=as.data.frame(matrix(,8,(2+nrow(Zz.MEAN.GCLrnfl.RETvitDIST.LIGHT))))
LIGHT.THICKNESSES[1,1]=FILENAME
LIGHT.THICKNESSES[1,2]="thicknesses_in_um"
LIGHT.THICKNESSES[1,3]="from_125_to_625um_from_OpticN"
LIGHT.THICKNESSES[2,1]="SPAN"
LIGHT.THICKNESSES[2,2]="MEAN_across_images"
for(x in 1:nrow(Zz.MEAN.GCLrnfl.RETvitDIST.LIGHT))
 {LIGHT.THICKNESSES[2,2+x]=paste("image_",Zz.MEAN.GCLrnfl.RETvitDIST.LIGHT[x,1],sep="")}
LIGHT.THICKNESSES[3,1]="WHOLE_RETINA"
LIGHT.THICKNESSES[4,1]="RNFL_THICK"
LIGHT.THICKNESSES[5,1]="INL_IPL_border_to_GCL_RNFL_border"
LIGHT.THICKNESSES[6,1]="ONL_OPL_border_to_INL_IPL_border"
LIGHT.THICKNESSES[7,1]="OLM_to_ONL_OPL_border"
LIGHT.THICKNESSES[8,1]="base_of_RPE_to_OLM"

LIGHT.THICKNESSES[3,2]=round(mean(na.rm=T,t(Zz.MEAN.RPE.VITREOUS.DIST.LIGHT[,2])),2)
LIGHT.THICKNESSES[3,3:ncol(LIGHT.THICKNESSES)]=round(t(Zz.MEAN.RPE.VITREOUS.DIST.LIGHT[,2]),2)

LIGHT.THICKNESSES[4,2]=round(mean(na.rm=T,t(Zz.MEAN.GCLrnfl.RETvitDIST.LIGHT[,2])),2)
LIGHT.THICKNESSES[4,3:ncol(LIGHT.THICKNESSES)]=round(t(Zz.MEAN.GCLrnfl.RETvitDIST.LIGHT[,2]),2)

LIGHT.THICKNESSES[5,2]=round(mean(na.rm=T,t(Zz.MEAN.INLipl.GCLrnfl.DIST.LIGHT[,2])),2)
LIGHT.THICKNESSES[5,3:ncol(LIGHT.THICKNESSES)]=round(t(Zz.MEAN.INLipl.GCLrnfl.DIST.LIGHT[,2]),2)

LIGHT.THICKNESSES[6,2]=round(mean(na.rm=T,t(Zz.MEAN.ONLopl.INLipl.DIST.LIGHT[,2])),2)
LIGHT.THICKNESSES[6,3:ncol(LIGHT.THICKNESSES)]=round(t(Zz.MEAN.ONLopl.INLipl.DIST.LIGHT[,2]),2)

LIGHT.THICKNESSES[7,2]=round(mean(na.rm=T,t(Zz.MEAN.OLM.ONLopl.DIST.LIGHT[,2])),2)
LIGHT.THICKNESSES[7,3:ncol(LIGHT.THICKNESSES)]=round(t(Zz.MEAN.OLM.ONLopl.DIST.LIGHT[,2]),2)

LIGHT.THICKNESSES[8,2]=round(mean(na.rm=T,t(Zz.MEAN.RPE.OML.DIST.LIGHT[,2])),2)
LIGHT.THICKNESSES[8,3:ncol(LIGHT.THICKNESSES)]=round(t(Zz.MEAN.RPE.OML.DIST.LIGHT[,2]),2)

#write(t(LIGHT.THICKNESSES),ncol=ncol(LIGHT.THICKNESSES),file=paste(FILENAME,"_light_thicknesses.txt",sep=""))
#write(t(DARK.THICKNESSES),ncol=ncol(DARK.THICKNESSES),file=paste(FILENAME,"_all_thicknesses.txt",sep=""))



#################&&&


Zz.MEAN.RPE.OML.DIST.RIGHT=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)
Zz.MEAN.RPE.OML.DIST.LEFT=cbind(IMAGE.INDEX.LIGHT,IMAGE.INDEX.LIGHT)

for(x in 1:nrow(Zz.MEAN.RPE.OML.DIST.RIGHT))
  {Zz.MEAN.RPE.OML.DIST.RIGHT[x,2]=mean(na.rm=TRUE,(R.OLM.POSITION.DARK[c(break.right:end.move),x]-R.RPE.POSITION.DARK[c(break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.RPE.OML.DIST.LEFT))
  {Zz.MEAN.RPE.OML.DIST.LEFT[x,2]=mean(na.rm=TRUE,(R.OLM.POSITION.LIGHT[c(start.move:break.left),x]-R.RPE.POSITION.LIGHT[c(start.move:break.left),x]))}

mean(Zz.MEAN.RPE.OML.DIST.RIGHT[,2])
mean(Zz.MEAN.RPE.OML.DIST.LEFT[,2])
#t.test(Zz.MEAN.RPE.OML.DIST.DARK[,2],Zz.MEAN.RPE.OML.DIST.LIGHT[,2])



Zz.MEAN.RPE.VITREOUS.DIST.RIGHT=Zz.MEAN.RPE.OML.DIST.DARK
Zz.MEAN.RPE.VITREOUS.DIST.LEFT=Zz.MEAN.RPE.OML.DIST.LIGHT
for(x in 1:nrow(Zz.MEAN.RPE.VITREOUS.DIST.RIGHT))
  {Zz.MEAN.RPE.VITREOUS.DIST.RIGHT[x,2]=mean(na.rm=TRUE,(R.VITREOUS.RETINA.POSITION.DARK[c(break.right:end.move),x]-R.RPE.POSITION.DARK[c(break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.RPE.VITREOUS.DIST.LEFT))
  {Zz.MEAN.RPE.VITREOUS.DIST.LEFT[x,2]=mean(na.rm=TRUE,(R.VITREOUS.RETINA.POSITION.LIGHT[c(start.move:break.left),x]-R.RPE.POSITION.LIGHT[c(start.move:break.left),x]))}

mean(Zz.MEAN.RPE.VITREOUS.DIST.RIGHT[,2])
mean(Zz.MEAN.RPE.VITREOUS.DIST.LEFT[,2])
#t.test(Zz.MEAN.RPE.VITREOUS.DIST.DARK[,2],Zz.MEAN.RPE.VITREOUS.DIST.LIGHT[,2])


Zz.MEAN.OLM.ONLopl.DIST.RIGHT=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)
Zz.MEAN.OLM.ONLopl.DIST.LEFT=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)

for(x in 1:nrow(Zz.MEAN.OLM.ONLopl.DIST.RIGHT))
  {Zz.MEAN.OLM.ONLopl.DIST.RIGHT[x,2]=mean(na.rm=TRUE,(R.ONL.OPL.POSITION.DARK[c(break.right:end.move),x]-R.OLM.POSITION.DARK[c(break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.OLM.ONLopl.DIST.LEFT))
  {Zz.MEAN.OLM.ONLopl.DIST.LEFT[x,2]=mean(na.rm=TRUE,(R.ONL.OPL.POSITION.DARK[c(start.move:break.left),x]-R.OLM.POSITION.DARK[c(start.move:break.left),x]))}

mean(Zz.MEAN.OLM.ONLopl.DIST.RIGHT[,2])
mean(Zz.MEAN.OLM.ONLopl.DIST.LEFT[,2])
#t.test(Zz.MEAN.OLM.ONLopl.DIST.RIGHT[,2],Zz.MEAN.OLM.ONLopl.DIST.LEFT[,2])


###
##

Zz.MEAN.ONLopl.INLipl.DIST.RIGHT=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)
Zz.MEAN.ONLopl.INLipl.DIST.LEFT=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)

for(x in 1:nrow(Zz.MEAN.ONLopl.INLipl.DIST.RIGHT))
  {Zz.MEAN.ONLopl.INLipl.DIST.RIGHT[x,2]=mean(na.rm=TRUE,(R.INL.IPL.POSITION.DARK[c(break.right:end.move),x]-R.ONL.OPL.POSITION.DARK[c(break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.ONLopl.INLipl.DIST.LEFT))
  {Zz.MEAN.ONLopl.INLipl.DIST.LEFT[x,2]=mean(na.rm=TRUE,(R.INL.IPL.POSITION.DARK[c(start.move:break.left),x]-R.ONL.OPL.POSITION.DARK[c(start.move:break.left),x]))}

mean(Zz.MEAN.ONLopl.INLipl.DIST.RIGHT[,2])
mean(Zz.MEAN.ONLopl.INLipl.DIST.LEFT[,2])
#t.test(Zz.MEAN.ONLopl.INLipl.DIST.RIGHT[,2],Zz.MEAN.ONLopl.INLipl.DIST.LEFT[,2])

###
##

Zz.MEAN.INLipl.GCLrnfl.DIST.RIGHT=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)
Zz.MEAN.INLipl.GCLrnfl.DIST.LEFT=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)

for(x in 1:nrow(Zz.MEAN.INLipl.GCLrnfl.DIST.RIGHT))
  {Zz.MEAN.INLipl.GCLrnfl.DIST.RIGHT[x,2]=mean(na.rm=TRUE,(R.RNFL.GCL.POSITION.DARK[c(break.right:end.move),x]-R.INL.IPL.POSITION.DARK[c(break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.INLipl.GCLrnfl.DIST.LEFT))
  {Zz.MEAN.INLipl.GCLrnfl.DIST.LEFT[x,2]=mean(na.rm=TRUE,(R.RNFL.GCL.POSITION.DARK[c(start.move:break.left),x]-R.INL.IPL.POSITION.DARK[c(start.move:break.left),x]))}

mean(Zz.MEAN.INLipl.GCLrnfl.DIST.RIGHT[,2])
mean(Zz.MEAN.INLipl.GCLrnfl.DIST.LEFT[,2])
#t.test(Zz.MEAN.INLipl.GCLrnfl.DIST.RIGHT[,2],Zz.MEAN.INLipl.GCLrnfl.DIST.LEFT[,2])


###
##

Zz.MEAN.GCLrnfl.RETvitDIST.RIGHT=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)
Zz.MEAN.GCLrnfl.RETvitDIST.LEFT=cbind(IMAGE.INDEX.DARK,IMAGE.INDEX.DARK)

for(x in 1:nrow(Zz.MEAN.GCLrnfl.RETvitDIST.RIGHT))
  {Zz.MEAN.GCLrnfl.RETvitDIST.RIGHT[x,2]=mean(na.rm=TRUE,(R.VITREOUS.RETINA.POSITION.DARK[c(break.right:end.move),x]-R.RNFL.GCL.POSITION.DARK[c(break.right:end.move),x]))}

for(x in 1:nrow(Zz.MEAN.GCLrnfl.RETvitDIST.LEFT))
  {Zz.MEAN.GCLrnfl.RETvitDIST.LEFT[x,2]=mean(na.rm=TRUE,(R.VITREOUS.RETINA.POSITION.DARK[c(start.move:break.left),x]-R.RNFL.GCL.POSITION.DARK[c(start.move:break.left),x]))}

mean(Zz.MEAN.GCLrnfl.RETvitDIST.RIGHT[,2])
mean(Zz.MEAN.GCLrnfl.RETvitDIST.LEFT[,2])
#t.test(Zz.MEAN.GCLrnfl.RETvitDIST.RIGHT[,2],Zz.MEAN.GCLrnfl.RETvitDIST.LEFT[,2])



RIGHT.THICKNESSES=as.data.frame(matrix(,8,(2+nrow(Zz.MEAN.GCLrnfl.RETvitDIST.RIGHT))))
RIGHT.THICKNESSES[1,1]=FILENAME
RIGHT.THICKNESSES[1,2]="thicknesses_in_um"
RIGHT.THICKNESSES[1,3]="from_125_to_625um_from_OpticN"
RIGHT.THICKNESSES[2,1]="SPAN"
RIGHT.THICKNESSES[2,2]="MEAN_across_images"
for(x in 1:nrow(Zz.MEAN.GCLrnfl.RETvitDIST.RIGHT))
 {RIGHT.THICKNESSES[2,2+x]=paste("image_",Zz.MEAN.GCLrnfl.RETvitDIST.RIGHT[x,1],sep="")}
RIGHT.THICKNESSES[3,1]="WHOLE_RETINA"
RIGHT.THICKNESSES[4,1]="RNFL_THICK"
RIGHT.THICKNESSES[5,1]="INL_IPL_border_to_GCL_RNFL_border"
RIGHT.THICKNESSES[6,1]="ONL_OPL_border_to_INL_IPL_border"
RIGHT.THICKNESSES[7,1]="OLM_to_ONL_OPL_border"
RIGHT.THICKNESSES[8,1]="base_of_RPE_to_OLM"

RIGHT.THICKNESSES[3,2]=round(mean(na.rm=T,t(Zz.MEAN.RPE.VITREOUS.DIST.RIGHT[,2])),2)
RIGHT.THICKNESSES[3,3:ncol(RIGHT.THICKNESSES)]=round(t(Zz.MEAN.RPE.VITREOUS.DIST.RIGHT[,2]),2)

RIGHT.THICKNESSES[4,2]=round(mean(na.rm=T,t(Zz.MEAN.GCLrnfl.RETvitDIST.RIGHT[,2])),2)
RIGHT.THICKNESSES[4,3:ncol(RIGHT.THICKNESSES)]=round(t(Zz.MEAN.GCLrnfl.RETvitDIST.RIGHT[,2]),2)

RIGHT.THICKNESSES[5,2]=round(mean(na.rm=T,t(Zz.MEAN.INLipl.GCLrnfl.DIST.RIGHT[,2])),2)
RIGHT.THICKNESSES[5,3:ncol(RIGHT.THICKNESSES)]=round(t(Zz.MEAN.INLipl.GCLrnfl.DIST.RIGHT[,2]),2)

RIGHT.THICKNESSES[6,2]=round(mean(na.rm=T,t(Zz.MEAN.ONLopl.INLipl.DIST.RIGHT[,2])),2)
RIGHT.THICKNESSES[6,3:ncol(RIGHT.THICKNESSES)]=round(t(Zz.MEAN.ONLopl.INLipl.DIST.RIGHT[,2]),2)

RIGHT.THICKNESSES[7,2]=round(mean(na.rm=T,t(Zz.MEAN.OLM.ONLopl.DIST.RIGHT[,2])),2)
RIGHT.THICKNESSES[7,3:ncol(RIGHT.THICKNESSES)]=round(t(Zz.MEAN.OLM.ONLopl.DIST.RIGHT[,2]),2)

RIGHT.THICKNESSES[8,2]=round(mean(na.rm=T,t(Zz.MEAN.RPE.OML.DIST.RIGHT[,2])),2)
RIGHT.THICKNESSES[8,3:ncol(RIGHT.THICKNESSES)]=round(t(Zz.MEAN.RPE.OML.DIST.RIGHT[,2]),2)


LEFT.THICKNESSES=as.data.frame(matrix(,8,(2+nrow(Zz.MEAN.GCLrnfl.RETvitDIST.LEFT))))
LEFT.THICKNESSES[1,1]=FILENAME
LEFT.THICKNESSES[1,2]="thicknesses_in_um"
LEFT.THICKNESSES[1,3]="from_125_to_625um_from_OpticN"
LEFT.THICKNESSES[2,1]="SPAN"
LEFT.THICKNESSES[2,2]="MEAN_across_images"
for(x in 1:nrow(Zz.MEAN.GCLrnfl.RETvitDIST.LEFT))
 {LEFT.THICKNESSES[2,2+x]=paste("image_",Zz.MEAN.GCLrnfl.RETvitDIST.LEFT[x,1],sep="")}
LEFT.THICKNESSES[3,1]="WHOLE_RETINA"
LEFT.THICKNESSES[4,1]="RNFL_THICK"
LEFT.THICKNESSES[5,1]="INL_IPL_border_to_GCL_RNFL_border"
LEFT.THICKNESSES[6,1]="ONL_OPL_border_to_INL_IPL_border"
LEFT.THICKNESSES[7,1]="OLM_to_ONL_OPL_border"
LEFT.THICKNESSES[8,1]="base_of_RPE_to_OLM"

LEFT.THICKNESSES[3,2]=round(mean(na.rm=T,t(Zz.MEAN.RPE.VITREOUS.DIST.LEFT[,2])),2)
LEFT.THICKNESSES[3,3:ncol(LEFT.THICKNESSES)]=round(t(Zz.MEAN.RPE.VITREOUS.DIST.LEFT[,2]),2)

LEFT.THICKNESSES[4,2]=round(mean(na.rm=T,t(Zz.MEAN.GCLrnfl.RETvitDIST.LEFT[,2])),2)
LEFT.THICKNESSES[4,3:ncol(LEFT.THICKNESSES)]=round(t(Zz.MEAN.GCLrnfl.RETvitDIST.LEFT[,2]),2)

LEFT.THICKNESSES[5,2]=round(mean(na.rm=T,t(Zz.MEAN.INLipl.GCLrnfl.DIST.LEFT[,2])),2)
LEFT.THICKNESSES[5,3:ncol(LEFT.THICKNESSES)]=round(t(Zz.MEAN.INLipl.GCLrnfl.DIST.LEFT[,2]),2)

LEFT.THICKNESSES[6,2]=round(mean(na.rm=T,t(Zz.MEAN.ONLopl.INLipl.DIST.LEFT[,2])),2)
LEFT.THICKNESSES[6,3:ncol(LEFT.THICKNESSES)]=round(t(Zz.MEAN.ONLopl.INLipl.DIST.LEFT[,2]),2)

LEFT.THICKNESSES[7,2]=round(mean(na.rm=T,t(Zz.MEAN.OLM.ONLopl.DIST.LEFT[,2])),2)
LEFT.THICKNESSES[7,3:ncol(LEFT.THICKNESSES)]=round(t(Zz.MEAN.OLM.ONLopl.DIST.LEFT[,2]),2)

LEFT.THICKNESSES[8,2]=round(mean(na.rm=T,t(Zz.MEAN.RPE.OML.DIST.LEFT[,2])),2)
LEFT.THICKNESSES[8,3:ncol(LEFT.THICKNESSES)]=round(t(Zz.MEAN.RPE.OML.DIST.LEFT[,2]),2)

write(t(LEFT.THICKNESSES),ncol=ncol(LEFT.THICKNESSES),file=paste(FILENAME,"_left_thicknesses.txt",sep=""))
write(t(RIGHT.THICKNESSES),ncol=ncol(RIGHT.THICKNESSES),file=paste(FILENAME,"_right_thicknesses.txt",sep=""))








################&&&




## now, to get the profile means!

###
image(RAS.FLATTENED.LIGHT.RETINA.cropped[,,1])
##
## flip 0s back to NAs
RAS.FLATTENED.LIGHT.RETINA.cropped[which(RAS.FLATTENED.LIGHT.RETINA.cropped==0)]=NA
mean.RAS.FLATTENED.LIGHT.RETINA.cropped=matrix(,200,(1+dim(RAS.FLATTENED.LIGHT.RETINA.cropped)[3]))
for(z in 1:dim(RAS.FLATTENED.LIGHT.RETINA.cropped)[3])
 {image=RAS.FLATTENED.LIGHT.RETINA.cropped[,,z];
  for(x in 1:ncol(image))
  {mean.RAS.FLATTENED.LIGHT.RETINA.cropped[x,(z+1)]=mean(na.rm=T,image[,x])}}

perc.thick=round(c(seq(-1100/179,(-1*100/179),100/179),seq(0,100,100/179),seq((100+100/179),(100+1100/179),100/179)),2)
perc.thick=perc.thick[(length(perc.thick)-1):2]

mean.RAS.FLATTENED.LIGHT.RETINA.cropped[,1]=perc.thick

plot(mean.RAS.FLATTENED.LIGHT.RETINA.cropped[,1],mean.RAS.FLATTENED.LIGHT.RETINA.cropped[,2],col="white")
for(x in 2:ncol(mean.RAS.FLATTENED.LIGHT.RETINA.cropped)) matlines(mean.RAS.FLATTENED.LIGHT.RETINA.cropped[,1],mean.RAS.FLATTENED.LIGHT.RETINA.cropped[,x],col=grey(x/ncol(mean.RAS.FLATTENED.LIGHT.RETINA.cropped)),lty=x)

#write(t(mean.RAS.FLATTENED.LIGHT.RETINA.cropped[nrow(mean.RAS.FLATTENED.LIGHT.RETINA.cropped):1,]),ncol=ncol(mean.RAS.FLATTENED.LIGHT.RETINA.cropped),file=paste(FILENAME,"_light_profiles.txt",sep=""))

###
image(RAS.FLATTENED.DARK.RETINA.cropped[,,1])
##
## flip 0s back to NAs
RAS.FLATTENED.DARK.RETINA.cropped[which(RAS.FLATTENED.DARK.RETINA.cropped==0)]=NA
mean.RAS.FLATTENED.DARK.RETINA.cropped=matrix(,200,(1+dim(RAS.FLATTENED.DARK.RETINA.cropped)[3]))
for(z in 1:dim(RAS.FLATTENED.DARK.RETINA.cropped)[3])
 {image=RAS.FLATTENED.DARK.RETINA.cropped[,,z];
  for(x in 1:ncol(image))
  {mean.RAS.FLATTENED.DARK.RETINA.cropped[x,(z+1)]=mean(na.rm=T,image[,x])}}

perc.thick=round(c(seq(-1100/179,(-1*100/179),100/179),seq(0,100,100/179),seq((100+100/179),(100+1100/179),100/179)),2)
perc.thick=perc.thick[(length(perc.thick)-1):2]
mean.RAS.FLATTENED.DARK.RETINA.cropped[,1]=perc.thick

plot(mean.RAS.FLATTENED.DARK.RETINA.cropped[,1],mean.RAS.FLATTENED.DARK.RETINA.cropped[,2],col="white")
for(x in 2:ncol(mean.RAS.FLATTENED.DARK.RETINA.cropped)) matlines(mean.RAS.FLATTENED.DARK.RETINA.cropped[,1],mean.RAS.FLATTENED.DARK.RETINA.cropped[,x],col=grey(x/ncol(mean.RAS.FLATTENED.DARK.RETINA.cropped)),lty=x)

write(t(mean.RAS.FLATTENED.DARK.RETINA.cropped[nrow(mean.RAS.FLATTENED.DARK.RETINA.cropped):1,]),ncol=ncol(mean.RAS.FLATTENED.DARK.RETINA.cropped),file=paste(FILENAME,"_all_profiles_c.txt",sep=""))

save.image(paste(FILENAME,"_c.Rdata",sep=""))


peakSIG=max(na.rm=T,DARK.flat.cropped[,,1])
image(Zz.distances.from.opticN.microns,(seq(-19,160,1)*y.RESOLUTION),DARK.flat.cropped[,,1],col=gray.colors(33),main="Thick red = basement membrane")
##
abline(h=0*y.RESOLUTION,col="red",lwd=2)
matlines(Retina.Points.cropped[,4],R.VITREOUS.RETINA.POSITION.DARK[,1]-R.RPE.POSITION.DARK[,1],col="red")
matlines(Retina.Points.cropped[,4],R.RNFL.GCL.POSITION.DARK[,1]-R.RPE.POSITION.DARK[,1],col="blue")
matlines(Retina.Points.cropped[,4],R.ONL.OPL.POSITION.DARK[,1]-R.RPE.POSITION.DARK[,1],col="blue")
matlines(Retina.Points.cropped[,4],R.OLM.POSITION.DARK[,1]-R.RPE.POSITION.DARK[,1],col="red",lty=2)
matlines(Retina.Points.cropped[,4],R.INL.IPL.POSITION.DARK[,1]-R.RPE.POSITION.DARK[,1],col="red")
savePlot(filename = paste(FILENAME,"__",floor(runif(1000, 1,101))[1],".png",sep=""),
         type = "png")






R.RPE.POSITION.LIGHT

