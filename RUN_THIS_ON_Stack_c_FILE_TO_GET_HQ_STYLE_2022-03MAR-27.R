###
### so, the retinal strip locked to the basement membrane is called DFC

image(DFC)

## and distances are stored in Zz.distances.from.opticN.microns

dim(DFC)
length(Zz.distances.from.opticN.microns)

## just use from +/- 350 to 630


LEFT.DFC=DFC[which(abs(Zz.distances.from.opticN.microns-(-630))==min(abs(Zz.distances.from.opticN.microns-(-630)))):
             which(abs(Zz.distances.from.opticN.microns-(-350))==min(abs(Zz.distances.from.opticN.microns-(-350)))),]

RIGHT.DFC=DFC[which(abs(Zz.distances.from.opticN.microns-(350))==min(abs(Zz.distances.from.opticN.microns-(350)))):
              which(abs(Zz.distances.from.opticN.microns-(630))==min(abs(Zz.distances.from.opticN.microns-(630)))),]


LEFT.DFC.profile=cbind((seq(-19,160,1)*y.RESOLUTION),(seq(-19,160,1)*y.RESOLUTION))
LEFT.DFC.profile[,2]<-NA
colnames(LEFT.DFC.profile)<-c("Distance.microns.from.BM","Signal")
RIGHT.DFC.profile=LEFT.DFC.profile
MEAN.DFC.profile=LEFT.DFC.profile

for(y in 1:ncol(LEFT.DFC))
 {LEFT.DFC.profile[y,2]=mean(na.rm=T,LEFT.DFC[,y]);
  RIGHT.DFC.profile[y,2]=mean(na.rm=T,RIGHT.DFC[,y]);
  MEAN.DFC.profile[y,2]=mean(na.rm=T,c(LEFT.DFC[,y],RIGHT.DFC[,y]))}

plot(LEFT.DFC.profile[,1],LEFT.DFC.profile[,2],type="l",ylim=c(0,250),main="data are anchored to the basement membrane (BM)", sub="left=black, right=red, mean(left&right)=blue",xlab="distance (microns) from basement membrane")
matlines(RIGHT.DFC.profile[,1],RIGHT.DFC.profile[,2],lwd=2,col="red")
matlines(MEAN.DFC.profile[,1],MEAN.DFC.profile[,2],lwd=2,col="blue")
abline(v=0)
text(0,250,"BM")
text(240,250,"vitreous")


savePlot(filename = paste(FILENAME,"_mean350to630profiles_",floor(runif(1000, 1,101))[1],".png",sep=""),
         type = "png")


RIGHT.LANDMARK.SET=NEW.LANDMARKS.by.col.position.flat.cropped[1:2,]
RIGHT.LANDMARK.SET[1,]<-NA
RIGHT.LANDMARK.SET[2,]<-NA
LEFT.LANDMARK.SET=RIGHT.LANDMARK.SET
MEAN.LANDMARK.SET=RIGHT.LANDMARK.SET

RIGHT.LANDMARK.SET.values=NEW.LANDMARKS.by.col.position.flat.cropped[which(abs(Zz.distances.from.opticN.microns-(350))==min(abs(Zz.distances.from.opticN.microns-(350)))):
                                                                     which(abs(Zz.distances.from.opticN.microns-(630))==min(abs(Zz.distances.from.opticN.microns-(630)))),]
LEFT.LANDMARK.SET.values=NEW.LANDMARKS.by.col.position.flat.cropped[which(abs(Zz.distances.from.opticN.microns-(-630))==min(abs(Zz.distances.from.opticN.microns-(-630)))):
                                                                    which(abs(Zz.distances.from.opticN.microns-(-350))==min(abs(Zz.distances.from.opticN.microns-(-350)))),]
RIGHT.LANDMARK.SET[1,1]=mean(na.rm=T,RIGHT.LANDMARK.SET.values[,1])
RIGHT.LANDMARK.SET[1,5]=mean(na.rm=T,RIGHT.LANDMARK.SET.values[,5])
RIGHT.LANDMARK.SET[1,6]=mean(na.rm=T,RIGHT.LANDMARK.SET.values[,6])
RIGHT.LANDMARK.SET[1,7]=mean(na.rm=T,RIGHT.LANDMARK.SET.values[,7])
RIGHT.LANDMARK.SET[1,8]=mean(na.rm=T,RIGHT.LANDMARK.SET.values[,8])
RIGHT.LANDMARK.SET[1,9]=mean(na.rm=T,RIGHT.LANDMARK.SET.values[,9])
RIGHT.LANDMARK.SET[1,10]=mean(na.rm=T,RIGHT.LANDMARK.SET.values[,10])
LEFT.LANDMARK.SET[1,1]=mean(na.rm=T,LEFT.LANDMARK.SET.values[,1])
LEFT.LANDMARK.SET[1,5]=mean(na.rm=T,LEFT.LANDMARK.SET.values[,5])
LEFT.LANDMARK.SET[1,6]=mean(na.rm=T,LEFT.LANDMARK.SET.values[,6])
LEFT.LANDMARK.SET[1,7]=mean(na.rm=T,LEFT.LANDMARK.SET.values[,7])
LEFT.LANDMARK.SET[1,8]=mean(na.rm=T,LEFT.LANDMARK.SET.values[,8])
LEFT.LANDMARK.SET[1,9]=mean(na.rm=T,LEFT.LANDMARK.SET.values[,9])
LEFT.LANDMARK.SET[1,10]=mean(na.rm=T,LEFT.LANDMARK.SET.values[,10])
MEAN.LANDMARK.SET=(RIGHT.LANDMARK.SET+LEFT.LANDMARK.SET)/2


### now, work through each and see if we can find a peak near OStips and the RPE.

RIGHT.TROUGH.EXTERIOR.TO.OLM=TROUGH.EXTERIOR.TO.OLM[1:2,]
RIGHT.TROUGH.EXTERIOR.TO.OLM[1,]<-NA
RIGHT.TROUGH.EXTERIOR.TO.OLM[2,]<-NA
RIGHT.ISOS=ISOS[1:2,]
RIGHT.ISOS[1,]<-NA
RIGHT.ISOS[2,]<-NA
RIGHT.GAP.JUST.INTERIOR.TO.RPE=GAP.JUST.INTERIOR.TO.RPE[1:2,]
RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,]<-NA
RIGHT.GAP.JUST.INTERIOR.TO.RPE[2,]<-NA
RIGHT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION=EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1:2,]
RIGHT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,]<-NA
RIGHT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[2,]<-NA


  x=1
  profile=cbind(RIGHT.DFC.profile[,1],seq(1,180,1),RIGHT.DFC.profile[,2]);
  profile=cbind(profile,profile[,3]-profile[c(2:nrow(profile),nrow(profile)),3]);
  plot(profile[,2],profile[,3])
  abline(v=(RIGHT.LANDMARK.SET[1,]),lwd=2)
  ## focus down:
  ## we know that none of the landmarks we're logging will be beyond the OLM nor the RPE PEAK
  profile2=profile[round(RIGHT.LANDMARK.SET[1,1]):round(RIGHT.LANDMARK.SET[1,10]),];
  matlines(profile2[,2],profile2[,3],col="red");
  profile2.spline=smooth.spline(profile2[,1],profile2[,3],df=nrow(profile2));
  profile2.spline.detail=seq(min(profile2[,1]),max(profile2[,1]),y.RESOLUTION/10);
  profile2.spline.detail=cbind(profile2.spline.detail,as.numeric(predict(profile2.spline,profile2.spline.detail)$y));
  matlines((profile2.spline.detail[,1]/y.RESOLUTION)+20,profile2.spline.detail[,2],col="red",lwd=2);
  ##
  ## first, get trough exterior to OLM
  if(!(is.na(RIGHT.LANDMARK.SET[1,9])))
    {search=profile[(RIGHT.LANDMARK.SET[1,9]-1):(RIGHT.LANDMARK.SET[1,9]+1),];
     RIGHT.TROUGH.EXTERIOR.TO.OLM[1,c(2,3,5)]=search[which(search[,3]==min(search[,3]))[1],c(2,1,3)];
     ## now, get ISOS (if present)(here, count as maximum in the proscribed area, even if not a true peak
     if( (!(is.na(RIGHT.LANDMARK.SET[1,7]))) | (!(is.na(RIGHT.LANDMARK.SET[1,8]))))
       {minny=min(na.rm=T,RIGHT.LANDMARK.SET[1,7:8]);
        maxxy=max(na.rm=T,RIGHT.LANDMARK.SET[1,7:8]);
        search=round(profile[(minny-1):(maxxy+1),]);
        RIGHT.ISOS[1,c(2,3,5)]=search[which(search[,3]==max(search[,3]))[1],c(2,1,3)];
        ## as long as I have a value for ISOS peak and position, I can now use half-height method to extract width of 
        ## trough exterior to OLM
        search.exterior=profile2.spline.detail[which(round(profile2.spline.detail[,1])==round(RIGHT.ISOS[1,3]))[1]:which(round(profile2.spline.detail[,1])==round(RIGHT.TROUGH.EXTERIOR.TO.OLM[1,3]))[length(which(round(profile2.spline.detail[,1])==round(RIGHT.TROUGH.EXTERIOR.TO.OLM[1,3])))],];
        search.interior=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(RIGHT.TROUGH.EXTERIOR.TO.OLM[1,3],2)):nrow(profile2.spline.detail),];
        #
        mid.ext=(max(na.rm=T,search.exterior[,2])+min(na.rm=T,search.exterior[,2]))/2;
        search.exterior[,2]=abs(search.exterior[,2]-mid.ext)
        ext.border=search.exterior[which(search.exterior[,2]==min(na.rm=T,search.exterior[,2]))[1],1];
        #
        mid.int=(max(na.rm=T,search.interior[,2])+min(na.rm=T,search.interior[,2]))/2;
        search.interior[,2]=abs(search.interior[,2]-mid.int)
        int.border=search.interior[which(search.interior[,2]==min(na.rm=T,search.interior[,2]))[1],1];
        ## 
        RIGHT.TROUGH.EXTERIOR.TO.OLM[x,4]=int.border-ext.border;
        ## to proceed, there needs to be an actual peak in the vicinity of the IS/OS 
        ## (expand search view, and look for at least one negative slope (walking from inside to outside the eye)
        ## (and at least one positive slope interior to that).
        search=profile[(minny-6):(maxxy+4),];
        find=which(search[,4]<0)
        ## 
        ## <new for 2022-MAR-27>: if there isn't a true peak in the vicinity of IS/OS, then we are forced to pick the place with the lowest slope
        if(length(find)==0) find=which(search[,4]==min(search[,4]))
        ## </new for 2022-MAR-27>: 
        if(length(find)>0)
          {seek=search[find[1]:nrow(search),];
           if(find[1]==nrow(search)) {seek=rbind(search[nrow(search),],search[nrow(search),]); seek[,1]<-NA};
           if(length(which(seek[,4]>0))>0)
             {## find the local minimum
              if( (!(is.na(RIGHT.LANDMARK.SET[1,5]))) | (!(is.na(RIGHT.LANDMARK.SET[1,6]))))
                {minny=min(na.rm=T,RIGHT.LANDMARK.SET[1,5:6]);
                 maxxy=max(na.rm=T,RIGHT.LANDMARK.SET[1,5:6]);
                 search=profile[(minny-1):(maxxy+1),];
                 THE.TROUGH=search[which(search[,3]==min(search[,3]))[1],];
                  #*
                 if(round(RIGHT.ISOS[1,3],2)==THE.TROUGH[1]) THE.TROUGH=search[which(search[,3]==min(search[,3]))[1]-1,];
                  #*
                 ## as long as I have a these values for the trough position and value, I can now use half-height method to extract width of 
                 ## ISOS
                 search.exteriorISOS=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(THE.TROUGH[1],2)):which(round(profile2.spline.detail[,1])==round(RIGHT.ISOS[1,3]))[1],];
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
                     #********************************************************
                 RIGHT.SAVE.FOR.LATER=cbind(search[,1],search[,3]);
                     #********************************************************

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
                    RIGHT.GAP.JUST.INTERIOR.TO.RPE[x,c(3,5)]=the.dip[c(1,2)];
                    RIGHT.GAP.JUST.INTERIOR.TO.RPE[x,6]=the.exterior.peak[3];
                    RIGHT.GAP.JUST.INTERIOR.TO.RPE[x,7]=the.interior.peak[3];
                    find.column=profile;
                    find.column[,1]<-abs(find.column[,1]-the.dip[1]);
                    RIGHT.GAP.JUST.INTERIOR.TO.RPE[x,2]=find.column[which(find.column[,1]==min(na.rm=T,find.column[,1]))[1],2];
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
                    RIGHT.GAP.JUST.INTERIOR.TO.RPE[x,8]=inner.dip.border;
                    RIGHT.GAP.JUST.INTERIOR.TO.RPE[x,9]=outer.dip.border;
                    RIGHT.GAP.JUST.INTERIOR.TO.RPE[x,4]=inner.dip.border-outer.dip.border;
                    RIGHT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[x,3]=search[1,1];
                    RIGHT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[x,2]=search[THEDIP+1,1]}}}}}}


LEFT.TROUGH.EXTERIOR.TO.OLM=TROUGH.EXTERIOR.TO.OLM[1:2,]
LEFT.TROUGH.EXTERIOR.TO.OLM[1,]<-NA
LEFT.TROUGH.EXTERIOR.TO.OLM[2,]<-NA
LEFT.ISOS=ISOS[1:2,]
LEFT.ISOS[1,]<-NA
LEFT.ISOS[2,]<-NA
LEFT.GAP.JUST.INTERIOR.TO.RPE=GAP.JUST.INTERIOR.TO.RPE[1:2,]
LEFT.GAP.JUST.INTERIOR.TO.RPE[1,]<-NA
LEFT.GAP.JUST.INTERIOR.TO.RPE[2,]<-NA
LEFT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION=EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1:2,]
LEFT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,]<-NA
LEFT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[2,]<-NA


  x=1
  profile=cbind(LEFT.DFC.profile[,1],seq(1,180,1),LEFT.DFC.profile[,2]);
  profile=cbind(profile,profile[,3]-profile[c(2:nrow(profile),nrow(profile)),3]);
  plot(profile[,2],profile[,3])
  abline(v=(LEFT.LANDMARK.SET[1,]),lwd=2)
  ## focus down:
  ## we know that none of the landmarks we're logging will be beyond the OLM nor the RPE PEAK
  profile2=profile[round(LEFT.LANDMARK.SET[1,1]):round(LEFT.LANDMARK.SET[1,10]),];
  matlines(profile2[,2],profile2[,3],col="red");
  profile2.spline=smooth.spline(profile2[,1],profile2[,3],df=nrow(profile2));
  profile2.spline.detail=seq(min(profile2[,1]),max(profile2[,1]),y.RESOLUTION/10);
  profile2.spline.detail=cbind(profile2.spline.detail,as.numeric(predict(profile2.spline,profile2.spline.detail)$y));
  matlines((profile2.spline.detail[,1]/y.RESOLUTION)+20,profile2.spline.detail[,2],col="red",lwd=2);
  ##
  ## first, get trough exterior to OLM
  if(!(is.na(LEFT.LANDMARK.SET[1,9])))
    {search=profile[(LEFT.LANDMARK.SET[1,9]-1):(LEFT.LANDMARK.SET[1,9]+1),];
     LEFT.TROUGH.EXTERIOR.TO.OLM[1,c(2,3,5)]=search[which(search[,3]==min(search[,3]))[1],c(2,1,3)];
     ## now, get ISOS (if present)(here, count as maximum in the proscribed area, even if not a true peak
     if( (!(is.na(LEFT.LANDMARK.SET[1,7]))) | (!(is.na(LEFT.LANDMARK.SET[1,8]))))
       {minny=min(na.rm=T,LEFT.LANDMARK.SET[1,7:8]);
        maxxy=max(na.rm=T,LEFT.LANDMARK.SET[1,7:8]);
        search=round(profile[(minny-1):(maxxy+1),]);
        LEFT.ISOS[1,c(2,3,5)]=search[which(search[,3]==max(search[,3]))[1],c(2,1,3)];
        ## as long as I have a value for ISOS peak and position, I can now use half-height method to extract width of 
        ## trough exterior to OLM
        search.exterior=profile2.spline.detail[which(round(profile2.spline.detail[,1])==round(LEFT.ISOS[1,3]))[1]:which(round(profile2.spline.detail[,1])==round(LEFT.TROUGH.EXTERIOR.TO.OLM[1,3]))[length(which(round(profile2.spline.detail[,1])==round(LEFT.TROUGH.EXTERIOR.TO.OLM[1,3])))],];
        search.interior=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(LEFT.TROUGH.EXTERIOR.TO.OLM[1,3],2)):nrow(profile2.spline.detail),];
        #
        mid.ext=(max(na.rm=T,search.exterior[,2])+min(na.rm=T,search.exterior[,2]))/2;
        search.exterior[,2]=abs(search.exterior[,2]-mid.ext)
        ext.border=search.exterior[which(search.exterior[,2]==min(na.rm=T,search.exterior[,2]))[1],1];
        #
        mid.int=(max(na.rm=T,search.interior[,2])+min(na.rm=T,search.interior[,2]))/2;
        search.interior[,2]=abs(search.interior[,2]-mid.int)
        int.border=search.interior[which(search.interior[,2]==min(na.rm=T,search.interior[,2]))[1],1];
        ## 
        LEFT.TROUGH.EXTERIOR.TO.OLM[x,4]=int.border-ext.border;
        ## to proceed, there needs to be an actual peak in the vicinity of the IS/OS 
        ## (expand search view, and look for at least one negative slope (walking from inside to outside the eye)
        ## (and at least one positive slope interior to that).
        search=profile[(minny-6):(maxxy+4),];
        find=which(search[,4]<0)
        ## <new for 2022-MAR-27>: if there isn't a true peak in the vicinity of IS/OS, then we are forced to pick the place with the lowest slope
        if(length(find)==0) find=which(search[,4]==min(search[,4]))
        ## </new for 2022-MAR-27>: 
        if(length(find)>0)
          {seek=search[find[1]:nrow(search),];
           if(find[1]==nrow(search)) {seek=rbind(search[nrow(search),],search[nrow(search),]); seek[,1]<-NA};
           if(length(which(seek[,4]>0))>0)
             {## find the local minimum
              if( (!(is.na(LEFT.LANDMARK.SET[1,5]))) | (!(is.na(LEFT.LANDMARK.SET[1,6]))))
                {minny=min(na.rm=T,LEFT.LANDMARK.SET[1,5:6]);
                 maxxy=max(na.rm=T,LEFT.LANDMARK.SET[1,5:6]);
                 search=profile[(minny-1):(maxxy+1),];
                 THE.TROUGH=search[which(search[,3]==min(search[,3]))[1],];
                  #*
                 if(round(LEFT.ISOS[1,3],2)==THE.TROUGH[1]) THE.TROUGH=search[which(search[,3]==min(search[,3]))[1]-1,];
                  #*
                 ## as long as I have a these values for the trough position and value, I can now use half-height method to extract width of 
                 ## ISOS
                 search.exteriorISOS=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(THE.TROUGH[1],2)):which(round(profile2.spline.detail[,1])==round(LEFT.ISOS[1,3]))[1],];
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
                     #********************************************************
                 LEFT.SAVE.FOR.LATER=cbind(search[,1],search[,3]);
                     #********************************************************

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
                    LEFT.GAP.JUST.INTERIOR.TO.RPE[x,c(3,5)]=the.dip[c(1,2)];
                    LEFT.GAP.JUST.INTERIOR.TO.RPE[x,6]=the.exterior.peak[3];
                    LEFT.GAP.JUST.INTERIOR.TO.RPE[x,7]=the.interior.peak[3];
                    find.column=profile;
                    find.column[,1]<-abs(find.column[,1]-the.dip[1]);
                    LEFT.GAP.JUST.INTERIOR.TO.RPE[x,2]=find.column[which(find.column[,1]==min(na.rm=T,find.column[,1]))[1],2];
                    interior.to.dip=dip.search[which(dip.search[,1]==the.dip[1]):nrow(dip.search),];
                    midway=(min(interior.to.dip[,2])+max(interior.to.dip[,2]))/2;
                    interior.to.dip[,2]<-abs(interior.to.dip[,2]-midway);
                    inner.dip.border=interior.to.dip[which(interior.to.dip[,2]==min(interior.to.dip[,2]))[1],1];
                    exterior.to.dip=dip.search[1:which(dip.search[,1]==the.dip[1]),];
			  ### <new for 2020_08AUG_12>
                    if(length(exterior.to.dip)==2) exterior.to.dip=rbind(exterior.to.dip,exterior.to.dip);
			  ### </new for 2020_08AUG_12>
			  ### plot(search[,1],search[,3])
			  ### matlines(exterior.to.dip[,1],exterior.to.dip[,2])
                    midway=(min(exterior.to.dip[,2])+max(exterior.to.dip[,2]))/2;
                    exterior.to.dip[,2]<-abs(exterior.to.dip[,2]-midway);
                    outer.dip.border=exterior.to.dip[which(exterior.to.dip[,2]==min(exterior.to.dip[,2]))[1],1];
                    LEFT.GAP.JUST.INTERIOR.TO.RPE[x,8]=inner.dip.border;
                    LEFT.GAP.JUST.INTERIOR.TO.RPE[x,9]=outer.dip.border;
                    LEFT.GAP.JUST.INTERIOR.TO.RPE[x,4]=inner.dip.border-outer.dip.border;
                    LEFT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[x,3]=search[1,1];
                    LEFT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[x,2]=search[THEDIP+1,1]}}}}}}


MEAN.TROUGH.EXTERIOR.TO.OLM=TROUGH.EXTERIOR.TO.OLM[1:2,]
MEAN.TROUGH.EXTERIOR.TO.OLM[1,]<-NA
MEAN.TROUGH.EXTERIOR.TO.OLM[2,]<-NA
MEAN.ISOS=ISOS[1:2,]
MEAN.ISOS[1,]<-NA
MEAN.ISOS[2,]<-NA
MEAN.GAP.JUST.INTERIOR.TO.RPE=GAP.JUST.INTERIOR.TO.RPE[1:2,]
MEAN.GAP.JUST.INTERIOR.TO.RPE[1,]<-NA
MEAN.GAP.JUST.INTERIOR.TO.RPE[2,]<-NA
MEAN.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION=EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1:2,]
MEAN.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,]<-NA
MEAN.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[2,]<-NA


  x=1
  profile=cbind(MEAN.DFC.profile[,1],seq(1,180,1),MEAN.DFC.profile[,2]);
  profile=cbind(profile,profile[,3]-profile[c(2:nrow(profile),nrow(profile)),3]);
  plot(profile[,2],profile[,3])
  abline(v=(MEAN.LANDMARK.SET[1,]),lwd=2)
  ## focus down:
  ## we know that none of the landmarks we're logging will be beyond the OLM nor the RPE PEAK
  profile2=profile[round(MEAN.LANDMARK.SET[1,1]):round(MEAN.LANDMARK.SET[1,10]),];
  matlines(profile2[,2],profile2[,3],col="red");
  profile2.spline=smooth.spline(profile2[,1],profile2[,3],df=nrow(profile2));
  profile2.spline.detail=seq(min(profile2[,1]),max(profile2[,1]),y.RESOLUTION/10);
  profile2.spline.detail=cbind(profile2.spline.detail,as.numeric(predict(profile2.spline,profile2.spline.detail)$y));
  matlines((profile2.spline.detail[,1]/y.RESOLUTION)+20,profile2.spline.detail[,2],col="red",lwd=2);
  ##
  ## first, get trough exterior to OLM
  if(!(is.na(MEAN.LANDMARK.SET[1,9])))
    {search=profile[(MEAN.LANDMARK.SET[1,9]-1):(MEAN.LANDMARK.SET[1,9]+1),];
     MEAN.TROUGH.EXTERIOR.TO.OLM[1,c(2,3,5)]=search[which(search[,3]==min(search[,3]))[1],c(2,1,3)];
     ## now, get ISOS (if present)(here, count as maximum in the proscribed area, even if not a true peak
     if( (!(is.na(MEAN.LANDMARK.SET[1,7]))) | (!(is.na(MEAN.LANDMARK.SET[1,8]))))
       {minny=min(na.rm=T,MEAN.LANDMARK.SET[1,7:8]);
        maxxy=max(na.rm=T,MEAN.LANDMARK.SET[1,7:8]);
        search=round(profile[(minny-1):(maxxy+1),]);
        MEAN.ISOS[1,c(2,3,5)]=search[which(search[,3]==max(search[,3]))[1],c(2,1,3)];
        ## as long as I have a value for ISOS peak and position, I can now use half-height method to extract width of 
        ## trough exterior to OLM
        search.exterior=profile2.spline.detail[which(round(profile2.spline.detail[,1])==round(MEAN.ISOS[1,3]))[1]:which(round(profile2.spline.detail[,1])==round(MEAN.TROUGH.EXTERIOR.TO.OLM[1,3]))[length(which(round(profile2.spline.detail[,1])==round(MEAN.TROUGH.EXTERIOR.TO.OLM[1,3])))],];
        search.interior=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(MEAN.TROUGH.EXTERIOR.TO.OLM[1,3],2)):nrow(profile2.spline.detail),];
        #
        mid.ext=(max(na.rm=T,search.exterior[,2])+min(na.rm=T,search.exterior[,2]))/2;
        search.exterior[,2]=abs(search.exterior[,2]-mid.ext)
        ext.border=search.exterior[which(search.exterior[,2]==min(na.rm=T,search.exterior[,2]))[1],1];
        #
        mid.int=(max(na.rm=T,search.interior[,2])+min(na.rm=T,search.interior[,2]))/2;
        search.interior[,2]=abs(search.interior[,2]-mid.int)
        int.border=search.interior[which(search.interior[,2]==min(na.rm=T,search.interior[,2]))[1],1];
        ## 
        MEAN.TROUGH.EXTERIOR.TO.OLM[x,4]=int.border-ext.border;
        ## to proceed, there needs to be an actual peak in the vicinity of the IS/OS 
        ## (expand search view, and look for at least one negative slope (walking from inside to outside the eye)
        ## (and at least one positive slope interior to that).
        search=profile[(minny-6):(maxxy+4),];
        find=which(search[,4]<0)
        ## <new for 2022-MAR-27>: if there isn't a true peak in the vicinity of IS/OS, then we are forced to pick the place with the lowest slope
        if(length(find)==0) find=which(search[,4]==min(search[,4]))
        ## </new for 2022-MAR-27>: 
        if(length(find)>0)
          {seek=search[find[1]:nrow(search),];
           if(find[1]==nrow(search)) {seek=rbind(search[nrow(search),],search[nrow(search),]); seek[,1]<-NA};
           if(length(which(seek[,4]>0))>0)
             {## find the local minimum
              if( (!(is.na(MEAN.LANDMARK.SET[1,5]))) | (!(is.na(MEAN.LANDMARK.SET[1,6]))))
                {minny=min(na.rm=T,MEAN.LANDMARK.SET[1,5:6]);
                 maxxy=max(na.rm=T,MEAN.LANDMARK.SET[1,5:6]);
                 search=profile[(minny-1):(maxxy+1),];
                 THE.TROUGH=search[which(search[,3]==min(search[,3]))[1],];
                  #*
                 if(round(MEAN.ISOS[1,3],2)==THE.TROUGH[1]) THE.TROUGH=search[which(search[,3]==min(search[,3]))[1]-1,];
                  #*
                 ## as long as I have a these values for the trough position and value, I can now use half-height method to extract width of 
                 ## ISOS
                 search.exteriorISOS=profile2.spline.detail[which(round(profile2.spline.detail[,1],2)==round(THE.TROUGH[1],2)):which(round(profile2.spline.detail[,1])==round(MEAN.ISOS[1,3]))[1],];
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
                     #********************************************************
                 MEAN.SAVE.FOR.LATER=cbind(search[,1],search[,3]);
                     #********************************************************

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
                    MEAN.GAP.JUST.INTERIOR.TO.RPE[x,c(3,5)]=the.dip[c(1,2)];
                    MEAN.GAP.JUST.INTERIOR.TO.RPE[x,6]=the.exterior.peak[3];
                    MEAN.GAP.JUST.INTERIOR.TO.RPE[x,7]=the.interior.peak[3];
                    find.column=profile;
                    find.column[,1]<-abs(find.column[,1]-the.dip[1]);
                    MEAN.GAP.JUST.INTERIOR.TO.RPE[x,2]=find.column[which(find.column[,1]==min(na.rm=T,find.column[,1]))[1],2];
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
                    MEAN.GAP.JUST.INTERIOR.TO.RPE[x,8]=inner.dip.border;
                    MEAN.GAP.JUST.INTERIOR.TO.RPE[x,9]=outer.dip.border;
                    MEAN.GAP.JUST.INTERIOR.TO.RPE[x,4]=inner.dip.border-outer.dip.border;
                    MEAN.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[x,3]=search[1,1];
                    MEAN.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[x,2]=search[THEDIP+1,1]}}}}}}



##################################
##################################
##################################

##################################
##################################
##################################

##################################
##################################
##################################

##################################
##################################
##################################

##################################
##################################
##################################



##################################
##################################
##################################

##################################
##################################
##################################

##################################
##################################
##################################

##################################
##################################
##################################

##################################
##################################
##################################






MEAN.GAP.JUST.INTERIOR.TO.RPE
MEAN.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION
MEAN.SAVE.FOR.LATER

LEFT.GAP.JUST.INTERIOR.TO.RPE
LEFT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION
LEFT.SAVE.FOR.LATER

RIGHT.GAP.JUST.INTERIOR.TO.RPE
RIGHT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION
RIGHT.SAVE.FOR.LATER

plot(MEAN.SAVE.FOR.LATER[,1],MEAN.SAVE.FOR.LATER[,2],ylim=c(80,220))
matlines(LEFT.SAVE.FOR.LATER[,1],LEFT.SAVE.FOR.LATER[,2],col="purple")
matlines(RIGHT.SAVE.FOR.LATER[,1],RIGHT.SAVE.FOR.LATER[,2],col="orange")



### OK, so, the "backup plan" if there is OSpeak identified, is to shave a couple points off of the SAVE.FOR.LATERs, 
### and try fitting a line anways. We'll identify the most-exterior trough, and move a step in from there.

CUT=max(MEAN.SAVE.FOR.LATER[,1])
if(CUT > max(LEFT.SAVE.FOR.LATER[,1])) CUT = max(LEFT.SAVE.FOR.LATER[,1])
if(CUT > max(RIGHT.SAVE.FOR.LATER[,1])) CUT = max(RIGHT.SAVE.FOR.LATER[,1])

MEAN.S=MEAN.SAVE.FOR.LATER
LEFT.S=LEFT.SAVE.FOR.LATER
RIGHT.S=RIGHT.SAVE.FOR.LATER

MEAN.S=MEAN.S[1:(which(MEAN.S[,1]==CUT)-1),]
LEFT.S=LEFT.S[1:(which(LEFT.S[,1]==CUT)-1),]
RIGHT.S=RIGHT.S[1:(which(RIGHT.S[,1]==CUT)-1),]

plot(MEAN.S[,1],MEAN.S[,2],ylim=c(80,220),pch=19,main="signal contours RPE (x=0) and trough between OStips and IS/OS",xlab="distance from RPE microns",ylab="Signal a.u.",sub="black=mean,orange=right,purp=left")
matpoints(LEFT.S[,1],LEFT.S[,2],col="purple",pch=19)
matpoints(RIGHT.S[,1],RIGHT.S[,2],col="orange",pch=19)
matlines(MEAN.S[,1],MEAN.S[,2],lwd=3)
matlines(LEFT.S[,1],LEFT.S[,2],col="purple",lwd=3)
matlines(RIGHT.S[,1],RIGHT.S[,2],col="orange",lwd=3)


INPUT=MEAN.S
### STEP#1
### use an iterative process where a linear best-fit is fit, the value with the lowest residual is removed, and then
### a new fit is made, on until there are only two points left.
INPUT2=INPUT
for(x in 1:(nrow(INPUT)-2))
 {theFIT=cbind(INPUT2[,c(1,2)],summary(lm(INPUT2[,2]~INPUT2[,1]))[3]$res);
  INPUT2=theFIT[-1*(which(theFIT[,3]==min(theFIT[,3]))[1]),]}


## <<new for 2021-FEB-08; is there is a value for GAP.JUST.INTERIOR.TO.RPE[1,8:9], require that one point remain on the interior to that dip, and one ##   point remain exterior>>
INPUT.GUIDE=mean(MEAN.GAP.JUST.INTERIOR.TO.RPE[1,8:9])
if(!(is.na(INPUT.GUIDE)))
 {INPUT2=INPUT;
  for(x in 1:(nrow(INPUT)-2))
  {theFIT=cbind(INPUT2[,c(1,2)],summary(lm(INPUT2[,2]~INPUT2[,1]))[3]$res);
   theFIT=cbind(theFIT,theFIT[,1]);
   theFIT[,4]<-ifelse(theFIT[,1]<INPUT.GUIDE,0,1);
   ## if there's only one value on each side, this is not to be removed
   if(length(which(theFIT[,4]==0))==1) theFIT[which(theFIT[,4]==0),3]=NA;
   if(length(which(theFIT[,4]==1))==1) theFIT[which(theFIT[,4]==1),3]=NA;
   INPUT2=theFIT[-1*(which(theFIT[,3]==min(na.rm=T,theFIT[,3]))[1]),]}}
## <</new for 2021-FEB-08>>


## <<new for 2021-FEB-16>>
#
# so, frustratingly, when there is no true local minimum, there are occassionally misses.
# when that happens, let's use the place with greatest (most positive) second derivative is the
if(is.na(INPUT.GUIDE))
 {WHERE=cbind(INPUT,c(INPUT[2:nrow(INPUT),2],NA),INPUT[,1]);
  WHERE[,4]<-WHERE[,3]-WHERE[,2];
  WHERE=cbind(INPUT,WHERE[,4],c(WHERE[2:nrow(INPUT),4],NA),INPUT[,1]);
  WHERE[,5]<-WHERE[,4]-WHERE[,3];
  ## and because of how we calculated, the "dip" is one row beyond the max value for the second derivative...
  INPUT.GUIDE=WHERE[(which(WHERE[,5]==max(na.rm=T,WHERE[,5]))[1]+1),1];
  rm(WHERE);
  INPUT2=INPUT;
  for(x in 1:(nrow(INPUT)-2))
  {theFIT=cbind(INPUT2[,c(1,2)],summary(lm(INPUT2[,2]~INPUT2[,1]))[3]$res);
   theFIT=cbind(theFIT,theFIT[,1]);
   theFIT[,4]<-ifelse(theFIT[,1]<INPUT.GUIDE,0,1);
   ## if there's only one value on each side, this is not to be removed
   if(length(which(theFIT[,4]==0))==1) theFIT[which(theFIT[,4]==0),3]=NA;
   if(length(which(theFIT[,4]==1))==1) theFIT[which(theFIT[,4]==1),3]=NA;
   INPUT2=theFIT[-1*(which(theFIT[,3]==min(na.rm=T,theFIT[,3]))[1]),]};
  rm(INPUT.GUIDE)}
## <</new for 2021-FEB-16>>


### STEP#2
### check and see if any refinement possible...
INPUT=cbind(INPUT,INPUT)
INPUT[,3]<- (summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[2,1])*INPUT[,1] + (summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[1,1])
INPUT[,4]<-INPUT[,2]-INPUT[,3]
###
### this should return with the value "integer(0)"
which(INPUT[,4]>0.1)
###
###

abline(lm(INPUT2[,2]~INPUT2[,1]))

### and the biggest devation from the straight line approximation is,,,
INPUT3=INPUT[c((-1)*(1:min(which(abs(round(INPUT[,4],2))<0.001))),(-1)*(max(which(abs(round(INPUT[,4],2))<0.001)):nrow(INPUT))),]
dip=which(INPUT3[,4]==min(INPUT3[,4]))


OUTPUT=as.data.frame(matrix(,13,15))
OUTPUT[1,1]=FILENAME
OUTPUT[1,2]="subretinal_space"
OUTPUT[1,5]="from_350_to_630um_from_OpticN"
OUTPUT[2,1]="______"
OUTPUT[2,2]="______"
OUTPUT[2,3]="______"
OUTPUT[2,4]="______"
OUTPUT[2,5]="______"
OUTPUT[2,6]="______"
OUTPUT[3,1]="next_lines:"
OUTPUT[3,2]="values_only_if_true signal_trough_idenfitied"
OUTPUT[4,1]="______"
OUTPUT[4,2]="RPE_peak_value"
OUTPUT[4,3]="RPE_peak_location_microns_from_BM"
OUTPUT[4,4]="OStips_peak_value"
OUTPUT[4,5]="OStips_peak_location_microns_from_BM"
OUTPUT[4,6]="trough_value"
OUTPUT[4,7]="trough_location_microns_from_BM"
OUTPUT[4,8]="trough_width_by_half_height_method"
OUTPUT[4,9]="peak_to_peak_line_estimate_at_trough"
OUTPUT[4,10]="difference_in_measured_value_minus_peak_to_peak_line_estimate_at_trough"
OUTPUT[4,11]="HQ_line_RPE_vicinity_value"
OUTPUT[4,12]="HQ_line_RPE_vicinity_location_microns_from_BM"
OUTPUT[4,13]="HQ_line_OStips_vicinity_value"
OUTPUT[4,14]="HQ_line_OStips_vicinity_location_microns_from_BM"
OUTPUT[4,15]="HQ_line_estimate_at_trough"
OUTPUT[4,16]="difference_in_measured_value_minus_HQ_line_estimate_at_trough"


OUTPUT[5,1]="left"
OUTPUT[6,1]="right"
OUTPUT[7,1]="mean"
OUTPUT[8,1]="______"
OUTPUT[8,2]="______"
OUTPUT[8,3]="______"
OUTPUT[8,4]="______"
OUTPUT[8,5]="______"
OUTPUT[8,6]="______"
OUTPUT[9,1]="next_lines:"
OUTPUT[9,2]="values_if_NO_true signal_trough_idenfitied"
OUTPUT[10,1]="______"
OUTPUT[10,2]="RPE_peak_value"
OUTPUT[10,3]="RPE_peak_location_microns_from_BM"
OUTPUT[10,4]="HQ_line_RPE_vicinity_value"
OUTPUT[10,5]="HQ_line_RPE_vicinity_location_microns_from_BM"
OUTPUT[10,6]="HQ_line_OStips_vicinity_value"
OUTPUT[10,7]="HQ_line_OStips_vicinity_location_microns_from_BM"
OUTPUT[10,8]="HQ_line_based_trough_estimate_value"
OUTPUT[10,9]="HQ_line_based_trough_estimate_location_microns_from_BM"
OUTPUT[10,10]="measured_value_at_HQ_line_based_trough_location_estimate"
OUTPUT[10,11]="difference_in_measured_value_minus_HQ_line_estimate_at_trough_location_estimate"
OUTPUT[11,1]="left"
OUTPUT[12,1]="right"
OUTPUT[13,1]="mean"


OUTPUT[7,11]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),2],2)
OUTPUT[7,12]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),1],2)

OUTPUT[7,13]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),2],2)
OUTPUT[7,14]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),1],2)

OUTPUT[13,5]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),1],2)
OUTPUT[13,4]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),2],2)

OUTPUT[13,7]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),1],2)
OUTPUT[13,6]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),2],2)

OUTPUT[13,9]=round(INPUT3[dip,1],2)
OUTPUT[13,8]=round(INPUT3[dip,3],2)
OUTPUT[13,10]=round(INPUT3[dip,2],2)
OUTPUT[13,11]=round(INPUT3[dip,4],2)

OUTPUT[13,2]=round(INPUT[1,2],2)
OUTPUT[13,3]=round(INPUT[1,1],2)

OUTPUT[7,2]=round(INPUT[1,2],2)
OUTPUT[7,3]=round(INPUT[1,1],2)



if(!(is.na(MEAN.GAP.JUST.INTERIOR.TO.RPE[1,2])))
 {newline=cbind(c(MEAN.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,3],MEAN.GAP.JUST.INTERIOR.TO.RPE[1,6]),
                c(MEAN.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,2],MEAN.GAP.JUST.INTERIOR.TO.RPE[1,7]));
  newline.param=summary(lm(newline[2,]~newline[1,]))$coef[,1];
  OUTPUT[7,2]=MEAN.GAP.JUST.INTERIOR.TO.RPE[1,6];
  OUTPUT[7,3]=MEAN.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,3];
  OUTPUT[7,4]=MEAN.GAP.JUST.INTERIOR.TO.RPE[1,7];
  OUTPUT[7,5]=MEAN.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,2];
  OUTPUT[7,6]=MEAN.GAP.JUST.INTERIOR.TO.RPE[1,5];
  OUTPUT[7,7]=MEAN.GAP.JUST.INTERIOR.TO.RPE[1,3];
  OUTPUT[7,8]=MEAN.GAP.JUST.INTERIOR.TO.RPE[1,4];
  OUTPUT[7,9]=newline.param[1]+newline.param[2]*MEAN.GAP.JUST.INTERIOR.TO.RPE[1,3];
  OUTPUT[7,10]=MEAN.GAP.JUST.INTERIOR.TO.RPE[1,5]-(newline.param[1]+newline.param[2]*MEAN.GAP.JUST.INTERIOR.TO.RPE[1,3]);
  OUTPUT[7,15]=summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[2,1]*MEAN.GAP.JUST.INTERIOR.TO.RPE[1,3] + summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[1,1]
  OUTPUT[7,16]=MEAN.GAP.JUST.INTERIOR.TO.RPE[1,5]-(summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[2,1]*MEAN.GAP.JUST.INTERIOR.TO.RPE[1,3] + summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[1,1])}










INPUT=RIGHT.S
### STEP#1
### use an iterative process where a linear best-fit is fit, the value with the lowest residual is removed, and then
### a new fit is made, on until there are only two points left.
INPUT2=INPUT
for(x in 1:(nrow(INPUT)-2))
 {theFIT=cbind(INPUT2[,c(1,2)],summary(lm(INPUT2[,2]~INPUT2[,1]))[3]$res);
  INPUT2=theFIT[-1*(which(theFIT[,3]==min(theFIT[,3]))[1]),]}



## <<new for 2021-FEB-08; is there is a value for GAP.JUST.INTERIOR.TO.RPE[1,8:9], require that one point remain on the interior to that dip, and one ##   point remain exterior>>
INPUT.GUIDE=mean(RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,8:9])
if(!(is.na(INPUT.GUIDE)))
 {INPUT2=INPUT;
  for(x in 1:(nrow(INPUT)-2))
  {theFIT=cbind(INPUT2[,c(1,2)],summary(lm(INPUT2[,2]~INPUT2[,1]))[3]$res);
   theFIT=cbind(theFIT,theFIT[,1]);
   theFIT[,4]<-ifelse(theFIT[,1]<INPUT.GUIDE,0,1);
   ## if there's only one value on each side, this is not to be removed
   if(length(which(theFIT[,4]==0))==1) theFIT[which(theFIT[,4]==0),3]=NA;
   if(length(which(theFIT[,4]==1))==1) theFIT[which(theFIT[,4]==1),3]=NA;
   INPUT2=theFIT[-1*(which(theFIT[,3]==min(na.rm=T,theFIT[,3]))[1]),]}}
## <</new for 2021-FEB-08>>



## <<new for 2021-FEB-16>>
#
# so, frustratingly, when there is no true local minimum, there are occassionally misses.
# when that happens, let's use the place with greatest (most positive) second derivative is the
if(is.na(INPUT.GUIDE))
 {WHERE=cbind(INPUT,c(INPUT[2:nrow(INPUT),2],NA),INPUT[,1]);
  WHERE[,4]<-WHERE[,3]-WHERE[,2];
  WHERE=cbind(INPUT,WHERE[,4],c(WHERE[2:nrow(INPUT),4],NA),INPUT[,1]);
  WHERE[,5]<-WHERE[,4]-WHERE[,3];
  ## and because of how we calculated, the "dip" is one row beyond the max value for the second derivative...
  INPUT.GUIDE=WHERE[(which(WHERE[,5]==max(na.rm=T,WHERE[,5]))[1]+1),1];
  rm(WHERE);
  INPUT2=INPUT;
  for(x in 1:(nrow(INPUT)-2))
  {theFIT=cbind(INPUT2[,c(1,2)],summary(lm(INPUT2[,2]~INPUT2[,1]))[3]$res);
   theFIT=cbind(theFIT,theFIT[,1]);
   theFIT[,4]<-ifelse(theFIT[,1]<INPUT.GUIDE,0,1);
   ## if there's only one value on each side, this is not to be removed
   if(length(which(theFIT[,4]==0))==1) theFIT[which(theFIT[,4]==0),3]=NA;
   if(length(which(theFIT[,4]==1))==1) theFIT[which(theFIT[,4]==1),3]=NA;
   INPUT2=theFIT[-1*(which(theFIT[,3]==min(na.rm=T,theFIT[,3]))[1]),]};
  rm(INPUT.GUIDE)}
## <</new for 2021-FEB-16>>


### STEP#2
### check and see if any refinement possible...
INPUT=cbind(INPUT,INPUT)
INPUT[,3]<- (summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[2,1])*INPUT[,1] + (summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[1,1])
INPUT[,4]<-INPUT[,2]-INPUT[,3]
###
### this should return with the value "integer(0)"
which(INPUT[,4]>0.1)
###
###




abline(lm(INPUT2[,2]~INPUT2[,1]),col="orange")

### and the biggest devation from the straight line approximation is,,,
INPUT3=INPUT[c((-1)*(1:min(which(abs(round(INPUT[,4],2))<0.001))),(-1)*(max(which(abs(round(INPUT[,4],2))<0.001)):nrow(INPUT))),]
dip=which(INPUT3[,4]==min(INPUT3[,4]))


OUTPUT[6,11]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),2],2)
OUTPUT[6,12]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),1],2)

OUTPUT[6,13]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),2],2)
OUTPUT[6,14]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),1],2)

OUTPUT[12,5]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),1],2)
OUTPUT[12,4]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),2],2)

OUTPUT[12,7]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),1],2)
OUTPUT[12,6]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),2],2)

OUTPUT[12,9]=round(INPUT3[dip,1],2)
OUTPUT[12,8]=round(INPUT3[dip,3],2)
OUTPUT[12,10]=round(INPUT3[dip,2],2)
OUTPUT[12,11]=round(INPUT3[dip,4],2)

OUTPUT[12,2]=round(INPUT[1,2],2)
OUTPUT[12,3]=round(INPUT[1,1],2)

OUTPUT[6,2]=round(INPUT[1,2],2)
OUTPUT[6,3]=round(INPUT[1,1],2)


if(!(is.na(RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,2])))
 {newline=cbind(c(RIGHT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,3],RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,6]),
                c(RIGHT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,2],RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,7]));
  newline.param=summary(lm(newline[2,]~newline[1,]))$coef[,1];
  OUTPUT[6,2]=RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,6];
  OUTPUT[6,3]=RIGHT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,3];
  OUTPUT[6,4]=RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,7];
  OUTPUT[6,5]=RIGHT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,2];
  OUTPUT[6,6]=RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,5];
  OUTPUT[6,7]=RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,3];
  OUTPUT[6,8]=RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,4];
  OUTPUT[6,9]=newline.param[1]+newline.param[2]*RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,3];
  OUTPUT[6,10]=RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,5]-(newline.param[1]+newline.param[2]*RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,3]);
  OUTPUT[6,15]=summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[2,1]*RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,3] + summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[1,1]
  OUTPUT[6,16]=RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,5]-(summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[2,1]*RIGHT.GAP.JUST.INTERIOR.TO.RPE[1,3] + summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[1,1])}







INPUT=LEFT.S
### STEP#1
### use an iterative process where a linear best-fit is fit, the value with the lowest residual is removed, and then
### a new fit is made, on until there are only two points left.
INPUT2=INPUT
for(x in 1:(nrow(INPUT)-2))
 {theFIT=cbind(INPUT2[,c(1,2)],summary(lm(INPUT2[,2]~INPUT2[,1]))[3]$res);
  INPUT2=theFIT[-1*(which(theFIT[,3]==min(theFIT[,3]))[1]),]}

## <<new for 2021-FEB-08; is there is a value for GAP.JUST.INTERIOR.TO.RPE[1,8:9], require that one point remain on the interior to that dip, and one ##   point remain exterior>>
INPUT.GUIDE=mean(LEFT.GAP.JUST.INTERIOR.TO.RPE[1,8:9])
if(!(is.na(INPUT.GUIDE)))
 {INPUT2=INPUT;
  for(x in 1:(nrow(INPUT)-2))
  {theFIT=cbind(INPUT2[,c(1,2)],summary(lm(INPUT2[,2]~INPUT2[,1]))[3]$res);
   theFIT=cbind(theFIT,theFIT[,1]);
   theFIT[,4]<-ifelse(theFIT[,1]<INPUT.GUIDE,0,1);
   ## if there's only one value on each side, this is not to be removed
   if(length(which(theFIT[,4]==0))==1) theFIT[which(theFIT[,4]==0),3]=NA;
   if(length(which(theFIT[,4]==1))==1) theFIT[which(theFIT[,4]==1),3]=NA;
   INPUT2=theFIT[-1*(which(theFIT[,3]==min(na.rm=T,theFIT[,3]))[1]),]}}

## <</new for 2021-FEB-08>>


## <<new for 2021-FEB-16>>
#
# so, frustratingly, when there is no true local minimum, there are occassionally misses.
# when that happens, let's use the place with greatest (most positive) second derivative is the
if(is.na(INPUT.GUIDE))
 {WHERE=cbind(INPUT,c(INPUT[2:nrow(INPUT),2],NA),INPUT[,1]);
  WHERE[,4]<-WHERE[,3]-WHERE[,2];
  WHERE=cbind(INPUT,WHERE[,4],c(WHERE[2:nrow(INPUT),4],NA),INPUT[,1]);
  WHERE[,5]<-WHERE[,4]-WHERE[,3];
  ## and because of how we calculated, the "dip" is one row beyond the max value for the second derivative...
  INPUT.GUIDE=WHERE[(which(WHERE[,5]==max(na.rm=T,WHERE[,5]))[1]+1),1];
  rm(WHERE);
  INPUT2=INPUT;
  for(x in 1:(nrow(INPUT)-2))
  {theFIT=cbind(INPUT2[,c(1,2)],summary(lm(INPUT2[,2]~INPUT2[,1]))[3]$res);
   theFIT=cbind(theFIT,theFIT[,1]);
   theFIT[,4]<-ifelse(theFIT[,1]<INPUT.GUIDE,0,1);
   ## if there's only one value on each side, this is not to be removed
   if(length(which(theFIT[,4]==0))==1) theFIT[which(theFIT[,4]==0),3]=NA;
   if(length(which(theFIT[,4]==1))==1) theFIT[which(theFIT[,4]==1),3]=NA;
   INPUT2=theFIT[-1*(which(theFIT[,3]==min(na.rm=T,theFIT[,3]))[1]),]};
  rm(INPUT.GUIDE)}
## <</new for 2021-FEB-16>>



### STEP#2
### check and see if any refinement possible...
INPUT=cbind(INPUT,INPUT)
INPUT[,3]<- (summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[2,1])*INPUT[,1] + (summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[1,1])
INPUT[,4]<-INPUT[,2]-INPUT[,3]
###
### this should return with the value "integer(0)"
which(INPUT[,4]>0.1)
###
###

abline(lm(INPUT2[,2]~INPUT2[,1]),col="purple")



### and the biggest devation from the straight line approximation is,,,
INPUT3=INPUT[c((-1)*(1:min(which(abs(round(INPUT[,4],2))<0.001))),(-1)*(max(which(abs(round(INPUT[,4],2))<0.001)):nrow(INPUT))),]
dip=which(INPUT3[,4]==min(INPUT3[,4]))


OUTPUT[5,11]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),2],2)
OUTPUT[5,12]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),1],2)

OUTPUT[5,13]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),2],2)
OUTPUT[5,14]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),1],2)

OUTPUT[11,5]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),1],2)
OUTPUT[11,4]=round(INPUT[min(which(abs(round(INPUT[,4],2))<0.001)),2],2)

OUTPUT[11,7]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),1],2)
OUTPUT[11,6]=round(INPUT[max(which(abs(round(INPUT[,4],2))<0.001)),2],2)

OUTPUT[11,9]=round(INPUT3[dip,1],2)
OUTPUT[11,8]=round(INPUT3[dip,3],2)
OUTPUT[11,10]=round(INPUT3[dip,2],2)
OUTPUT[11,11]=round(INPUT3[dip,4],2)

OUTPUT[11,2]=round(INPUT[1,2],2)
OUTPUT[11,3]=round(INPUT[1,1],2)

OUTPUT[5,2]=round(INPUT[1,2],2)
OUTPUT[5,3]=round(INPUT[1,1],2)


if(!(is.na(LEFT.GAP.JUST.INTERIOR.TO.RPE[1,2])))
 {newline=cbind(c(LEFT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,3],LEFT.GAP.JUST.INTERIOR.TO.RPE[1,6]),
                c(LEFT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,2],LEFT.GAP.JUST.INTERIOR.TO.RPE[1,7]));
  newline.param=summary(lm(newline[2,]~newline[1,]))$coef[,1];
  OUTPUT[5,2]=LEFT.GAP.JUST.INTERIOR.TO.RPE[1,6];
  OUTPUT[5,3]=LEFT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,3];
  OUTPUT[5,4]=LEFT.GAP.JUST.INTERIOR.TO.RPE[1,7];
  OUTPUT[5,5]=LEFT.EXTRA.STUFF.RPE.PEAK.POSITION.AND.LOCAL.MAX.OS.TIP.POSITION[1,2];
  OUTPUT[5,6]=LEFT.GAP.JUST.INTERIOR.TO.RPE[1,5];
  OUTPUT[5,7]=LEFT.GAP.JUST.INTERIOR.TO.RPE[1,3];
  OUTPUT[5,8]=LEFT.GAP.JUST.INTERIOR.TO.RPE[1,4];
  OUTPUT[5,9]=newline.param[1]+newline.param[2]*LEFT.GAP.JUST.INTERIOR.TO.RPE[1,3];
  OUTPUT[5,10]=LEFT.GAP.JUST.INTERIOR.TO.RPE[1,5]-(newline.param[1]+newline.param[2]*LEFT.GAP.JUST.INTERIOR.TO.RPE[1,3]);
  OUTPUT[5,15]=summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[2,1]*LEFT.GAP.JUST.INTERIOR.TO.RPE[1,3] + summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[1,1]
  OUTPUT[5,16]=LEFT.GAP.JUST.INTERIOR.TO.RPE[1,5]-(summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[2,1]*LEFT.GAP.JUST.INTERIOR.TO.RPE[1,3] + summary(lm(INPUT2[,2]~INPUT2[,1]))$coef[1,1])}

write(t(OUTPUT),ncol=ncol(OUTPUT),file=paste(FILENAME,"_mean_profile_subretinal_space_calc_.txt",sep=""),sep="\t")



savePlot(filename = paste("_mean_profile_subretinal_space_calc_",FILENAME,"__",floor(runif(1000, 1,101))[1],".png",sep=""),
         type = "png")









