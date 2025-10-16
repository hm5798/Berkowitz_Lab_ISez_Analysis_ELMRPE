###

#image(RAS.FLATTENED.DARK.RETINA.cropped[,,1])
#play=RAS.FLATTENED.DARK.RETINA.cropped
#play[1:458,,1]=0
#image(play[,,1])

## ok, so 1:458 will average left, 
## 459:917 will average right

##
## flip 0s back to NAs
RAS.FLATTENED.DARK.RETINA.cropped[which(RAS.FLATTENED.DARK.RETINA.cropped==0)]=NA
mean.RAS.FLATTENED.DARK.RETINA.cropped.left=matrix(,200,(1+dim(RAS.FLATTENED.DARK.RETINA.cropped)[3]))
mean.RAS.FLATTENED.DARK.RETINA.cropped.right=mean.RAS.FLATTENED.DARK.RETINA.cropped.left
z=1
 image.left=RAS.FLATTENED.DARK.RETINA.cropped[1:458,,z];
 image.right=RAS.FLATTENED.DARK.RETINA.cropped[459:917,,z];
 for(x in 1:ncol(image))
  {mean.RAS.FLATTENED.DARK.RETINA.cropped.left[x,(z+1)]=mean(na.rm=T,image.left[,x]);
   mean.RAS.FLATTENED.DARK.RETINA.cropped.right[x,(z+1)]=mean(na.rm=T,image.right[,x])}

perc.thick=round(c(seq(-1100/179,(-1*100/179),100/179),seq(0,100,100/179),seq((100+100/179),(100+1100/179),100/179)),2)
perc.thick=perc.thick[(length(perc.thick)-1):2]
mean.RAS.FLATTENED.DARK.RETINA.cropped.left[,1]=perc.thick
mean.RAS.FLATTENED.DARK.RETINA.cropped.right[,1]=perc.thick


mean.RAS.FLATTENED.DARK.RETINA.cropped=mean.RAS.FLATTENED.DARK.RETINA.cropped.left
mean.RAS.FLATTENED.DARK.RETINA.cropped[,3]=mean.RAS.FLATTENED.DARK.RETINA.cropped.right[,2]

plot(mean.RAS.FLATTENED.DARK.RETINA.cropped[,1],mean.RAS.FLATTENED.DARK.RETINA.cropped[,2],col="white",sub="left is red")
matlines(mean.RAS.FLATTENED.DARK.RETINA.cropped[,1],mean.RAS.FLATTENED.DARK.RETINA.cropped[,2],col="red",lty=1)
matlines(mean.RAS.FLATTENED.DARK.RETINA.cropped[,1],mean.RAS.FLATTENED.DARK.RETINA.cropped[,3],col="blue",lty=1)

mean.RAS.FLATTENED.DARK.RETINA.cropped2=t(mean.RAS.FLATTENED.DARK.RETINA.cropped)
mean.RAS.FLATTENED.DARK.RETINA.cropped2=cbind(mean.RAS.FLATTENED.DARK.RETINA.cropped2[,1],mean.RAS.FLATTENED.DARK.RETINA.cropped2)
mean.RAS.FLATTENED.DARK.RETINA.cropped2[,1]=c("percent_thick","left","right")

write(mean.RAS.FLATTENED.DARK.RETINA.cropped2,ncol=ncol(mean.RAS.FLATTENED.DARK.RETINA.cropped),file=paste(FILENAME,"_left_and_right_profiles_c.txt",sep=""))
