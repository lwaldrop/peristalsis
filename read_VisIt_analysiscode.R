# Clears any previous data
rm(list=ls())


s1<- "short_JFreq2.0"
freq<-c(2.0)
num <-ceiling(10*freq)
dia<-(0.1)
rho <- 1000.0 #kg/m^3, 
#Pressure is non-dimensinalized by dividing by rho*(dia*freq), where freq = 1 Hz.
#Speed is non-dimensionalized by dividing by dia and multiplying by freq (=1)

###### LOADS EVERYTHING ######


# Sets working directory
dirname<-paste("/Volumes/LaCie_New/IBAMR/peristalsis_Wo1/shortracetrack/",s1,"/graphs/",sep="")
setwd(dirname)

# Displays the current working directory.
getwd()


### Vena Cava ###
place<-"vc"

measureit<-"avg.curve"
base<-paste(s1, place, measureit, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
datavcavg <- read.table(base, header=FALSE, sep="")
datavcavg$V3<-datavcavg$V2/dia
summary(datavcavg)

measureit<-"max.curve"
base<-paste(s1, place, measureit, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
datavcmax <- read.table(base, header=FALSE, sep="")
datavcmax$V3<-datavcmax$V2/dia
summary(datavcmax)

measureit<-"avgp.curve"
base<-paste(s1, place, measureit, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
datavcavgp <- read.table(base, header=FALSE, sep="")
datavcavgp$V3<-datavcavgp$V2/(dia^2*rho)
summary(datavcavgp)

### Aorta ###
place<-"aorta"

measureit<-"avg.curve"
base<-paste(s1, place, measureit, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
dataaortaavg <- read.table(base, header=FALSE, sep="")
dataaortaavg$V3<-dataaortaavg$V2/dia
summary(dataaortaavg)

measureit<-"max.curve"
base<-paste(s1, place, measureit, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
dataaortamax <- read.table(base, header=FALSE, sep="")
dataaortamax$V3<-dataaortamax$V2/dia
summary(dataaortamax)

measureit<-"avgp.curve"
base<-paste(s1, place, measureit, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
dataaortaavgp <- read.table(base, header=FALSE, sep="")
dataaortaavgp$V3<-dataaortaavgp$V2/(dia^2*rho)
summary(dataaortaavgp)

diff<-dataaortaavgp$V3-datavcavgp$V3


### Point speed ###
place<-"point"

measureit<-"Um.curve"
base<-paste(s1, place, measureit, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
datapointUm <- read.table(base, header=FALSE, sep="")
datapointUm$V3<-datapointUm$V2/dia
summary(datapointUm)

measureit<-"Ux.curve"
base<-paste(s1, place, measureit, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
datapointUx <- read.table(base, header=FALSE, sep="")
datapointUx$V3<-datapointUx$V2/dia
summary(datapointUx)

measureit<-"Uy.curve"
base<-paste(s1, place, measureit, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
datapointUy <- read.table(base, header=FALSE, sep="")
datapointUy$V3<-datapointUy$V2/(dia^2*rho)
summary(datapointUy)


###### AVERAGE SPEED ######
var_avg<-"Non-dimensional Average Speed" #Change the variable name here to automatically plot it. 
#Makes a basic line plot of the data in the VisIt curve file.
name1<-"avg_speed_"
end<-".eps"
nameit<-paste(name1,s1,end,sep="")
setEPS()
postscript(nameit,width=6, height=4.5)
plot(V3~V1,data=datavcavg,type="l",col="red",
xlab=list("Time",cex=1),ylab=list(var_avg,cex=1),ylim=c(min(dataaortaavg$V3),max(dataaortaavg$V3)))
lines(V3~V1,data=dataaortaavg,type="l",col="blue")
legend("topleft",legend=c("Vena cava","Aorta"),col=c("blue","red"),lty=c(1,1),pch=c(NA,NA))
dev.off()

###### MAX SPEED ######

var2<-"Non-dimensional Max Speed" #Change the variable name here to automatically plot it. 


#### MAX FLUID FLOW VC - FIND PEAKS ####

#Makes a basic line plot of the data in the VisIt curve file.
plot(V3~V1,data=datavcmax,type="l",col="red",
xlab=list("Time",cex=1),ylab=list(var2,cex=1),ylim=c(0,max(datavcmax$V3)))

g2<-identify(x=datavcmax$V1,y=datavcmax$V3,pos=FALSE,n=num)

#Uses the index numbers from the selected peaks (g2) to retrieve the time data (x) and mean data(y).
datavc.x<-datavcmax$V1[g2]
datavc.y<-datavcmax$V3[g2]
#Plots the data again, this time plots the select points in red on top of the regular plot
plot(V3~V1,data=datavcmax,type="l",col="black",
xlab=list("Time (s)",cex=1),ylab=list(var2,cex=1),ylim=c(0,max(datavcmax$V3)))
lines(datavc.x,datavc.y,type="o",lty=0,col="red",pch=21)



#### MAX FLUID FLOW AORTA - FIND PEAKS ####

#Makes a basic line plot of the data in the VisIt curve file.
plot(V3~V1,data=dataaortamax,type="l",col="red",
xlab=list("Time",cex=1),ylab=list(var2,cex=1),ylim=c(0,max(dataaortamax$V3)))

g2<-identify(x=dataaortamax$V1,y=dataaortamax$V3,pos=FALSE,n=num)

#Uses the index numbers from the selected peaks (g2) to retrieve the time data (x) and mean data(y).
dataaorta.x<-dataaortamax$V1[g2]
dataaorta.y<-dataaortamax$V3[g2]
#Plots the data again, this time plots the select points in red on top of the regular plot
plot(V3~V1,data=dataaortamax,type="l",col="black",
xlab=list("Time (s)",cex=1),ylab=list(var2,cex=1),ylim=c(0,max(dataaortamax$V3)))
lines(dataaorta.x,dataaorta.y,type="o",lty=0,col="red",pch=21)


name1<-"max_speed_"
end<-".eps"
nameit<-paste(name1,s1,end,sep="")
setEPS()
postscript(nameit,width=6, height=4.5)
plot(V3~V1,data=datavcmax,type="l",col="red",
xlab=list("Time",cex=1),ylab=list(var2,cex=1),ylim=c(min(dataaortamax$V3),max(dataaortamax$V3)))
lines(V3~V1,data=dataaortamax,type="l",col="blue")
legend("topleft",legend=c("Vena cava","Aorta"),col=c("blue","red"),lty=c(1,1),pch=c(NA,NA))
dev.off()


###### POINT SPEED ######


var3<-"Non-dimensional Speed" #Change the variable name here to automatically plot it. 
#Makes a basic line plot of the data in the VisIt curve file.
name1<-"point_speed_"
name2<-".eps"
nameit<-paste(name1,s1,name2,sep="")

setEPS()
postscript(nameit,width=6, height=4.5)
plot(V3~V1,data=datapointUm,type="l",col="red",
xlab=list("Time",cex=1),ylab=list(var3,cex=1),ylim=c(min(datapointUm$V3),max(datapointUm$V3)))
lines(V3~V1,data=datapointUx,col="blue")
lines(V3~V1,data=datapointUy,col="green")
legend("topright",legend=c("Um","Ux","Uy"),col=c("red","blue","green"),lty=c(1,1,1),pch=c(NA,NA,NA))
dev.off()



#### Pressure #####


var4<-"Non-dimensional Pressure" #Change the variable name here to automatically plot it. 
#Makes a basic line plot of the data in the VisIt curve file.
name1<-"pressure_"
name2<-".eps"
nameit<-paste(name1,s1,name2,sep="")
setEPS()
postscript(nameit,width=6, height=4.5)
plot(V3~V1,data=dataaortaavgp,type="l",col="red",
xlab=list("Time",cex=1),ylab=list(var4,cex=1),ylim=c(min(datavcavgp$V3),max(dataaortaavgp$V3)))
lines(V3~V1,data=datavcavgp,col="blue")
legend("bottomright",legend=c("Aorta","Vena cava"),col=c("red","blue"),lty=c(1,1))
dev.off()


var5<-"Non-dimensional Pressure difference" #Change the variable name here to automatically plot it. 
#Makes a basic line plot of the data in the VisIt curve file.
name1<-"pdiff_"
name2<-".eps"
nameit<-paste(name1,s1,name2,sep="")
setEPS()
postscript(nameit,width=6, height=4.5)

plot(dataaortaavgp$V1,diff,type="l",col="green",
xlab=list("Time",cex=1),ylab=list(var5,cex=1),ylim=c(min(diff),max(diff)))

dev.off()


#### MEAN VALUES  #####

mean(datavcavg$V3) 	# average speed Vena Cava
mean(datavcmax$V3)	# average max speed Vena Cava
max(datavcmax$V3)
mean(datavc.y) 	# Mean peak speeds, Vena Cava
mean(dataaortaavg$V3) 	# average speed Aorta
mean(dataaortamax$V3)	# average max speed Aorta
max(dataaortamax$V3)
mean(dataaorta.y) 	# Mean peak speeds, Aorta
mean(datapointUm$V3)	# point.Um
mean(datapointUx$V3)	# point.Ux
mean(dataaortaavgp$V3)	# aorta.avgp
mean(datavcavgp$V3)	# vc.avgp
mean(diff)		# pressure difference




