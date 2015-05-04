# Clears any previous data
rm(list=ls())

# Sets working directory
setwd("/Volumes/LaCie_New/IBAMR/peristalsis_Wo1/base_code/graphs")

# Displays the current working directory.
getwd()

s1<- "basecode"
freq<-c(1.0)

#### NOTE: In the Wo=1 set of simulations, the data analysis was turnicated so that many of the variables in this script (originally written for the Wo=10 sereis of simulations) were not analyzed because we anticipate not using them. These lines are commented out in case they will be used in the future! ####

###### AVERAGE SPEED ######

base<-"upper_avg.curve"
base1<-paste(s1, base, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
data1 <- read.table(base1, header=FALSE, sep="")
summary(data1)

var1<-"Average Speed (m/s)" #Change the variable name here to automatically plot it. 
#Makes a basic line plot of the data in the VisIt curve file.
name1<-"avg_speed_"
name2<-".eps"
nameit<-paste(name1,s1,name2,sep="")
setEPS()
postscript(nameit,width=6, height=4.5)
plot(V2~V1,data=data1,type="l",col="red",
xlab=list("Time (s)",cex=1),ylab=list(var1,cex=1),ylim=c(min(data1$V2),max(data1$V2)))
dev.off()

###### MAX SPEED ######


base<-"upper_max.curve"
base1<-paste(s1, base, sep = "_")
data21 <- read.table(base1, header=FALSE, sep="")
summary(data21)

var2<-"Max Speed (m/s)" #Change the variable name here to automatically plot it. 
#Makes a basic line plot of the data in the VisIt curve file.

plot(V2~V1,data=data21,type="l",col="black",
xlab=list("Time (s)",cex=1),ylab=list(var2,cex=1),ylim=c(min(data21$V2),max(data21$V2)))

num <-ceiling(10*freq)

#### MAX FLUID FLOW - FIND PEAKS ####
g2<-identify(x=data21$V1,y=data21$V2,pos=FALSE,n=num)

#Uses the index numbers from the selected peaks (g2) to retrieve the time data (x) and mean data(y).
data.x<-data21$V1[g2]
data.y<-data21$V2[g2]
#Plots the data again, this time plots the select points in red on top of the regular plot
plot(V2~V1,data=data21,type="l",col="black",
xlab=list("Time (s)",cex=1),ylab=list(var2,cex=1),ylim=c(min(data21$V2),max(data21$V2)))
lines(data.x,data.y,type="o",lty=0,col="red",pch=21)



#Makes a basic line plot of the data in the VisIt curve file.
name1<-"max_speed_"
name2<-".eps"
nameit<-paste(name1,s1,name2,sep="")
setEPS()
postscript(nameit,width=6, height=4.5)
plot(V2~V1,data=data21,type="l",col="red",
xlab=list("Time (s)",cex=1),ylab=list(var2,cex=1),ylim=c(min(data21$V2),max(data21$V2)))
dev.off()


###### POINT SPEED ######

base<-"point_Um.curve"
base31<-paste(s1, base, sep = "_")
base<-"point_Ux.curve"
base32<-paste(s1, base, sep = "_")
#Reads in data from the ".curve" files produced by VisIt.
data31 <- read.table(base31, header=FALSE, sep="")
summary(data31)
data32 <- read.table(base32, header=FALSE, sep="")
summary(data32)

var3<-"Speed (m/s)" #Change the variable name here to automatically plot it. 
#Makes a basic line plot of the data in the VisIt curve file.
name1<-"point_speed_"
name2<-".eps"
nameit<-paste(name1,s1,name2,sep="")
setEPS()
postscript(nameit,width=6, height=4.5)
plot(V2~V1,data=data31,type="l",col="red",
xlab=list("Time (s)",cex=1),ylab=list(var3,cex=1),ylim=c(min(data32$V2),max(data31$V2)))
lines(V2~V1,data=data32,col="blue")
legend("bottomright",legend=c(base31,base32),col=c("red","blue"),lty=c(1,1))
dev.off()



#### Pressure #####

base<-"aorta_avgp.curve"
base41<-paste(s1, base, sep = "_")
base<-"vc_avgp.curve"
base42<-paste(s1, base, sep = "_")
data41 <- read.table(base41, header=FALSE, sep="")
summary(data41)
data42 <- read.table(base42, header=FALSE, sep="")
summary(data42)

diff<-data41$V2-data42$V2

var4<-"Pressure (Pa)" #Change the variable name here to automatically plot it. 
#Makes a basic line plot of the data in the VisIt curve file.
name1<-"pressure_"
name2<-".eps"
nameit<-paste(name1,s1,name2,sep="")
setEPS()
postscript(nameit,width=6, height=4.5)
plot(V2~V1,data=data41,type="l",col="red",
xlab=list("Time (s)",cex=1),ylab=list(var4,cex=1),ylim=c(min(data42$V2),max(data41$V2)))
lines(V2~V1,data=data42,col="blue")
legend("bottomright",legend=c(base41,base42),col=c("red","blue"),lty=c(1,1))
dev.off()


var5<-"Pressure difference (Pa)" #Change the variable name here to automatically plot it. 
#Makes a basic line plot of the data in the VisIt curve file.
name1<-"pdiff_"
name2<-".eps"
nameit<-paste(name1,s1,name2,sep="")
setEPS()
postscript(nameit,width=6, height=4.5)

plot(data41$V1,diff,type="l",col="red",
xlab=list("Time (s)",cex=1),ylab=list(var5,cex=1),ylim=c(min(diff),max(diff)))

dev.off()


#### MEAN VALUES  #####

mean(data1$V2) 	# upper.avg
mean(data21$V2)	# upper.max
mean(data.y) 	# Mean peak speeds
mean(data31$V2)	# point.Um
mean(data32$V2)	# point.Ux
mean(data41$V2)	# aorta.avgp
mean(data42$V2)	# vc.avgp
mean(diff)		# pressure difference




