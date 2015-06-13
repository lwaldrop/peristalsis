# Clears any previous data
rm(list=ls())

s1<- "short_RFreq0.5"
freq<-c(0.5)
num <-ceiling(10*freq)
dia<-(0.1)
rho <- 1000.0 #kg/m^3, 
#Pressure is non-dimensinalized by dividing by rho*(dia*freq), where freq = 1 Hz.
#Speed is non-dimensionalized by dividing by dia and multiplying by freq (=1)



# Sets working directory
dirname<-paste("/Volumes/LaCie_New/IBAMR/peristalsis_Wo1/no_racetrack/",s1,sep="")
setwd(dirname)

# Displays the current working directory.
getwd()


flowdata<- read.table("finalflow.csv", header=FALSE, sep=",")
names(flowdata)<-c("time","Vol.Flow.Rate","Max.Speed")
flowdata$Max.Speed<-flowdata$Max.Speed/dia
flowdata$Vol.Flow.Rate<-flowdata$Vol.Flow.Rate/dia^2

###### AVERAGE VOL FLOW RATE ######
var_avg<-"Non-dimensional Volume Flow Rate" #Change the variable name here to automatically plot it. 
#Makes a basic line plot of the data in the VisIt curve file.
name1<-"vol_flow_"
end<-".eps"
nameit<-paste(name1,s1,end,sep="")
setEPS()
postscript(nameit,width=6, height=4.5)
plot(Vol.Flow.Rate~time,data=flowdata,type="l",col="red",
xlab=list("Time",cex=1),ylab=list(var_avg,cex=1),ylim=c(min(flowdata$Vol.Flow.Rate),max(flowdata$Vol.Flow.Rate)))
dev.off()


###### MAX SPEED ######
var2<-"Non-dimensional Max Speed" #Change the variable name here to automatically plot it. 

#### MAX FLUID FLOW VC - FIND PEAKS ####

#Makes a basic line plot of the data in the VisIt curve file.
plot(Max.Speed~time,data=flowdata,type="l",col="red",
xlab=list("Time",cex=1),ylab=list(var2,cex=1),ylim=c(min(flowdata$Max.Speed),max(flowdata$Max.Speed)))

g2<-identify(x=flowdata$time,y=flowdata$Max.Speed,pos=FALSE,n=num)

#Uses the index numbers from the selected peaks (g2) to retrieve the time data (x) and mean data(y).
datamax.x<-flowdata$time[g2]
datamax.y<-flowdata$Max.Speed[g2]
#Plots the data again, this time plots the select points in red on top of the regular plot
plot(Max.Speed~time,data=flowdata,type="l",col="black",
xlab=list("Time (s)",cex=1),ylab=list(var2,cex=1),ylim=c(min(flowdata$Max.Speed),max(flowdata$Max.Speed)))
lines(datamax.x,datamax.y,type="o",lty=0,col="red",pch=21)


#### MEAN VALUES  #####
mean(flowdata$Max.Speed)	# average Max speed 
mean(datamax.y) 	# Mean peak speeds
mean(flowdata$Vol.Flow.Rate) 	# average Volume Flow Rate
