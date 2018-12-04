#This script is for the section 4.1 BUM simulation
#The p-values are from beta-uniform distribution
library(ggplot2)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(parallel)
#set the working directory
path=rstudioapi::getActiveDocumentContext()$path
setwd(dirname(path))

#Parallel computing(Otherwise the simulation would be too slow)
cl=makeCluster(detectCores()-1)

#The asymptotic critical level of the equal local level test
critical_asym<-function(n,alpha){
  -log(1-alpha)/(2*log(log(n))*log(n))
}

#compute the p-value for KS test
#parameter:
#P:P-value obtained from individual tests
ks.p<-function(p){
  ks=ks.test(p,"punif")
  return(ks$p.value)
}
#simulation study of the power of KS test
#The samples are simulated from Beta(a,b) distribution
#The proportion of Beta(a,b) is controlled by pi
#The simulation time is determined by simN
computePower_KS<-function(pi,simN,snpN,a,b){
  total=simN*snpN
  n.uniform=rbinom(1,total,1-pi)
  n.beta=total-n.uniform
  sim=c(rbeta(n.beta,a,b),runif(n.uniform))
  sim=sample(sim,total)
  sim=matrix(sim,simN,snpN)
  result=apply(sim,1,ks.p)
  power=mean(result<0.05)
  return(power)
}

#simulation study of the power of asym FHC test
#The samples are simulated from Beta(a,b) distribution
#The proportion of Beta(a,b) is controlled by pi
#The simulation time is determined by simN
computePower_asymFHC<-function(pi,simN,snpN,a,b){
  total=simN*snpN
  n.uniform=rbinom(1,total,1-pi)
  n.beta=total-n.uniform
  sim=c(rbeta(n.beta,a,b),runif(n.uniform))
  sim=sample(sim,total)
  gc()
  sim=matrix(sim,simN,snpN)
  gc()
  critical=critical_asym(snpN,0.05)
  sim=t(apply(sim,1,sort))
  gc()
  for(i in 1:snpN){
    sim[,i]=pbeta(sim[,i],i,snpN+1-i)
  }
  gc()
  hc.min.Ind=apply(sim,1,min)
  power=mean(hc.min.Ind<critical)
  return(power)
}


clusterExport(cl=cl,"ks.p")
clusterExport(cl=cl,"computePower_KS")
clusterExport(cl=cl,"critical_asym")
clusterExport(cl=cl,"computePower_asymFHC")


#The file name
parm=c("025","05","075","125","150","175")
#The parameter of the Beta distribution
#Alpha
a_list=c(0.25,0.5,0.75,1,1,1)
#Beta
b_list=c(1,1,1,1.25,1.5,1.75)
N_list=c(100,1000)
simN=10000

mydata=c()
set.seed(1)
for(snpN in N_list){
#Create the file name for reading
  f_list=paste0("data\\power_data",snpN,"_beta_",parm,".txt")
  
  
  
  #Read the power data and append the power for KS and asym FHC
  for(i in 1:length(a_list)){
    power_data=read.table(file=f_list[i],header = T)
    power_data$KS=NULL
    power_data$asymHC=NULL
    write.table(power_data,file=f_list[i], append = FALSE, sep = " ", dec = ".",
                row.names = FALSE, col.names = TRUE)
    pi=power_data$pi
    a=a_list[i]
    b=b_list[i]
    
    result=parSapplyLB(cl=cl,pi,computePower_KS,simN=simN,snpN=snpN,a=a,b=b)
    power_data$KS=result
    
    result1=parSapplyLB(cl=cl,pi,computePower_asymFHC,simN=simN,snpN=snpN,a=a,b=b)
    power_data$asymFHC=result1
    
    power_data$type=paste0("Beta(",a,",",b,")")
    
    power_data$N=snpN
    mydata=rbind(mydata,power_data)
    print(paste0("Beta(",a,",",b,")"))
  }
}



#===================================(N=1000)========================================
mydata1=mydata[mydata$N==1000,]
mydata1=melt(data=mydata1,id=c("pi","type"),measure=c("HC","FHC","KS","asymFHC"))

mydata1$color.group="HC"
mydata1$color.group[mydata1$variable=="FHC"]="Exact FHC"
mydata1$color.group[mydata1$variable=="KS"]="KS"
mydata1$color.group[mydata1$variable=="asymFHC"]="Asym FHC"

mydata1$color.group=factor(mydata1$color.group,levels=c("KS","HC","Asym FHC","Exact FHC"))

tick=c(0,0.05,0.25,0.5,0.75,1)

type_list=c("Beta(0.25,1)","Beta(0.5,1)","Beta(0.75,1)","Beta(1,1.25)","Beta(1,1.5)","Beta(1,1.75)")
xsize=c(0.2,0.25,1,1,1)
for(i in 1:length(type_list)){
  type=type_list[i]
  g=ggplot(data=mydata1[mydata1$type==type,],aes(x=pi,y=value,color=color.group))+
    labs(x=expression(pi),y="Power")+geom_line(linetype ="solid",alpha=1)+
    geom_hline(yintercept=0.05,linetype="dashed",alpha=0.5)+scale_y_continuous( breaks = tick, labels =tick,limits = c(0,1))+
    labs(color = "Statistic")+scale_x_continuous(limits = c(0,xsize[i]))+
    theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"),
          legend.position = c(0.99, 0.1),legend.justification = c(1, 0),legend.background = element_blank(),
          legend.box.background = element_rect(colour = "white"),legend.key.size = unit(0.4,"in"),axis.title=element_text(size=15),axis.text=element_text(size=12))+
    scale_color_manual(values=c("black", "blue", "green","red"))
  print(g+ggtitle(type))
  filename=gsub("[.(),]","",type,fixed=F)
  #ggsave(paste0("pvalueSim1000",filename,".png"), g,width=6.5,height=5.17,units = "in")
}



#===================================(N=100)========================================
mydata1=mydata[mydata$N==100,]
mydata1=melt(data=mydata1,id=c("pi","type"),measure=c("HC","FHC","KS","asymFHC"))

mydata1$color.group="HC"
mydata1$color.group[mydata1$variable=="FHC"]="Exact FHC"
mydata1$color.group[mydata1$variable=="KS"]="KS"
mydata1$color.group[mydata1$variable=="asymFHC"]="Asym FHC"

mydata1$color.group=factor(mydata1$color.group,levels=c("KS","HC","Asym FHC","Exact FHC"))

tick=c(0,0.05,0.25,0.5,0.75,1)
lgpos=matrix(c(0.99,0.1,0.01,0.99),2,2,byrow=T)
lgjust=matrix(c(1,0,0,1),2,2,byrow=T)
type_list=c("Beta(0.25,1)","Beta(0.5,1)","Beta(0.75,1)","Beta(1,1.25)","Beta(1,1.5)","Beta(1,1.75)")
xsize=c(1,1,1,1,1)
for(i in 1:length(type_list)){
  type=type_list[i]
  g=ggplot(data=mydata1[mydata1$type==type,],aes(x=pi,y=value,color=color.group))+
    labs(x=expression(pi),y="Power")+geom_line(linetype ="solid",alpha=1)+
    geom_hline(yintercept=0.05,linetype="dashed",alpha=0.5)+scale_y_continuous( breaks = tick, labels =tick,limits = c(0,1))+
    labs(color = "Statistic")+scale_x_continuous(limits = c(0,xsize[i]))+
    theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"),
          legend.position = lgpos[(i>=3)+1,],legend.justification = lgjust[(i>=3)+1,],legend.background = element_blank(),
          legend.box.background = element_rect(colour = "white"),legend.key.size = unit(0.4,"in"),axis.title=element_text(size=15),axis.text=element_text(size=12))+
    scale_color_manual(values=c("black", "blue", "green","red"))
  
  print(g+ggtitle(type))
  filename=gsub("[.(),]","",type,fixed=F)
  #ggsave(paste0("pvalueSim100",filename,".png"), g,width=6.5,height=5.17,units = "in")
}

