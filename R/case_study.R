path=rstudioapi::getActiveDocumentContext()$path
setwd(dirname(path))
source("commonFunction\\func_statistics.r")
mydata=read.csv("data\\Yang_YoungOnsetHypertension_Illumina550_SBAS_rawpv.csv")


#======================plot===========================
library("ggplot2")
theGene="NAV2"
gene=mydata[mydata$Symbol==theGene,]
gene1=gene
gene1$CLR_N_BMI_pv=-log(gene1$CLR_N_BMI_pv)
critical1=-log(0.05)
critical2=-log(0.05/nrow(gene1))
g1=ggplot(gene1,aes(x=dbSNP_RS_ID,y=CLR_N_BMI_pv))+geom_point()+
  theme(axis.text.x = element_text(color="white"),panel.background = element_blank(),
        axis.line = element_line(colour = "black"),legend.key.size=unit(0.4,"in"),
        axis.title=element_text(size=15),axis.text=element_text(size=12))+
  labs(x="SNP",y="-log(P-value)")+
geom_hline(yintercept = critical1,linetype="dashed")+
  geom_hline(yintercept = critical2,linetype="twodash")+
scale_y_continuous(limits=c(0,critical2+1),breaks = c(critical1,critical2), labels =c("unadj","adj"))
g1
#ggsave("exampleManhattan.png", g1,width=6.5,height=5.17)

g2=ggplot(gene,aes(x=CLR_N_BMI_pv))+stat_ecdf(color="red")+
  geom_line(mapping=aes(x=x,y=y),data=data.frame(x=c(0,1),y=c(0,1)),linetype="dashed")+
  labs(x="x",y="Probability")+
  theme(panel.grid.major = element_line(colour = "grey80"),panel.background = element_blank(),
        axis.line = element_line(colour = "black"),legend.key.size=unit(0.4,"in"),
        axis.title=element_text(size=15),axis.text=element_text(size=12))+
  coord_cartesian(xlim = c(0, 1)) 
g2
#ggsave("exampleEcdf.png", g2,width=6.5,height=5.17)


#=========================Compute KS,HC,FHC statistics and p-values============================
theGene="NAV2"
pvalue=mydata[mydata$Symbol==theGene,"CLR_N_BMI_pv"]
#The KS statistic and p-value
p.ks=ks.test(pvalue,"punif")$p.value
#FHC statistic, the p-value of FHC needs to be computed from MATLAB
fhc=computeFHC_ind(0,0,pvalue)
fhc
#HC statistic, the p-value of HC needs to be computed from MATLAB
hc=computeHC(0,0,pvalue)
hc

#p-value for the asym local level test
snpN=length(pvalue)
localLevel=c()
sp=sort(pvalue)
for(i in 1:snpN){
  localLevel[i]=pbeta(sp[i],i,snpN+1-i)
}
p=1-exp(-2*log(log(snpN))*log(snpN)*min(localLevel))
p
