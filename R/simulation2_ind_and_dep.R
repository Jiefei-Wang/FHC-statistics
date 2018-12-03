#This script is for the section 4.2 correlated data simulation. 
#Both independent and dependent data are considered.
library(mvtnorm)
library(cobs)
library(reshape2)
library(ggplot2)
library(doParallel)
library(reshape2)
#set the working directory
path=rstudioapi::getActiveDocumentContext()$path
setwd(dirname(path))

#Read the necessary functions
source("commonFunction\\func.r")
source("commonFunction\\func_permutation_test.r")
source("commonFunction\\func_statistics.r")




#Parallel computing
cl=makeCluster(10)
doParallel::registerDoParallel(cl)

#set the parameters
simN=20000
snpN=100
sampleN=50
set.seed(1)

#Read the critical value for HC, FHC and asymptotic HC
creaticalList_HC=read.csv("commonFunction\\critical_value_HC.csv",header = T)
creaticalList_FHC=read.csv("commonFunction\\critical_value_fHC.csv",header = T)
creaticalList_asyHC=computeHCASYCV(snpN,0.05)


#The result will be stored in record
record=c()


#The effect size, number of significant value, number of CS block
effect.size_list=c(0,0.2)
num_effect_list=seq(0,snpN,4)[(snpN/4+1):1]
cov.block_list=c(snpN,snpN/4)

for(i1 in 1:length(effect.size_list)){
  for(i2 in 1:length(num_effect_list)){
    for(i3 in 1:length(cov.block_list)){
      effect.size=effect.size_list[i1]
      num_effect=num_effect_list[i2]
      if(sum(effect.size_list==0)!=0)
      if((effect.size==0&&num_effect!=0)||(effect.size!=0&&num_effect==0)){
        next
      }
      covMat=getCovMat(snpN,type="block_cs",parm=c(0.1,0.9,cov.block_list[i3]))
      
      packageList=c("mvtnorm","cobs")
      #Run the parallel computing to compute the power of HC, FHC
      result=foreach(icount(simN),.combine= resultCmb,.multicombine=TRUE,.inorder=FALSE,.packages =packageList) %dopar% {
        effect=rep(0,snpN)
        effect[sample(1:snpN,num_effect)]=effect.size
        control=rmvnorm(n = sampleN, rep(0,snpN), covMat)
        case=rmvnorm(n = sampleN, effect, covMat)
        
        #permute the p-value from the t-test
        permuted.p.list=permutation_func(func=computeP,case=case,control=control,simulate_time=1000,pivot=F,sample_replacement=F)
        permuted.sp.list=t(apply(permuted.p.list,1,sort))
        empirical.func=apply(permuted.sp.list,2,ecdf)
        
        #P-value from the t-test
        recordP=computeP(control,case)
        
        #Compute HC, KS, FHC statistic.
        recordHC=computeHC(case,control)
        recordKS=computeKS(case,control)
        recordFHC=computeFHC_permute(case,control,p=recordP,parms=empirical.func)

        stats=c(recordHC,recordKS,recordFHC)
        stat_func=c(computeHC,computeKS,computeFHC_permute)
        parms=list(list(),list(),empirical.func)
      
        #Permute the HC,KS,FHC statistics
        permute_result=permutation(stat_func=stat_func
                                   ,parms=parms
                                   ,case=case,control=control,simulate_time=1000,pivot = F,sample_replacement=F)
        
        #Get the permutation p-value from HC, FHC and KS stiatistics
        permuteP=get_pvalue(stats=stats
                            ,permute_result=permute_result,type=rep("u",ncol(permute_result)))
        
        #independent case
        #Compute the exact/asym p-value of KS,HC,FHC statistics
        ks=ks.test(recordP,"punif")
        exactKS=ks=ks$p.value
        
        exactFHC=computeFHC_ind(case,control,recordP)
        
        exactTest=c(KS=exactKS<0.05,Exact_HC=(recordHC>creaticalList_HC[snpN,2]),Exact_FHC=(exactFHC<creaticalList_FHC[snpN,2]))
        asyTest=c(recordHC>creaticalList_asyHC)
        list(recordHC=recordHC,recordFHC=recordFHC,permuteP=permuteP,exactTest=exactTest,asyTest=asyTest,recordP=recordP)
      }
      #Power of the permutation test
      permuteresult=colMeans(result$permuteP<0.05, na.rm = T)
      names(permuteresult)=c("Permute_HC","Permute_KS","Permute_FHC")
      #Power of the exact test
      exactResult=colMeans(result$exactTest)
      #Power of the asymptotic test
      asyResult=mean(result$asyTest)
      
      currentResult=c(effect.size=effect.size,num_effect=num_effect,covMat=cov.block_list[i3],permuteresult,exactResult,asyHCResult=asyResult)
      
      record=rbind(record,currentResult)

      print(currentResult)
    }
  }
}





#===================================plot========================================
record1=as.data.frame(record)


record2=melt(data=as.data.frame(record1),id=c("num_effect","effect.size","covMat"),
             measure=c("Permute_HC","Permute_FHC","Exact_HC","Exact_FHC"))

record2$color.group="Exact HC"
record2$color.group[record2$variable=="Permute_HC"]="Permute HC"
record2$color.group[record2$variable=="Exact_FHC"]="Exact FHC"
record2$color.group[record2$variable=="Permute_FHC"]="Permute FHC"

record2$color.group=factor(record2$color.group,levels=c("Exact HC","Permute HC","Exact FHC","Permute FHC"))
record2$num_effect=record2$num_effect/snpN
tick=c(0,0.05,0.25,0.5,0.75,1)
ind=(record2$effect.size==0.2|record2$num_effect==0)&record2$covMat==snpN
g1=ggplot(data=record2[ind,],aes(x=num_effect,y=value,color=color.group,linetype=color.group))+
  geom_line(alpha=0.5,size=1)+
  geom_hline(yintercept=0.05,alpha=0.5)+scale_y_continuous( breaks = tick, labels =tick,limits=c(0,1))+
  labs(x=expression(pi),y="Power")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"),
        legend.position = c(0.02, 1),legend.margin=margin(),legend.justification = c(0, 1),legend.background = element_blank(),
        legend.box.background = element_rect(colour = "white"),legend.text = element_text(size = 12),
        legend.key.size=unit(0.4,"in"),
        axis.title=element_text(size=15),axis.text=element_text(size=12))+
  scale_color_manual(name ="Statistic",values=c("red","red","black","black"))+
  scale_linetype_manual(name ="Statistic",values=c("solid","dashed","solid","dashed"))
g1
#ggsave(paste0("dataIndSim",".png"), g1,width=6,height=4)



ind=(record2$effect.size==0.2|record2$num_effect==0)&record2$covMat==snpN/4
g2=ggplot(data=record2[ind,],aes(x=num_effect,y=value,color=color.group,linetype=color.group))+
  geom_line(alpha=0.5,size=1)+
  geom_hline(yintercept=0.05,alpha=0.5)+scale_y_continuous( breaks = tick, labels =tick,limits=c(0,1))+
  labs(x=expression(pi),y="Power")+
  theme(plot.title = element_text(hjust = 0.5),panel.background = element_blank(),axis.line = element_line(colour = "black"),
        legend.position = c(0.02, 1),legend.margin=margin(),legend.justification = c(0, 1),legend.background = element_blank(),
        legend.box.background = element_rect(colour = "white"),legend.text = element_text(size = 12),
        legend.key.size=unit(0.4,"in"),
        axis.title=element_text(size=15),axis.text=element_text(size=12))+
  scale_color_manual(name ="Statistic",values=c("red","red","black","black"))+
  scale_linetype_manual(name ="Statistic",values=c("solid","dashed","solid","dashed"))                                                                          
g2
#ggsave(paste0("dataDepSim",".png"), g2,width=6,height=4)
