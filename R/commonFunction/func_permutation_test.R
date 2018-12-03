#The permutation function, it will permutate cases and controls and pass it to the stat_func, 
#the stat_func is the function to compute the statistics
#parms:
#stat_func: a list of functions that compute statistics from one permutation, each function will give one and only one statistics
#parms: the parameters that will be passed to stat_func
#case: The case data
#control: The control data
#pivot: if the function needs to remove the mean effect of the samples before the permutation, 
#if yes, the mean of the samples will be set to 0
#simulate_time: The simulation times
#sample_replacement: If a sample can appear more than once in the permutation, if yes, it is equivalent to bootstrap
#chunk.size: an option to control the parallel computing, please do not modify it.
#Return: A simulate_time by number of stat_func matrix, each row represents a result from one permutation
permutation<-function(cl=NULL,stat_func,parms,case,control,pivot=FALSE,simulate_time=1000,sample_replacement=F,chunk.size=ceiling(simulate_time/length(cl))){
  func_num=length(stat_func)
  
  if(pivot){
    case=scale(case, scale=FALSE)
    control=scale(control, scale=FALSE)
  }
  #prepare the data
  total=rbind(case,control)
  caseN=dim(case)
  controlN=dim(control)
  totalN=dim(total)
  
  if(is.null(cl))
    result_permute=lapply( 1:simulate_time, sampling,stat_func=stat_func,parms=parms,data=total,totalN=totalN,caseN=caseN,type=sample_replacement)
  else
    result_permute=parLapply(cl, 1:simulate_time, sampling,stat_func=stat_func,parms=parms,data=total,totalN=totalN,caseN=caseN,type=sample_replacement)
  
  result_matrix=matrix(NA,simulate_time,func_num)
  for(i in 1:simulate_time){
    result_matrix[i,]=result_permute[[i]]
  }
  
  return(result_matrix)
}


#Compute the p-values for one permutation
#The simInd does not have any meaning, just an index to fit the requirement for the parLapply function
#parms:
#data: The cases and controls
#totalN:Total number of sampels
#caseN: The number of cases
#type: Same as sample_replacement, a parameter to switch between permutation and bootstrap
sampling<-function(simInd,stat_func,parms,data,totalN,caseN,type){
  
  index=sample(1:totalN[1],totalN[1],replace = type)
  case_ind=index[1:caseN[1]]
  control_ind=index[-(1:caseN[1])]
  case_permu=data[case_ind,]
  control_permu=data[control_ind,]
  p=computeP2(control_permu,case_permu)
  func_num=length(stat_func)
  result=rep(NA,func_num)
  
  for(i in 1:func_num){
    if(length(parms[[i]])==0)
      result[i]=stat_func[[i]](case_permu,control_permu,p)
    else{
      result[i]=stat_func[[i]](case_permu,control_permu,p,parms[[i]])
    }
  }
  return(result)
}

#Compute the p-values for the the permutation result
#parms:
#stats:The observed statistic
#permute_result: The permuted statistic
#type: either 'u' or 'd' to indicate upper probability of lower probability
get_pvalue<-function(stats,permute_result,type){
  type=(type=='u')-1
  stats_number=length(stats)
  if(stats_number==1)
    result=mean(permute_result>stats)+type
  else
    result=colMeans(sweep(permute_result, 2, stats, `-`)>0)+type
  result=abs(result)
  return(result)
}



#The permutation function, it will permutate cases and controls and pass it to the stat_func, 
#the stat_func is the function to compute the p-values
#parms:
#stat_func: the function computes p-values for one permutation, this function will return a vector of p-values,
#it is used to compute the p-values for the t-test
#parms: the parameters that will be passed to stat_func(Currently not useful)
#case: The case data
#control: The control data
#pivot: if the function needs to remove the mean effect of the samples before the permutation, 
#if yes, the mean of the samples will be set to 0
#simulate_time: The simulation times
#sample_replacement: If a sample can appear more than once in the permutation, if yes, it is equivalent to bootstrap
#chunk.size: an option to control the parallel computing, please do not modify it.
#Return: A simulate_time by number of variables in cases and controls, each row represents a result from one permutation
permutation_func<-function(cl=NULL,func,parms=c(),case,control,pivot=FALSE,simulate_time=200,sample_replacement=F,chunk.size=ceiling(simulate_time/length(cl))){
  if(pivot){
    case=scale(case, scale=FALSE)
    control=scale(control, scale=FALSE)
  }
  
  total=rbind(case,control)
  caseN=dim(case)
  controlN=dim(control)
  totalN=dim(total)
  
  if(is.null(cl))
    result_permute=lapply( 1:simulate_time, sampling_func,func=func,parms=parms,data=total,totalN=totalN[1],caseN=caseN[1],type=sample_replacement)
  else
    result_permute=parLapply(cl, 1:simulate_time, sampling_func,func=func,parms=parms,data=total,totalN=totalN[1],caseN=caseN[1],type=sample_replacement)
  result_matrix=matrix(NA,simulate_time,length(result_permute[[1]]))
  
  for(i in 1:simulate_time){
    result_matrix[i,]=result_permute[[i]]
  }
  
  return(result_matrix)
}


#Compute the p-values for one permutation
#The simInd does not have any meaning, just an index to fit the requirement for the parLapply function
#parms:
#data: The cases and controls
#totalN:Total number of sampels
#caseN: The number of cases
#type: Same as sample_replacement, a parameter to switch between permutation and bootstrap
sampling_func<-function(simInd,func,parms=parms,data,totalN,caseN,type){
  
  index=sample(1:totalN,totalN,replace = type)
  case_ind=index[1:caseN]
  control_ind=index[-(1:caseN)]
  case_permu=data[case_ind,]
  control_permu=data[control_ind,]
  p=computeP2(control_permu,case_permu)
  controlN=totalN-caseN
  if(length(parms)>0)
    result=func(control_permu,case_permu,p,parms)
  else
    result=func(control_permu,case_permu,p)
  return(result)
}


