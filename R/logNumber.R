.logNum <- setClass("logNum", slots = c(sign = "integer", x = "ANY"))

logNum<-function(x){
  sign <- as.integer(sign(x))
  x = log(abs(x))
  .logNum(sign = sign, x= x)
}

logNumFromLogValue<-function(logX){
  sign <- 1L
  .logNum(sign = sign, x= logX)
}


logNumInternal <-function(sign, x){
  sign <- as.integer(sign)
  .logNum(sign = sign, x= x)
}


logNumAddOps<-function(e1,e2){
  xSign <- e1@sign
  ySign <- e2@sign
  if(xSign==0){
    return(e2)
  }
  if(ySign==0){
    return(e1)
  }
  x <- e1@x
  y <- e2@x
  if(xSign != ySign){
    if(xSign==-1){
      e1@sign = 1L
      return(e2-e1)
    }else{
      e2@sign = 1L
      return(e1-e2)
    }
  }
  if(x>y)
    z=logNumAdd(x,y)
  else
    z=logNumAdd(y,x)
  
  logNumInternal(xSign, z)
}


logNumMinusOps <-function(e1,e2){
  xSign <- e1@sign
  ySign <- -e2@sign
  if(xSign==0){
    e2@sign = -e2@sign
    return(e2)
  }
  if(ySign==0){
    return(e1)
  }
  x <- e1@x
  y <- e2@x
  if(xSign == ySign){
    e2@sign = -e2@sign
    return(e1+e2)
  }
  
  if(x>=y){
    z <- logNumMinus(x,y)
    zSign <- xSign
  }else{
    z <- logNumMinus(y,x)
    zSign <- -xSign
  }
  logNumInternal(zSign, z)
}

logNumMultipleOps <-function(e1,e2){
  zSign <- e1@sign*e2@sign
  z <- e1@x + e2@x
  logNumInternal(zSign, z)
}
logNumDivisionOps <-function(e1,e2){
  if(ySign==0){
    e1@x = - e2@x
    return(e1)
  }
  zSign <- e1@sign*e2@sign
  z <- e1@x - e2@x
  logNumInternal(zSign, z)
}

logNumPowerOps<- function(e1,e2){
  e1@x <- e1@x * e2
  e1
}

showLogNum <- function(object) {
  xSign= object@sign
  x <- format(object@x, digits=5)
  value <- format(as(object,"numeric"), digits=5)
  if(xSign == 1||xSign==0)
    signSymbol <- "+"
  else
    signSymbol <- "-"
    
  cat(paste0("Log transformed number. Sign:",signSymbol,"\n"))
  cat(paste0("Log Value: ", x, " , Real value: ", value))
}




logNumAdd <- function(x,y){
  x + log(exp(y-x)+1)
}


logNumMinus <- function(x,y){
  x + log(1 - exp(y-x))
}


setMethod("+", signature = c(e1 = "logNum", e2 = "logNum"),logNumAddOps)
setMethod("-", signature = c(e1 = "logNum", e2 = "logNum"),logNumMinusOps)
setMethod("*", signature = c(e1 = "logNum", e2 = "logNum"),logNumMultipleOps)
setMethod("/", signature = c(e1 = "logNum", e2 = "logNum"),logNumDivisionOps)
setMethod("^", signature = c(e1 = "logNum", e2 = "numeric"),logNumPowerOps)
setMethod("show", signature = "logNum",showLogNum)


setAs("logNum","numeric", function(from)from@sign * exp(from@x))

# .logNumVector <- setClass("logNumVector", slots = c(data = "list", length = "integer"))
# 
# 
# logNumVector<-function(x){
#   data <- lapply(x, function(x) logNum(x))
#   .logNumVector(data = data, length = length(x))
# }
# 
# logNumVectorFromLogValue<-function(logX){
#   data <- lapply(logX, function(x) logNumFromLogValue(x))
#   .logNumVector(data = data, length = length(logX))
# }
# 
# 
# 
# showLogNumVector<-function(object){
#   result <- sapply(object@data, function(x)as(x,"numeric"))
#   cat("Log transformed number vector.\n")
#   print(result)
# }
# 
# logNumberVectorAdd<-function(e1,e2){
#   if(length(e1)!=length(e2)){
#     stop("Unequal vector length")
#   }
#   result <- sapply(1:length(e1), function(i) e1@data[[i]]+e2@data[[i]])
#   e1@data <- result
#   e1
# }
# logNumberVectorMinus<-function(e1,e2){
#   if(missing(e2)){
#     e1@data=lapply(e1@data, function(x) {
#       x@sign=-x@sign
#       x
#       }
#     )
#     return(e1)
#   }
#   if(length(e1)!=length(e2)){
#     stop("Unequal vector length")
#   }
#   result <- sapply(1:length(e1), function(i) e1@data[[i]]-e2@data[[i]])
#   e1@data <- result
#   e1
# }
# 
# logNumberVectorMultiple<-function(e1,e2){
#   if(length(e1)==1&&length(e1)!=length(e2)){
#     e1@data <- sapply(seq_along(e2),function(i) e1@data)
#     e1@length <- length(e2)
#   }
#   if(length(e2)==1&&length(e1)!=length(e2)){
#     e2@data <- sapply(seq_along(e1),function(i) e2@data)
#     e2@length <- length(e1)
#   }
#   
#   if(length(e1)!=length(e2)){
#     stop("Unequal vector length")
#   }
#   
#   result <- sapply(1:length(e1), function(i) e1@data[[i]]*e2@data[[i]])
#   e1@data <- result
#   e1
# }
# logNumberVectorDivision<-function(e1,e2){
#   if(length(e1)!=length(e2)){
#     stop("Unequal vector length")
#   }
#   result <- sapply(1:length(e1), function(i) e1@data[[i]]/e2@data[[i]])
#   e1@data <- result
#   e1
# }
# logNumberVectorPower<-function(e1,e2){
#   if(length(e1)==1&&length(e1)!=length(e2)){
#     e1@data <- sapply(seq_along(e2),function(i) e1@data)
#     e1@length <- length(e2)
#   }
#   if(length(e1)!=length(e2)){
#     stop("Unequal vector length")
#   }
#   
#   result <- sapply(1:length(e1), function(i) e1@data[[i]]^e2[i])
#   e1@data <- result
#   e1
# }
# 
logNumberVectorSubset <-function(x,i){
  x@length <- length(i)
  x@data<- x@data[i]
  x
}
logNumberVectorSubsetAssign <-function(x,i,value){
  x@data[i] <- value@data
  x
}
# 
# setMethod("show", signature = "logNumVector",showLogNumVector)
# setMethod("+", signature = c(e1 = "logNumVector", e2 = "logNumVector"),logNumberVectorAdd)
# setMethod("-", signature = c(e1 = "logNumVector", e2 = "logNumVector"),logNumberVectorMinus)
# setMethod("-", signature = c(e1 = "logNumVector", e2 = "missing"),logNumberVectorMinus)
# setMethod("*", signature = c(e1 = "logNumVector", e2 = "logNumVector"),logNumberVectorMultiple)
# setMethod("/", signature = c(e1 = "logNumVector", e2 = "logNumVector"),logNumberVectorDivision)
# setMethod("^", signature = c(e1 = "logNumVector", e2 = "numeric"),logNumberVectorPower)
# setMethod("[", signature = "logNumVector",logNumberVectorSubset)
# setMethod("[<-", signature = "logNumVector",logNumberVectorSubsetAssign)
# setMethod("length", signature = "logNumVector",function(x) x@length)
# 
# setAs("logNumVector","numeric", function(from)sapply(1:length(from),function(i)as(from@data[[i]],"numeric")))


.logNumVector <- setClass("logNumVector", slots = c(data = "list", length = "integer"), contains = "Vector")
setMethod("length", signature = "logNumVector",function(x) x@length)
setMethod("[", signature = "logNumVector",logNumberVectorSubset)
setMethod("[<-", signature = "logNumVector",logNumberVectorSubsetAssign)

logNumVector<-function(x){
  data <- as.list(x)
  .logNumVector(data = data, length = length(x))
}








