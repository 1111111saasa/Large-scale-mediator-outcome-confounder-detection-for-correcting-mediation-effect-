rm(list=ls())
gc()
library(doParallel)
library(foreach)
library(lars)
library(MASS)
#library(bindata)




##########################################
######## Settings: Parallel ##############
#######################################
ns =100        ###simulation times
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

##########################################
######## Settings: Parameters ############
##########################################
n =1600  ##sample size
B=30    ##refinement times
p =800  ##dimension
ratio=0.5 ##splitting ratio
pi10 = 0.03  ##portion of alpha!=0,beta=0
pi01 = 0.03
pi00 = 1-pi10-pi01*2
pi = c(pi00,pi10,pi01,1-sum(pi00,pi10,pi01))
nonzero = (sum(pi[1:3])*p+1):p
nonzero_beta=(sum(pi[1:2])*p+1):p
q = 0.2 ##FDR level
tau = 3 ##tau=c(0.1,0.3,0.7)signal strength
corr=0.5


Process <- function(n,p,pi,nonzero,nonzero_beta,q,tau,B,correlation=corr){
 
  
#---------------------------data generating---------------------------------- 
  
  ######################################################
  #######generate the AR matrix of covariance###########
  ######################################################
  AR_matrix<- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    rho^exponent
  }
  
  
  #######2022609 Model settings similar to Dai et al (2020,JASA)########
  #####Model: Mk=ak*X+ek;Y=r*X+sum_{k=1}^p bk*Mk+eps,k=1,...,p########
  #####Inputs#####
  #pi: the proportion for each hypothesis; (pi00,pi10,pi01,pi11)
  ModelMultiReg <- function(n,p,pi,cor,tau){
    
    ####generate X####
    X = rnorm(n) #X~N(0,1)
    C =rnorm(n,sd=0.5)
    #X = rbinom(n,1,0.2) #X~Ber(0.2)
    theta_1=c(rep(1,p/2),rep(-1,p/2))*0.2
    theta_0=-0.5
    
    
    ####generate alpha and beta######
    m = pi*p
    coeff = data.frame(alpha=rep(0,p),beta=rep(0,p),C_coeff=rep(0,p))
    
    aInd = (m[1]+1):(m[1]+m[2])
    if(m[4]!=0){aInd = c(aInd,(sum(m[1:3])+1):p)}
    coeff$alpha[aInd] =c(runif(n=length(aInd)/2,min=0.6/3,max=1.2/3),  runif(n=length(aInd)/2,min=-1.2/3,max=-0.6/3))*6
    
    bInd = (sum(m[1:2])+1):p
    
    coeff$beta[bInd] = c(rep(-1,length(bInd)/2),rep(1,length(bInd)/2))*1.2/3*tau#*runif(length(bInd),0.5,2.5) #beta_j ~ 0.3Unif(0.5,2.5)
    coeff$C_coeff=rnorm(p,mean=0.2,sd=1)
    ####generate M: a matrix--n*p#####
    #Em =mvrnorm(n,rep(0,times=p),matrix) #Em ~ N(0,1)
    # Em = matrix(rt(n*p, df = 3), n, p)*sqrt(1/3)
    # Em = matrix(rchisq(n*p, df=1)-1,n,p)*sqrt(2)
    # Em = matrix(rexp(n*p)-1,n,p) # Em ~ exp(1)
    prob_M=X%*%t(coeff$alpha)+theta_1+C%*%t(coeff$C_coeff)
    M=matrix(data=NA,nrow =n ,ncol=p)
    M[,1]=exp(prob_M[,1])/(1+exp(prob_M[,1]))
    M[,1]=sapply(1:n,function(x){rbinom(n=1,size=1,prob=M[x,1])}) 
    #M=sapply(1:n,function(x){rmvbin(n=1,margprob =exp(prob_M[x,])/(1+exp(prob_M[x,])) , bincorr=matrix)})
    for (i in 2:p) {
      prob_M_temp=prob_M[,i]+cor*M[,i-1]
      prob_M_temp=exp(prob_M_temp)/(1+exp(prob_M_temp))
      M[,i]=sapply(1:n,function(x){rbinom(n=1,size=1,prob=prob_M_temp[x])}) 
    }
    #M=lapply(1:p,function(x){rbinom(n=n,size=1,prob=exp(prob_M[,x])/(1+exp(prob_M[,x])))  })
    #M=sapply(1:p, function(x){M[[x]]})
    
    ####generate Y: a matrix--n*p#####
    betaX = 1
    #beta_j=0 for M_j
    #EpY0 = rnorm(n) #Var(Y|x)=1
    prob=betaX*X+theta_0+M%*%coeff$beta
    prob=exp(prob)/(1+exp(prob))
    # Y = betaX*X + M%*%coeff$beta + EpY0
    Y=sapply(1:n,function(x){rbinom(n=1,size=1,prob=prob[x])})    
    return(list(X=X,Y=Y,M=M))
  }
  

#---------------------------DIMOC-----------------------------------------------
  
  ######################################################
  ## Estimate alpha & beta for Y=aX+sum_{k=1}^p M_k+e ##
  ######################################################

  EstimateCoef2 <- function(X,M,Y,p=ncol(M),nonbeta){
    
    EsAlpha= sapply(1:p,function(t){
      temp = summary(glm(M[,t]~X,family=binomial(link="logit"),control=list(maxit=100)))
      #temp$coefficients[,2] = temp$coefficients[,2]/temp$sigma #no need to estimate sigma
      return(temp$coefficients[2,])
    })
    
    EsAlpha = data.frame(t(EsAlpha)) #4 columns: 1-estimate coefficient; 2- std. error; 3-t value; 4-p value
    
    ###Initialize EsBeta###
    EsBeta = EsAlpha
    EsBeta$Estimate = 0;
    EsBeta$Std..Error = 1;
    EsBeta$z.value = 0;
    EsBeta$Pr...z.. = 1; ###
    
    lmYXM = summary(glm(Y~X+M[,nonbeta],family=binomial(link="logit"),control=list(maxit=100)))
    #lmYXM$coefficients[,2] = lmYXM$coefficients[,2]/lmYXM$coefficients[,2]
    EsBeta[nonbeta, ] = lmYXM$coefficients[-c(1,2),]
    
    return(list(EsAlpha=EsAlpha,EsBeta=EsBeta))  
  }


  ###########################################################
  ########     Function for linear regression      ##########
  ###########################################################  
  regression=function(y,x){
    glm(as.factor(y)~x,family=binomial(link="logit"),control=list(maxit=100))
  }
  
  ############################################################
  ############   Calculate the empirical FDR, power ##########
  ############################################################
  fdr_ap_fun <- function(S01,S02,W_rej){
    S_free <- 2*((S01<S02)+0)
    S_null <- (S01>S02)+0+S_free        ##here S_null=0 for alpha=beta;S_null=1 for alpha=0,beta!=0;S_null=2 for alpha!=0,beta=0 
    false_rej <- which((S_null==0)&(W_rej==1))
    true_rej <- which((S_null==1)&(W_rej==1))
    power <- length(true_rej)/max(sum(S_null==1),1)
    FDR <- length(false_rej)/max(sum(W_rej),1)
    fdr_ap_result <- list(FDR,power)
    return(fdr_ap_result)
  }
  
  
  
  ############################################################
  ####### Given the L2 threshold, it returns L1 threshold#####
  ############################################################  
  W1_thre_new_fun <- function(W1,W2,L2,alpha,options,p){
    W1_abs_sort = sort(abs(W1))
    W1_abs_sort <- W1_abs_sort[which(W1_abs_sort!=0)]
    W_frame <- data.frame(ind=c(1:p),W1=W1,W2=W2)
    
    if(options=='+'){
      Ta = sapply(W1_abs_sort,function(x){
        (1+dim(W_frame[(W_frame$W1<=(-x))&(W_frame$W2<=L2),])[1] +  dim(W_frame[(W_frame$W1>=(x))&(W_frame$W2<=L2),])[1] -
           dim(W_frame[((W_frame$W2>=(-L2))&(W_frame$W2<0)&(W_frame$W1>=x)),])[1] - dim(W_frame[(W_frame$W2<=0)&(W_frame$W1>=x),])[1] )/
          max((dim(W_frame[((W_frame$W2<=(L2))&(W_frame$W1>=x)),])[1]),1)
      })
    }else{
      Ta = sapply(W1_abs_sort,function(x){
        (dim(W_frame[(W_frame$W1<=(-x))&(W_frame$W2<=L2),])[1] +  dim(W_frame[(W_frame$W1>=(x))&(W_frame$W2<=L2),])[1] -
           dim(W_frame[((W_frame$W2>=(-L2))&(W_frame$W2<0)&(W_frame$W1>=x)),])[1] - dim(W_frame[(W_frame$W2<=0)&(W_frame$W1>=x),])[1] )/
          max((dim(W_frame[((W_frame$W2<=(L2))&(W_frame$W1>=x)),])[1]),1)
      })
    }
    if(min(Ta)>alpha){
      W1_thre <- 1000000
      estimate_DTFDR=100
    }else{
      W1_thre = min(W1_abs_sort[which(Ta<=alpha)])
      ind=which(W1_abs_sort==min(W1_abs_sort[which(Ta<=alpha)]))
      estimate_DTFDR=Ta[[ind]]
    }
    
    
    det=W_frame[(W_frame$W1>=W1_thre) & (W_frame$W2<=L2) ,]$ind
    rej = rep(0,p)
    rej[det] = 1
    W_rej_list <- list(W1_thre,L2,rej) #return L1,L2,reject set
    return(list(W_rej_list=W_rej_list,estimate_DTFDR=estimate_DTFDR))
  }
  
  
  ############################################################
  # It inputs each value of W2 as the L2 threshold into the ##
  # above function and returns the best threshold pairs.   ###
  ############################################################
  opt_L1_L2_fun <- function(W1,W2,alpha,options,p){
    L2_can <- sort(abs(W2))
    L2_can <- L2_can[which(L2_can!=0)]
    num_rej <- sapply(L2_can,function(x){
      sum(W1_thre_new_fun(W1,W2,x,alpha,options,p)$W_rej_list[[3]]) 
    })
    
    L2 <- max(L2_can[which(num_rej==max(num_rej))])
    L1_L2_rej_list <- W1_thre_new_fun(W1,W2,L2,alpha,options,p)  
    return(L1_L2_rej_list)
  }
  
  #####################################################
  ###### main function for detecting confounders  #####
  #####################################################   
  DIMOC=function(X,M,Y,n,p,q,pi,B,ratio=0.5,nonzero_beta){
    result_rej2<-rep(0,p)
    All_result <- vector("list",length = B)
    for (b in 1:B)
    {
      #######DIMOC#############
      data=list("Y"=Y,'covariate'=cbind(X,M))
      S01=rep(0,time=p);S02=rep(0,time=p)
      support_S02=c((pi[1]*p+1):(sum(pi[1:2])*p),(sum(pi[1:3])*p+1):p)      ##beta
      S02[support_S02]=1
      support_S01=(sum(pi[1:2])*p+1):p       ##alpha
      S01[support_S01]=1
      
      ###Construct filter for beta
      D_1= sample(1:n,floor(n/2))
      D_2= (1:n)[-D_1]  
      support_1=c(1,nonzero_beta+1)
      
      Y1=data$Y[D_2]
      x1=data$covariate[D_2,1]
      M1=data$covariate[D_2,setdiff(support_1,1)]
      reg_coefficients1=glm(as.factor(Y1)~x1+M1,family=binomial(link="logit"),control=list(maxit=100))$coefficients[-c(1,2)]
      #reg_coefficients1=glm(as.factor(Y1)~x1+M1,family=binomial(link="logit"),control=list(maxit=100))
      
      W_11=rep(0,time=p+1)
      W_11[setdiff(support_1,1)]=reg_coefficients1
      W_11=W_11[-1]
      
      
      Y2=data$Y[D_1]
      x2=data$covariate[D_1,1]
      M2=data$covariate[D_1,setdiff(support_1,1)]
      reg_coefficients2=glm(as.factor(Y2)~x2+M2,family=binomial(link="logit"),control=list(maxit=100))$coefficients[-c(1,2)]
      W_12=rep(0,time=p+1)
      W_12[setdiff(support_1,1)]=reg_coefficients2
      W_12=W_12[-1]
      
      W_1=W_11*W_12
      
      ###alpha
      regression_list21=apply(data$covariate[D_1,-1],2,regression,x=data$covariate[D_1,1])                 #fit the data with  model(1)
      regression_list22=apply(data$covariate[D_2,-1],2,regression,x=data$covariate[D_2,1])
      coefficients21=lapply(regression_list21,function(x) x$coefficients[2])
      coefficients22=lapply(regression_list22,function(x) x$coefficients[2])
      W_2=unlist(coefficients21)*unlist(coefficients22)
      
      #choose L1,L2
      L1_L2_rej <- opt_L1_L2_fun(W_1,W_2,q,"",p)$W_rej_list
      estimate_DTFDR=opt_L1_L2_fun(W_1,W_2,q,"",p)$estimate_DTFDR
      L1 <- L1_L2_rej[[1]]
      L2 <- L1_L2_rej[[2]]
      W_rej <- L1_L2_rej[[3]]
      
      
      emp_FDR_power <- fdr_ap_fun(S01,S02,W_rej)
      emp_FDR <- emp_FDR_power[[1]]
      emp_power <- emp_FDR_power[[2]]
      
      result_DT<- list("W1_DT"=W_1,"W2_DT"=W_2,"L1_DT"=L1,"L2_DT"=L2,"reject_indexDT"=W_rej,
                       "FDR_DT"=emp_FDR,"Power_DT"=emp_power,"Support1_DT"=support_1
      )
      result_rej_b2<-result_DT[[5]]
      result_rej2<-cbind(result_rej2,result_rej_b2)
      All_result[[b]] <- result_DT
    }
    result_rej2<-result_rej2[,-1]
    sum_rej2<-rowSums(result_rej2)
    bag_rej2<-(sum_rej2>=ceiling(B/2))+0
    
    num<-sapply(1:B, function(x){
      sum((bag_rej2+result_rej2[,x]==2)+0) + sum((bag_rej2+result_rej2[,x]==0)+0)
    })
    final_decide<-which(num==max(num)) 
    final_decide<-final_decide[1]
    
    result_list <- All_result[[final_decide]] #list(emp_FDR,emp_power,result_rej)
    return(result_list) 
  }
  
  
  
  Dat = ModelMultiReg(n,p,pi,corr,tau)     ###generate the data
  
  Res=DIMOC(Dat$X,Dat$M,Dat$Y,n,p,q,pi,B,ratio=0.5,nonzero_beta)
  
 
  return(list(Dat = Dat,Res = Res))
  closeAllConnections()
}

##########################################
######## Simulation: Process #############
##########################################

set.seed(1000)

nonzero=(sum(pi[1:3])*p+1):p

print("running...")
start_time=proc.time()

RES0 = foreach(i=1:ns,.packages = c("MASS","lars")) %dopar% Process(n,p,pi,nonzero,nonzero_beta,q,tau,B)    
cat("模拟 time:", (proc.time()-start_time)[3][[1]],"秒","\n")
saveRDS(RES0, file=paste0("./Simulation/pro=_",pi01,"_corr=_",corr,"_n=_",n,"_p=_",p,"_t=_",tau,"_ratio=_",ratio,"_refine=_",B,"_mediation.RData"))

stopCluster(cl)

print("Down")

#saveRDS(RES0,"n=600,p=1000,t=3.5.RData")

##########################################
###### Summary of Simulation Results #####
##########################################
# fdp1=lapply(RES0,function(x){x$FDP1})
# mean(unlist(fdp1))
