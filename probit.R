# bayesian-probit-model

#load()
data <- read.table("C:/GAUSS/sample/data/uregdata.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
# 1:height 2:intercept 3:weight 4:foot size 5:dammy(man=1,woman=0)
y=as.matrix(data[,1])
x=as.matrix(data[,2:ncol(data)])
y <- ifelse(y>mean(y),y,mean(y)) 
n=nrow(y)
k=ncol(x)
library(bayesm)
library(msm)
# model:   y = x*beta + e   e ~ N(0,sigma) when y>=0   
#          y = b0 + x1*b1 + x2*b2 + x3*b3 + e  
# truncated interval is 0.
truncated_point <- min(y)
norm_truncated_point <- 0
#  OLS estimate  
  b_ols=chol2inv(chol(crossprod(x)))%*%t(x)%*%y 
  res=y-x%*%b_ols
  s_ols=t(res)%*%res/(n-k)
              
              
#  mcmc iteration  #
  nblow=1000     # burn in#
  smcmc=1000
  
#  prior  
 #  beta ~ N(u0,v0)  
    u0=matrix(0,k,1)      # mean
    v0=diag(k)*1000     # variance
    v0i=solve(v0)
    v0iu0=v0i%*%u0
  
 #  sigma ~ IG(s0,r0)  
    s0=2;       # scale
    r0=2;       # shape
    rn=r0+n;	  # posterior shape parameter
    
# initial value and box for gibbs sampler  ***/
    beta=matrix(0,k,1)            # initial value of beta
    betag=matrix(0,smcmc,k)       # keeping mcmc sample
    
    sigma=1                     # initial value of sigma
    sigmai=1/sigma              # inverse of sigma
    sigmag=matrix(0,smcmc,1)      # keeping mcmc sample    
    
    result<- list(3)
    
  uprobit=function(y,x,beta,sigmai,v0i,u0,s0,rn,truncated_point,norm_truncated_point){
     
      vec<- (y==truncated_point)
      
      mu <- x%*%beta
     
      y0<-y
      y0[y==truncated_point]<-rtnorm(sum(vec),x[vec,]%*%beta, sd=5.2,lower = -Inf,upper = truncated_point) 
      #  beta : p(beta|sigma,y,x) ~ N(bm,bv)  
      bv=chol2inv(chol(v0i+sigmai*(t(x)%*%x) ))		# conditional posterior variance
      bm=bv%*%(v0i%*%u0+sigmai*t(x)%*%y0)		# conditional posterior mean
      beta=t(chol(bv))%*%rnorm(k,0,1)+bm		# generate mcmc sample from Normal	
      #beta <- as.matrix(beta)
      #print(beta)
      #  sigma : p(sigma|beta,y,x) ~ IG(sn,rn)  
      res=y0-x%*%beta					# residual vector
      sn=s0+t(res)%*%res					# conditional posterior scale parameter
      sigma <- (sn/(rgamma(1,shape=(rn/2),scale = 2)))	# generate mcmc sample from Inverted Gamma	
      sigma <- as.numeric(sigma)
      sigmai <- 1/sigma
      
      return(list(beta,sigma,sigmai))
    }
    
# do mcmc  
 
    for (i1 in c(1:nblow)) {
        imcmc=as.numeric(i1)
        result=uprobit(y,x,beta,sigmai,v0i,u0,s0,rn,truncated_point,norm_truncated_point)
        beta=result[[1]]
        sigma=result[[2]]
        sigmai=result[[3]]
      }    
    
    
    for (i1 in c(1:smcmc) ) {
      imcmc=as.numeric(i1)
      
      result=uprobit(y,x,beta,sigmai,v0i,u0,s0,rn,truncated_point,norm_truncated_point)
      beta=result[[1]]
      sigma=result[[2]]
      sigmai=result[[3]]
      
      betag[imcmc,]=as.matrix(t(beta))     
      
      sigmag[imcmc,]=as.matrix(sigma)

    }
    
# posterior  
    library(matrixStats)      # Matrice
    betam=colMeans(betag)     # posterior mean
    betas=colSds(betag)         # posterior sd
    
    sigmam=colMeans(sigmag)
    sigmas=colSds(sigmag)    
    
# convergence check & parameter distribution   
    par(mfrow=c(2,2))
    plot(betag[,1],type = "l")
    hist(betag[,1])
    plot(betag[,2],type = "l")
    hist(betag[,2])


   
    
