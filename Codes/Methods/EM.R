
# BASIC EM ALGORITHM
# Example 4.2 on page 91 in Givens and Hoeting
#
# Key idea:  use the EM algorithm to estimate the parameters in the
# peppered moth example described on page 91.

#The observed data (phenotype counts for each variety)
nc=85  #Number of observed carbonaria phenotype
ni=196 #number of observed insularia phenotype
nt=341 #Number of observed typica phenotype
n=nc+ni+nt #total sample size

#Initialize the parameter estimate vectors
#One vector for each phenotype
pct=rep(0,20)
pit=rep(0,20)
ptt=rep(0,20)
pct[1]=1/3  #Initial estimates of the parameters
pit[1]=1/3
ptt[1]=1/3

#This for loop carries out EM for the peppered moth example.
for (t in 2:20) {
  #Equations (4.5)-(4.9), expectation step
  denom1=(pct[t-1]^2+2*pct[t-1]*pit[t-1]+2*pct[t-1]*ptt[t-1])
  ncct=nc*pct[t-1]^2/denom1
  ncit=2*nc*pct[t-1]*pit[t-1]/denom1
  nctt=nc-ncct-ncit
  niit=ni*pit[t-1]^2/(pit[t-1]^2+2*pit[t-1]*ptt[t-1])
  nitt=ni-niit
  nttt=nt
  
  #Equations (4.13-4.15), maximization step
  pct[t]=(2*ncct+ncit+nctt)/(2*n)
  pit[t]=(2*niit+ncit+nitt)/(2*n)
  ptt[t]=(2*nttt+nctt+nitt)/(2*n)
}

#Examine the output
cbind(pct,pit,ptt)

#equation (4.16) on page 94
rcc=sqrt( (diff(pct)^2+diff(pit)^2)/(pct[-20]^2+pit[-20]^2) )
rcc=c(0,rcc)  #adjusts the length to make the table below

#convergence diagnostics defined below (4.16) on page 94
d1=(pct[-1]-pct[20])/(pct[-20]-pct[20])
d1=c(d1,0)
d2=(pit[-1]-pit[20])/(pit[-20]-pit[20])
d2=c(d2,0)

#Table like Table 4.1 on page 94
print(cbind(pct,pit,rcc,d1,d2)[1:9,],digits=5)



####################################################3
## Em algorithm for two-component Poisson mixture
## the input if x, a numeric vector that contains observations from a
## two-component Poisson mixtures.

## generate two-compoment Poisson mixture observations
mu1=2.5
mu2=5.2
p1=0.4
set.seed(1)
group=rbinom(1000,1,p1)
x=rpois(1000,mu1)*group+rpois(1000,mu2)*(1-group)

EM2Pois(x)

EM2Pois=function(x){
  n=length(x)
  z=matrix(0,n,2)   ##define an indicator matrix
  mu1=0            #mu1,mu2,loglike,pi are used to store the EM results that achieves the maximum likleihood from k=0 to current iteration
  mu2=0
  loglike=0
  pi=0

  ##choose random initial values of the Poisson mixture parameters.
  ## repeat for 10 times to search for global maximum.
  for (k in 1:10){
    cat("k=    ",k,"\n")
    p=runif(1)
    u1=runif(1)*max(x)
    u2=runif(1)*max(x)
    current=sum(log(p*dpois(x,u1)+(1-p)*dpois(x,u2)))
    last=current-1

    while(current-last>10^(-8)){
      last=current
      ## E step
      z[,1]=p*dpois(x,u1)/(p*dpois(x,u1)+(1-p)*dpois(x,u2))
      z[,2]=1-z[,1]

      ## M step
      u1=sum(x*z[,1])/sum(z[,1])
      u2=sum(x*z[,2])/sum(z[,2])
      p=sum(z[,1])/n
      current=sum(log(p*dpois(x,u1)+(1-p)*dpois(x,u2)))
      cat("current=   ",current, u1,"\t",u2,"\t",p,"\n")
    }
    if(k==1){
      like=current
      mu1=u1
      mu2=u2
      pi=p  
    } else if(like<current){
      like=current
      mu1=u1
      mu2=u2
      pi=p
      cat(mu1,"    ", mu2,"    ",pi,"\n")
    }
  }
  list(mu1=mu1,mu2=mu2,pi=pi,like=like)
}
