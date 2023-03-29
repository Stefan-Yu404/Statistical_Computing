## R codes for chapter 9

###########################################
## bootstrap, Jackknife, and permutation
###########################################

setwd("~/stat344F22/")

x=c(0.44,  0.06, -0.54,  0.40, -0.06,  0.42, -1.52,  1.09, -0.79,  0.32)
xbar=mean(x)

##nonparatric bootstrap for variance of sample mean E.g.2 in boostrap lectures
B=10000
average=numeric(B)
for (i in 1:B){
  y=sample(x,10,replace=TRUE)
  average[i]=mean(y)
}

var(average) #bootstrap estimate of v(x_bar)
var(x)/10 # estimate calcualted from s^2/n


## bootstrapping regression for Copper-nickle alloy data in table 9.2 on page 257
## Example 9.3, 9.4, 9.5, 9.6, 9.7

dat <- read.table("datasets/alloy.dat", header=T) #available from canvas under boostrap
data.size <- dim(dat)[1]
x=dat[,1]
y=dat[,2]


## Example 9.3
# use linear regression to fit the model Y = beta0 + beta1*x + epsilon
temp0 <- lm(y~x)             # linear regression of y on x
summary(temp0)
beta0 <- temp0$coef[1]       # estimate for beta0
beta1 <- temp0$coef[2]       # estimate for beta1
theta <- beta1/beta0        # estimate of beta1/beta1
theta
#-0.1850722 

vb <- vcov(temp0)          # variance-covariance matrix of (beta0, beta1)




# use (9.16) on page 263 to calculate V(F^hat), estimate of var(beta1/beta0)
vfhat <- (beta1/beta0)^2*(vb[2,2]/beta1^2+vb[1,1]/beta0^2-2*vb[1,2]/(beta0*beta1))
vfhat # 6.981006e-05 

##bootstrap for estimation of v(beta1/beta0)
B <- 10000
n <- length(x)
yfitted=temp$fitted          # size of original sample
sigma=sqrt(sum(temp$residual^2)/temp$df.residual)  #temp$residual is the residual, and temp$df.redidual is degree of freedom
##boostrap estimate for v(\beta_1)

beta1_B=numeric(B)

for(i in 1:B) {             # loop begins
  tempRes <- rnorm(n, 0,sigma);  # bootstrapped residuals
  tempy <- yfitted + tempRes;  # boostrapped y
  temp <- lm(tempy~x);
  beta1_B[i] <- temp$coef[2]     # estimate for beta1
}  
  var(beta1_B)


###boostrap for variabce of beta_1/beta_0

thetastar0 <- numeric(B)    # bootstrapped estimates for theta1/theta0
Rfstar0 <- numeric(B)        # bootstrapped (theta^hat-theta)/sqrt(V(F^hat))

set.seed(2)
for(i in 1:B) {             # loop begins
  tempRes <- rnorm(n, 0,sigma);  # bootstrapped residuals
  tempy <- yfitted + tempRes;  # boostrapped y
  temp <- lm(tempy~x);
  tempb0 <- temp$coef[1]     # estimate for beta0
  tempb1 <- temp$coef[2]     # estimate for beta1
  thetastar0[i] <- tempb1/tempb0;
  tempV <- vcov(temp)        # variance-covariance matrix of (beta0, beta1), for Example 9.7
  tempV2 <- (tempb1/tempb0)^2*(tempV[2,2]/tempb1^2+tempV[1,1]/tempb0^2-2*tempV[1,2]/(tempb0*tempb1));  #variance for the b1/b0
  Rfstar0[i] <- (tempb1/tempb0-theta)/sqrt(tempV2);
}
summary(thetastar0)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.2188 -0.1905 -0.1850 -0.1849 -0.1792 -0.1538 
var(thetastar0)
hist(thetastar0, xlab="Boostrapped beta1/beta0 based on residuals", nclass=40)


## bootstrapping the residuals,nonparametric
set.seed(2)
B <- 10000
n <- length(x)
yfitted=temp$fitted          # size of original sample
residual <- y - yfitted     # temp$fitted is the fitted values 
thetastar0 <- numeric(B)    # bootstrapped estimates for theta1/theta0
Rfstar0 <- numeric(B)        # bootstrapped (theta^hat-theta)/sqrt(V(F^hat))
for(i in 1:B) {             # loop begins
  tempRes <- sample(residual, size=n, replace=T);  # bootstrapped residuals
  tempy <- yfitted + tempRes;  # boostrapped y
  temp <- lm(tempy~x);
  tempb0 <- temp$coef[1]     # estimate for beta0
  tempb1 <- temp$coef[2]     # estimate for beta1
  thetastar0[i] <- tempb1/tempb0;
  tempV <- vcov(temp)        # variance-covariance matrix of (beta0, beta1), for Example 9.7
  tempV2 <- (tempb1/tempb0)^2*(tempV[2,2]/tempb1^2+tempV[1,1]/tempb0^2-2*tempV[1,2]/(tempb0*tempb1));  #variance for the b1/b0
  Rfstar0[i] <- (tempb1/tempb0-theta)/sqrt(tempV2);
}
summary(thetastar0)
var(thetastar0)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.2111 -0.1903 -0.1851 -0.1849 -0.1797 -0.1576 
hist(thetastar0, xlab="Boostrapped beta1/beta0 based on residuals", nclass=40)


# bootstrapping the cases
set.seed(2)
B <- 10000
n <- length(x)              # size of original sample
thetastar <- rep(0, B)      # bootstrapped estimates for theta=theta1/theta0
Rfstar <- rep(0, B)         # bootstrapped (theta^hat-theta)/sqrt(V(F^hat))
for(i in 1:B) {             # loop begins
  xyindex <- sample(seq(1:n), size=n, replace=T);
  tempx <- x[xyindex];
  tempy <- y[xyindex];
  temp <- lm(tempy~tempx);
  tempb0 <- temp$coef[1]     # estimate for beta0
  tempb1 <- temp$coef[2]     # estimate for beta1
  thetastar[i] <- tempb1/tempb0;  # estimate for theta
  tempV <- vcov(temp)        # variance-covariance matrix of (beta0, beta1), for Example 9.7
  tempV2 <- (tempb1/tempb0)^2*(tempV[2,2]/tempb1^2+tempV[1,1]/tempb0^2-2*tempV[1,2]/(tempb0*tempb1));  #variance for the b1/b0
  Rfstar[i] <- (tempb1/tempb0-theta)/sqrt(tempV2);
}
summary(thetastar)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.2406 -0.1907 -0.1851 -0.1863 -0.1806 -0.1515 
hist(thetastar, xlab="Boostrap estimates of theta=beta1/beta0", nclass=40)

## Example 9.4
bias <- mean(thetastar) - theta      # mean value of thetahat^* - thetahat
bias   #-0.001192565
theta - bias                         # bias-corrected bootstrap estimate
#-0.1838797 

## Example 9.5
quantile(thetastar, c(0.025, 0.975)) # a 95% bootstrapped confidence interval for beta1/beta0
#      2.5%      97.5% 
#-0.2052362 -0.1732016  

## Example 9.6
thetai <- rep(0, n)         # thetahat_(-i), thetahat computed omitting the ith observation
for(i in 1:n) {
  tempx <- x[-i];
  tempy <- y[-i];
  temp <- lm(tempy~tempx)$coef;
  thetai[i] <- temp[2]/temp[1];
}
psii <- mean(thetai) - thetai                # psi_i, (9.13) on page 262
a <- sum(psii^3)/(6*(sum(psii^2))^(3/2))     # a, (9.12) on page 262
a # 0.04859325
b <- qnorm( sum(thetastar <= theta)/B )      # b, line right above (9.12) on page 262
b # 0.007018617
alpha <- 0.05               # for 95% confidence interval
beta1 <- pnorm(b+(b+qnorm(alpha/2))/(1-a*(b+qnorm(alpha/2))))     # beta_1, (9.10) on page 262
beta1  # 0.03781242
beta2 <- pnorm(b+(b+qnorm(1-alpha/2))/(1-a*(b+qnorm(1-alpha/2)))) # beta_1, (9.11) on page 262
beta2  # 0.9854408
quantile(thetastar, c(beta1, beta2)) # a BCa 95% bootstrapped confidence interval
# 3.781242%  98.54408% 
#-0.2029387 -0.1720287 

## Example 9.7, confidence interval by bootstrap t, using cases
hist(Rfstar, xlab="R(X^*, F^hat)", nclass=80) # histogram of R(X^*, F^hat)
xi <- quantile(Rfstar, c(0.025, 0.975)) # quantiles of Ghat^*
#     2.5%     97.5% 
#-2.183922  2.086920 
c(theta-sqrt(vfhat)*xi[2], theta-sqrt(vfhat)*xi[1])
# -0.2025089 -0.1668250 

#confidence interval by bootstrap t,  based on boostrapping residuals
hist(Rfstar0, xlab="R(X^*, F^hat)", nclass=80) # histogram of R(X^*, F^hat)
xi <- quantile(Rfstar0, c(0.025, 0.975)) # quantiles of Ghat^*
#     2.5%     97.5% 
#-2.233093  2.225583
c(theta-sqrt(vfhat)*xi[2], theta-sqrt(vfhat)*xi[1])
# -0.2036675 -0.1664142 


