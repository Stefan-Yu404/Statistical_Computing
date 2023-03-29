##jackknife
set.seed(0)
x=rnorm(100,0,1)
s2=var(x) ##unbiased
s2_b=var(x)*(100-1)/100
s2
s2_b

s2_j=numeric(100)
for (i in 1:100){
  x_temp=x[-i]
  s2_j[i]=var(x)*(99-1)/99
}
jackknife=100*s2_b-99/100*sum(s2_j)


##permutation test
set.seed(0)
x=rnorm(100,0,1)
y=rnorm(80,1,2)

diff0=mean(x)-mean(y)

diff_perm=numeric(10000)
for (i in 1:10000){
  xy=sample(c(x,y),length(x)+length(y),replace=FALSE)
  diff_perm[i]=mean(xy[1:10])-mean(xy[11:18])
}
hist(diff_perm,breaks=20)
abline(v=diff0)
pvalue=sum(diff_perm<diff0)/10000
pvalue