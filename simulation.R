######################################################################################
#### Creating function to generate time varying spatio-temporal covariance function ######
######################################################################################
library(MASS)
library(mvtnorm)
library(fields)
my.matern<-function(h,c,sigma,nu)
{
  h[h==0]<-1e-10
  num1<-(sigma^2)*(2^(1-nu))/gamma(nu)
  num2<-(h*c)^nu
  num3<-besselK(x=(h*c), nu=nu)
  return(num1*num2*num3)
}



###### First we create a psi function #######
psi<-function(t,a,alpha,beta)
{
  temp<-(a*(t^alpha)+1)^beta
  return(temp)
}

nT = 50
tseq<- 1:nT#seq(0,10,length.out = nT)

#### range and smoothness for different timepoints
t0 = 25
c_vec = ifelse(tseq>t0, 5,5)
nu_vec = ifelse(tseq>t0, 0.5,0.5)


#### Plotting psi #####
#plot(tseq,1/psi(t=tseq^2,a=1,alpha=1,beta = 0.5),main=bquote(psi~"(t)"),xlab="t",ylab=bquote(psi~"(t)"))
#lines(tseq,1/psi(t=tseq^2,a=1,alpha=1,beta = 0.5))
#lines(tseq,1/psi(t=tseq^2,a=1,alpha=1,beta = 0.6))
#lines(tseq,1/psi(t=tseq^2,a=1.5,alpha=1,beta = 0.5))


##### Creating function for f_cvec  #######
f_cvec<-function(t,c_vec,a,alpha,beta)
{
  c_mean<-mean(c_vec)
  len<-length(c_vec)
  al.creat<-function(a1,a2)
  {
    return(((1/(a1^2))+(1/(a2^2)))/2)
  }
  d2<-outer(c_vec,c_vec,al.creat)
  tmat<-outer(t,t,function(t1,t2) (abs(t1-t2))^2)
  d1<-psi(t=tmat,a=a,alpha = alpha,beta=beta)/(c_mean^2)
  d3<-matrix(psi(t=0,a=a,alpha = alpha,beta=beta)/(c_mean^2),nrow=nrow(d1),ncol=ncol(d1))
  return(1/sqrt(d1+d2-d3))
}



zeta<-function(t,nu_vec,c_vec,a,alpha,beta,d=2,delta)
{
  n2<-outer(nu_vec,nu_vec,function(a,b) gamma((a+b)/2))
  f<-(f_cvec(t=t,c_vec = c_vec,a=a,alpha=alpha,beta=beta))
  n1<-f^d
  dig<-(diag(f))^(d/2)
  #d1<-dig%*%dig
  d1<-outer(dig,dig,function(a,b) a*b)
  d2<-outer(nu_vec,nu_vec,function(a,b) sqrt(gamma(a)*gamma(b)))
  #n2<-outer(c_vec,c_vec,function(a,b) ((1/a)^(d/4))+((1/b)^(d/4)))
  #d2<-(f_c_vec(t=t,c_vec = c_vec,a=a,alpha=alpha,beta=beta))^(d/2)
  tmat<-outer(t,t,function(t1,t2) (abs(t1-t2))^2)
  d3<-psi(t=tmat,a=a,alpha = alpha,beta=delta) # purely temporal corr
  
  return((n1*n2)/(d1*d2*d3))
}



mya <-0.5 #0.05
myalpha<-0.5 #0.5
mybeta<-0.7 #0.7
mydelta<-0

myf_c<-f_cvec(t=tseq,c_vec = c_vec,a=mya,alpha=myalpha, beta=mybeta)

myzeta<-zeta(t=tseq,nu_vec = nu_vec,c_vec = c_vec,a=mya,alpha=myalpha, beta= mybeta, delta=mydelta)

mynus<-outer(nu_vec,nu_vec,function(a,b) (a+b)/2)
#diag(myzeta)
####### Creating irregular spatial grid #########

nx = 5
n = nx^2
grid = expand.grid(1:nx, 1:nx)
xy = matrix(runif(n*2, -0.4, 0.4), n, 2)
loc = (grid-0.5+xy)/nx

##########################################################
##### Creating  spatio-temporal covariance function ######
##########################################################
sl<-length(loc[,1])
tl<-length(tseq)
spcov<-matrix(NA,ncol = sl*tl,nrow=sl*tl)
dist.sp<-rdist(loc)



######### Optimized cov computing ##############
uniq.dist<-sort(unique(c(dist.sp)))

##### Unique indexing #####
mat.index<-NULL
for(i in 1:length(uniq.dist))
{
  mat.index[[i]]<-cbind(which(dist.sp==uniq.dist[i],arr.ind = T),i)
}
#cbind(which(dist.sp==uniq.dist[2],arr.ind = T),2)
mat.index<-do.call(rbind,mat.index)
spcov2<-matrix(NA,ncol = sl*tl,nrow=sl*tl)
temp<-matrix(NA,ncol = sl,nrow = sl)


system.time(for(i in 1:tl)
{
  for(j in i:tl)
  {
    tmp1<-myzeta[i,j]*my.matern(h=uniq.dist,
                                c=myf_c[i,j],
                                nu=mynus[i,j],
                                sigma=1)
    temp[mat.index[,-3]]<-tmp1[mat.index[,3]]
    spcov2[((i-1)*sl+1):((i-1)*sl+sl),((j-1)*sl+1):((j-1)*sl+sl)]<-spcov2[((j-1)*sl+1):((j-1)*sl+sl),((i-1)*sl+1):((i-1)*sl+sl)]<-temp
  }
}
)




## simulating data from nonstationary covariance
#### range and smoothness for different timepoints


st_covmat_func =  function(c1,c2,nu, a, alpha, beta, delta){
  t0 = 25
  tseq = 1:tl
  c_vec = ifelse(tseq>t0, c2,c1)
  nu_vec = ifelse(tseq>t0, nu,nu)
  
  myf_c<-f_cvec(t=tseq,c_vec = c_vec,a,alpha, beta)
  myzeta<-zeta(t=tseq,nu_vec = nu_vec,c_vec = c_vec,a=a,alpha=alpha, beta=beta, delta=delta)
  mynus<-outer(nu_vec,nu_vec,function(a,b) (a+b)/2)
  
  spcov<-matrix(NA,ncol = sl*tl,nrow=sl*tl)
  temp<-matrix(NA,ncol = sl,nrow = sl)
  
  for(i in 1:tl)
  {
    for(j in i:tl)
    {
      tmp1<-myzeta[i,j]*my.matern(h=uniq.dist,
                                  c=myf_c[i,j],
                                  nu=mynus[i,j],
                                  sigma=1)
      temp[mat.index[,-3]]<-tmp1[mat.index[,3]]
      spcov[((i-1)*sl+1):((i-1)*sl+sl),((j-1)*sl+1):((j-1)*sl+sl)]<-spcov[((j-1)*sl+1):((j-1)*sl+sl),((i-1)*sl+1):((i-1)*sl+sl)]<-temp
    }
  }
  return(spcov)
}

spcov2 = st_covmat_func(c1=5,c2=5,nu=0.5, a=mya, alpha=myalpha, beta = mybeta, delta= mydelta)
check<-chol(spcov2)
#rm(check)
library(mvtnorm)
set.seed(124)
nsim<-1
#mysim<-mvrnorm(n=nsim,mu=rep(0,times=sl*tl),Sigma = spcov2)
mysim<-rmvnorm(n=nsim,mean = c(rep(0, times=sl*tl/2), rep(0, times=sl*tl/2)), sigma=spcov2,method="chol")
Z<-c(mysim)


quilt.plot(loc, Z[1:25],  nx=7, ny=7,col = viridis::viridis(100), zlim=c(-3,3))
quilt.plot(loc, Z[626:650],  nx=7, ny=7,col = viridis::viridis(100), zlim=c(-3,3))


mean(Z[601:625])


#### Fitting full likelihood########


logL_st = function(par,Z,X,sl,tl,k, hypo){
  a = mya; alpha = myalpha; beta = mybeta
  sigma= 1;  nu = 0.5;  delta = mydelta
  
  if(hypo==0){
    c = par[1]#;a = par[2]; alpha= par[3]; beta=par[4]
    c1=c2=c
  } else if (hypo==1){
    c1 = par[1]; c2=par[2]#; a = par[3]; alpha= par[4]; beta=par[5]
  }
  
  
  st_cov_mat= matrix(NA,ncol = sl*tl,nrow=sl*tl)
  temp = matrix(NA,ncol = sl,nrow=sl)
  tseq = 1:tl
  
  c_vec = ifelse(tseq>k, c2,c1)
  nu_vec = rep(nu, times = nT)
  
  
  myzeta<-zeta(t=tseq,nu_vec = nu_vec,c_vec = c_vec,a=a,alpha=alpha, beta = beta,delta =delta)
  myf_c<-f_cvec(t=tseq,c_vec = c_vec,a =a,alpha=alpha, beta=beta)
  mynus<-outer(nu_vec,nu_vec,function(a,b) (a+b)/2)
  
  for(i in 1:tl) {
    for(j in i:tl) {
      tmp1<-myzeta[i,j]*my.matern(h=uniq.dist,
                                  c=myf_c[i,j],
                                  nu=mynus[i,j],
                                  sigma=1)
      temp[mat.index[,-3]]<-tmp1[mat.index[,3]]
      st_cov_mat[((i-1)*sl+1):((i-1)*sl+sl),((j-1)*sl+1):((j-1)*sl+sl)]<-st_cov_mat[((j-1)*sl+1):((j-1)*sl+sl),((i-1)*sl+1):((i-1)*sl+sl)]<-temp
    }
  }
  
  chol_st<-chol(st_cov_mat)
  inv_st = chol2inv(chol_st)
  #beta_gls = solve(t(X)%*%inv_st%*%X)%*%t(X)%*%inv_st%*%Z 
  
  #loglik = -0.5*2*sum(log(diag(chol_st))) -0.5*t(Z-X%*%beta_gls)%*%inv_st%*%(Z-X%*%beta_gls)
  loglik = -0.5*2*sum(log(diag(chol_st))) -0.5*t(Z)%*%inv_st%*%(Z)
  return(-loglik)
}

X = rep(1,times =sl*tl)

start.time = Sys.time()
fit.optim = optim(par=c(5,0.6,0.5,0.7),logL_st,Z=Z,X=X,sl=sl,tl=tl, k= tl, hypo=0)
end.time = Sys.time()
end.time-start.time
fit.optim





### log likelihood using markov property


sp_cov_mat = function(t1,t2, sigma, c1, c2, nu,a,alpha,beta, delta, k ){
  
  sp_mat<-matrix(NA,ncol = sl,nrow = sl)
  
  c_vec = ifelse(tseq>k, c2,c1)
  nu_vec = ifelse(tseq>k, nu,nu)
  c_mean = mean(c_vec)
  
  term1 = gamma((nu_vec[t1]+nu_vec[t2])/2)/ sqrt(gamma(nu_vec[t1])*gamma(nu_vec[t2]))
  term2 = 1/sqrt((c_vec[t1]^2)*(c_vec[t2]^2))
  term3 = fc_t1t2 = 1/(psi(t=(t1-t2)^2,a=a,alpha = alpha,beta=beta)/(c_mean^2) +
                         ((1/(c_vec[t1]^2))+(1/(c_vec[t2]^2)))/2 - psi(t=0,a=a,alpha = alpha,beta=beta)/(c_mean^2))
  term4 = 1/(psi(t=(t1-t2)^2,a=a,alpha = alpha,beta=delta)) #pure temporal cor
  
  tmp1<-term1*term2*term3*term4*my.matern(h=uniq.dist,
                                          c=fc_t1t2^(1/2),
                                          nu=(nu_vec[t1]+nu_vec[t2])/2,
                                          sigma=1)
  sp_mat[mat.index[,-3]]<-tmp1[mat.index[,3]]
  return(sp_mat)
}




logL_markov_mean = function(par, Z,sl,tl, k, hypo){
  
  a = mya; alpha = myalpha; beta = mybeta
  sigma= 1; delta = mydelta; c1=c2 =c; nu=0.5
  
  if(hypo==0){
    museg1= museg2 = par[1]; c = par[2]
  } else if (hypo==1){
    museg1 = par[1]; museg2 = par[2];c = par[3]
  }
  mu = ifelse(tseq>k,museg2, museg1)
  #mu=0
  #mu1= mu2 = mu*rep(1 , times=sl)
  
  logL=0
  for (i in tl:1) {
    
    if(i==1){
      mu1  = mu[i]*rep(1 , times=sl)
      Sigma11 = sp_cov_mat(t1 = i, t2=i, sigma,c1,c2,nu,a,alpha,beta,delta=mydelta, k)
      
      Z1 = Z[1:sl]
      chol_st<-chol(Sigma11)
      inv_st = chol2inv(chol_st)
      loglik = -0.5*2*sum(log(diag(chol_st))) -0.5*t(Z1-mu1)%*%inv_st%*%(Z1-mu1)
      logL = logL + loglik
    } else{
      mu1  = mu[i]*rep(1 , times=sl)
      mu2 = mu[i-1]*rep(1 , times=sl)
      Sigma11 = sp_cov_mat(t1 = i, t2=i, sigma,c1,c2,nu,a,alpha,beta,delta=mydelta, k)
      Sigma22 = sp_cov_mat(t1 = i-1, t2=i-1, sigma,c1,c2,nu,a,alpha,beta,delta=mydelta, k)
      Sigma12 = Sigma21 =  sp_cov_mat(t1 = i, t2=i-1, sigma,c1,c2,nu,a,alpha,beta,delta=mydelta, k)
      
      Z1 = Z[(sl*(i-1)+1):(sl*i)]; Z2 = Z[(sl*(i-2)+1):(sl*(i-1))]
      mu_cond = mu1+ Sigma12%*%solve(Sigma22)%*%(Z2- mu2)
      Sigma_cond = Sigma11- Sigma12%*%solve(Sigma22)%*%Sigma21
      chol_st<-chol(Sigma_cond)
      inv_st = chol2inv(chol_st)
      
      loglik = -0.5*2*sum(log(diag(chol_st))) -0.5*t(Z1-mu_cond)%*%inv_st%*%(Z1-mu_cond)
      logL = logL + loglik 
    }
  }
  
  return(-logL)
}




logL_markov_cov = function(par, Z,sl,tl, k, hypo){
  
  sigma= 1; delta = mydelta;  nu=0.5; a= mya; alpha =myalpha; beta = mybeta
  if(hypo==0){
    c = par[1] #; a = par[2]; alpha= par[3]; beta = par[4]
    c1=c; c2=c
  } else if (hypo==1){
    c1= par[1]; c2=par[2]#;  a = par[3]; alpha = par[4]; beta = par[5]
  }
  
  museg1= museg2=0
  mu = ifelse(tseq > k,museg2, museg1)
  
  logL=0
  for (i in tl:1) {
    
    if(i==1){
      mu1  = mu[i]*rep(1 , times=sl)
      Sigma11 = sp_cov_mat(t1 = i, t2=i, sigma,c1,c2,nu,a,alpha,beta,delta=mydelta, k)
      
      Z1 = Z[1:sl]
      chol_st<-chol(Sigma11)
      inv_st = chol2inv(chol_st)
      loglik = -0.5*2*sum(log(diag(chol_st))) -0.5*t(Z1-mu1)%*%inv_st%*%(Z1-mu1)
      logL = logL + loglik
    } else{
      mu1  = mu[i]*rep(1 , times=sl)
      mu2 = mu[i-1]*rep(1 , times=sl)
      Sigma11 = sp_cov_mat(t1 = i, t2=i, sigma,c1,c2,nu, a,alpha,beta,delta=mydelta,k)
      Sigma22 = sp_cov_mat(t1 = i-1, t2=i-1, sigma,c1,c2,nu,a,alpha,beta,delta=mydelta, k)
      Sigma12 = Sigma21 =  sp_cov_mat(t1 = i, t2=i-1, sigma,c1,c2,nu,a,alpha,beta,delta=mydelta, k)
      
      Z1 = Z[(sl*(i-1)+1):(sl*i)]; Z2 = Z[(sl*(i-2)+1):(sl*(i-1))]
      mu_cond = mu1+ Sigma12%*%solve(Sigma22)%*%(Z2- mu2)
      Sigma_cond = Sigma11- Sigma12%*%solve(Sigma22)%*%Sigma21
      chol_st<-chol(Sigma_cond)
      inv_st = chol2inv(chol_st)
      
      loglik = -0.5*2*sum(log(diag(chol_st))) -0.5*t(Z1-mu_cond)%*%inv_st%*%(Z1-mu_cond)
      logL = logL + loglik 
    }
  }
  
  return(-logL)
}



start.time = Sys.time()
fit.optim_markov = optim(par=c(1, 0.4,0.6),logL_markov_mean,Z=Z,sl=sl,tl=tl, k=tl, hypo=0)
end.time = Sys.time()
end.time-start.time
fit.optim_markov

fit.optim_markov_mu = optim(par=c(1,0.5,0.5,0.7),logL_markov_mean,Z=Z,sl=sl,tl=tl, k=tl, hypo=0)

optim(par=c(0,2),logL_markov,Z=Z,sl=sl,tl=tl, k=25, hypo=1)


##covariance change
## likelihood ratio 
start.time =  Sys.time()
optimfull = optim(par=c(5,0.56, 0.51, 0.7),logL_markov_cov, Z=Z,sl=sl,tl=tl, k=tl, hypo=0)
loglikfull = -optimfull$value

lr_ts = vector(length=tl-1)
for (k in 2:(tl-2)) {
  
  optimalt = optim(par=c(5,3,0.51, 0.51, 0.7),logL_markov_cov,Z=Z,sl=sl,tl=tl, k=k, hypo=1)
  loglikalt = -optimalt$value
  
  lr_ts[k]  = 2*(loglikalt-loglikfull) 
}
end.time =  Sys.time()
end.time-start.time

plot(lr_ts, type="o", xlab="t", ylab = "LR", main="LRT (single changepoint in spatial dependence)")
which.max(lr_ts)


## mean change
optimfull = optim(par=c(1,0.5,0.7),logL_markov_mean,Z=Z,sl=sl,tl=tl, k=tl, hypo=0)
loglikfull = -optimfull$value
lr_ts = vector(length=tl-1)
for (k in 2:(tl-2)) {
  
  optimalt = optim(par=c(0.5,1,0.5,0.7),logL_markov_mean,Z=Z,sl=sl,tl=tl, k=k, hypo=1)
  loglikalt = -optimalt$value
  
  lr_ts[k]  = 2*(loglikalt-loglikfull) 
}
end.time =  Sys.time()
end.time-start.time

plot(lr_ts, type="o", xlab="t", ylab = "LR", main="LRT (single changepoint in spatial mean)")
which.max(lr_ts)








### Testing for change in spatial process in the true model
## keeping nuisance parameters fixed

spcov = st_covmat_func(c1=5,c2=3,nu=0.5, a=mya, alpha=myalpha, beta = mybeta, delta= mydelta)

nsim<-1
mysim<-rmvnorm(n=nsim,mean = c(rep(0, times=sl*tl/2), rep(0, times=sl*tl/2)), sigma=spcov,method="chol")
Z<-c(mysim)




##mean change
start.time =  Sys.time()
optimfull = optim(par=c(1, 1, 0.5, 0.5),logL_markov_mean,Z=Z,sl=sl,tl=tl, k=tl, hypo=0)
loglikfull = -optimfull$value

lr_ts = vector(length=tl-1)
for (k in 2:(tl-2)) {
  
  optimalt = optim(par=c(0.5,1, 1, 0.5, 0.5),logL_markov_mean,Z=Z,sl=sl,tl=tl, k=k, hypo=1)
  loglikalt = -optimalt$value
  
  lr_ts[k]  = 2*(loglikalt-loglikfull) 
}
end.time =  Sys.time()
end.time-start.time

plot(lr_ts, type="o", xlab="t", ylab = "LR", main="LRT (single changepoint in spatial mean under misspecification)")
which.max(lr_ts)




##cov change
start.time =  Sys.time()
optimfull = optim(par=c(5),logL_markov_cov,Z=Z,sl=sl,tl=tl, k=tl, hypo=0)
loglikfull = -optimfull$value

lr_ts = vector(length=tl-1)
for (k in 2:(tl-2)) {
  
  optimalt = optim(par=c(5,3),logL_markov_cov,Z=Z,sl=sl,tl=tl, k=k, hypo=1)
  loglikalt = -optimalt$value
  
  lr_ts[k]  = 2*(loglikalt-loglikfull) 
}
end.time =  Sys.time()
end.time-start.time

plot(lr_ts, type="o", xlab="t", ylab = "LR", main="LRT (single changepoint in spatial dep)")
which.max(lr_ts)











### Repeating simulation for change in dependence   ####


### Our method-markov
mya = 0.25
## Function for computing LR test statistic
LR_func = function(Z, tl){
  optimnull = optim(par=c(5),logL_markov_cov,Z=Z,sl=sl,tl=tl, k=tl, hypo=0)
  loglikfull = -optimnull$value
  
  lr_ts = vector(length=tl-1)
  for (k in 3:(tl-2)) {
    
    optimalt = optim(par=c(5,5),logL_markov_cov,Z=Z,sl=sl,tl=tl, k=k, hypo=1)
    loglikalt = -optimalt$value
    
    lr_ts[k]  = 2*(loglikalt-loglikfull) 
  }
  return(lr_ts)
}


library(doParallel)
registerDoParallel(cores=22)

mya = 0.901; myalpha = 0.75; mybeta = 0.61; mydelta =0
sim_size = 20
ss = Sys.time()
LR_null = foreach(s = 1:sim_size, .combine =cbind) %dopar% {
  spcov = st_covmat_func(c1=5,c2=5,nu=0.5, a=mya, alpha=myalpha, beta = mybeta, delta= mydelta)
  #spcov = st_covmat_func(c1=1/140,c2=1/140,nu=0.5, a=0.901, alpha=0.75, beta = 0.61, delta= 0)
  mysim<-rmvnorm(n=1,mean = c(rep(0, times=sl*tl/2), rep(0, times=sl*tl/2)), sigma=spcov,method="chol")
  Z<-c(mysim)
  
  Lr = (LR_func(Z,tl))
  return(Lr)
}
ee =Sys.time()
ee-ss
apply(LR_null, 2, max)
th_markov_dep = quantile(apply(LR_null, 2, max), 0.95)
th_markov_dep = 8.759

set.seed(47)
data_sim = list()
for (s in 1:sim_size) {
  spcov = st_covmat_func(c1=1/0.4,c2=1/0.45,nu=0.5, a=mya, alpha=myalpha, beta = mybeta, delta= mydelta)
  mysim<-rmvnorm(n=1,mean = c(rep(0, times=sl*tl/2), rep(0, times=sl*tl/2)), sigma=spcov,method="chol")
  data_sim[[s]] = c(mysim)
}

system.time({Eval_markov = foreach(s = 1:sim_size, .combine = rbind) %dopar% {
 Z<-data_sim[[s]]
  Lr = LR_func(Z,tl)
  Lr_max = max(Lr)
  if(Lr_max > th_markov_dep){
    cpt_markov = which.max(Lr)} else {cpt_markov = NA}
  return(cbind(cpt_markov, Lr_max, matrix(Lr, nrow = 1)))
}
})
mean(Eval_markov[,"Lr_max"] > th_markov_dep) #power
sum(abs(Eval_markov[,"cpt_markov"] - t0)<=2, na.rm=T)/sim_size 

plot(Eval_markov[3,-c(1,2)], type="o")
Eval_markov[,c(1,2)]

### comparison with full likelihood 
LR_full_func = function(Z){
  # nlm is faster than optim
  optimnull = nlm(logL_st,p=c(5),Z=Z,X=X,sl=sl,tl=tl, k=tl, hypo=0, steptol=1e-2,gradtol=1e-2)
  loglikfull = -optimnull$minimum
  
  #optimnull = optim(par=c(5),logL_st,Z=Z,X=X,sl=sl,tl=tl, k=tl, hypo=0)
  #loglikfull = -optimnull$value
  
  lr_ts = vector(length=tl-1)
  for (k in 3:(tl-2)) {
    
    optimalt = nlm(logL_st,p=c(5,3.5),Z=Z,X=X,sl=sl,tl=tl, k=k, hypo=1, steptol=1e-2,gradtol=1e-2)
    loglikalt = -optimalt$minimum
    #optimalt = optim(par=c(5,4),logL_st,Z=Z,X=X,sl=sl,tl=tl, k=k, hypo=1)
    #loglikalt = -optimalt$value
    
    lr_ts[k]  = 2*(loglikalt-loglikfull) 
  }
  return(lr_ts)
}


sim_size = 20
LR_full_null = foreach(s = 1:sim_size, .combine =cbind) %dopar% {
  spcov = st_covmat_func(c1=5,c2=5,nu=0.5, a=mya, alpha=myalpha, beta = mybeta, delta= mydelta)
  mysim<-rmvnorm(n=1,mean = c(rep(0, times=sl*tl/2), rep(0, times=sl*tl/2)), sigma=spcov,method="chol")
  Z<-c(mysim)
  Lr = LR_full_func(Z)
  return(Lr)
}


apply(LR_full_null, 2, max)
th_full_dep = quantile(apply(LR_full_null, 2, max), 0.95)
th_full_dep = 8.1

start.time =  Sys.time()
Eval_full = foreach(s = 1:sim_size, .combine = rbind) %dopar% {
  Z<-data_sim[[s]]
  Lr = LR_full_func(Z)
  Lr_max = max(Lr)
  if(Lr_max > th_full_dep){
    cpt_full = which.max(Lr)} else {cpt_full = NA}
  return(cbind(cpt_full, Lr_max, matrix(Lr,nrow=1)))
}
end.time =  Sys.time()
end.time-start.time

mean(Eval_full[,"Lr_max"] > th_full_dep) #power
sum(abs(Eval_full[,"cpt_full"] - t0)<=2, na.rm=T)/sim_size #accuracy TDR
Eval_full[,c(1,2)]




plot(Eval_markov[4,-c(1,2)], ylim=c(0,40), ylab="LR", type="o", col= "green")
lines(Eval_full[4,-c(1,2)], type="o",col="blue")
legend("topright", legend=c("Markov lik", "Full lik"),
       col=c("green", "blue"), lty=1, lwd=2)




### comparison with independent segments with pairwise likelihood
sim_size = 20

d =2 #two dimensional spatial coordinates
phi<-function(h,c,nu)
{
  h[h==0]<-1e-10
  num1<-(2^(1-nu))/gamma(nu)
  num2<-((h^(1/2))*c)^nu
  num3<-besselK(x=((h^(1/2))*c), nu=nu)
  return(num1*num2*num3)
}



st_cov = function( h,t, sigma  , c  , nu, a , alpha , beta ){
  psi_fun = psi(abs(t)^2,a , alpha, beta)
  phi_fun = phi((h^2/ psi_fun),c, nu)
  return((sigma^2*phi_fun) / (psi_fun^(d/2)))
  
}

tseq = 1:tl
tmat<-outer(tseq,tseq,function(t1,t2) (abs(t1-t2)))
hmat_st =tmat_st= matrix(NA,ncol = sl*tl,nrow=sl*tl)
for(i in 1:tl)
{
  for(j in i:tl)
  {
    hmat_st[((i-1)*sl+1):((i-1)*sl+sl),((j-1)*sl+1):((j-1)*sl+sl)] = hmat_st[((j-1)*sl+1):((j-1)*sl+sl),((i-1)*sl+1):((i-1)*sl+sl)]<-dist.sp
    tmat_st[((i-1)*sl+1):((i-1)*sl+sl),((j-1)*sl+1):((j-1)*sl+sl)] = tmat_st[((j-1)*sl+1):((j-1)*sl+sl),((i-1)*sl+1):((i-1)*sl+sl)]<-tmat[i,j]
  }
}

logL_pair = function(par,Z,sl,tl){
  logL = 0
  sigma= 1; c = par[1]; a = mya; alpha= myalpha; beta = mybeta; beta_gls=0; nu=0.5
  
  hmat_st = hmat_st[1:(sl*tl), 1:(sl*tl)]; tmat_st = tmat_st[1:(sl*tl), 1:(sl*tl)]
 
  hmat_st[lower.tri(hmat_st, diag = T)]=NA; tmat_st[lower.tri(tmat_st, diag = T)]=NA # to consider unique pairs
  Index_lag = which( hmat_st<0.4 & tmat_st<3, arr.ind=TRUE) # restricting pairs according to spatial lag and temporal lag
  
  for (i in 1:dim(Index_lag)[1]) {
    biv_cov_mat = sigma^2*diag(2)
    biv_cov_mat[1,2] = biv_cov_mat[2,1] = st_cov(h = hmat_st[Index_lag[i,1],Index_lag[i,2]], 
                                                 t= tmat_st[Index_lag[i,1],Index_lag[i,2]], 
                                                 sigma, c , nu, a, alpha, beta) 
    inv_biv = solve(biv_cov_mat)
    Z_pair = c(Z[Index_lag[i,1]], Z[Index_lag[i,2]]); X_pair = matrix(1, nrow=2, ncol=1)
    logL = logL + -0.5*log(det(biv_cov_mat)) -0.5*t(Z_pair-X_pair%*%beta_gls)%*%inv_biv%*%(Z_pair-X_pair%*%beta_gls)
    
  }
  return(-logL)
}



LR_pair_func = function(Z){
  optimnull = optim(par=c(5),logL_pair,Z=Z,sl=sl,tl=tl)
  loglikfull = -optimnull$value
  
  lr_ts = vector(length=tl-1)
  for (k in 3:(tl-2)) {
    
    Z1 = Z[1:(sl*k)]
    optim1 = optim(par=c(5),logL_pair,Z=Z1,sl=sl,tl=k)
    loglik1 = -optim1$value
    #par_ini_seg1 = optim1$par
    
    Z2 = Z[-c(1:(sl*k))]
    optim2 = optim(par=c(4),logL_pair,Z=Z2,sl=sl,tl=tl-k)
    loglik2 = -optim2$value
    #par_ini_seg2 = optim2$par
    
    lr_ts[k]  = 2*(loglik1+loglik2-loglikfull) 
  }
  return(lr_ts)
}
sim_size = 20
LR_pair_null = foreach(s = 1:sim_size, .combine =cbind) %dopar% {
  spcov = st_covmat_func(c1=5,c2=5,nu=0.5, a=mya, alpha=myalpha, beta = mybeta, delta= mydelta)
  mysim<-rmvnorm(n=1,mean = c(rep(0, times=sl*tl/2), rep(0, times=sl*tl/2)), sigma=spcov,method="chol")
  Z<-c(mysim)
  Lr = LR_pair_func(Z)
  return(Lr)
}


#apply(LR_pair_null, 2, max)
th_pair = quantile(apply(LR_pair_null, 2, max), 0.90)
th_pair =500

start.time =  Sys.time()
Eval_pair = foreach(s = 1:sim_size, .combine = rbind) %dopar% {
  Z = data_sim[[s]]
  Lr = LR_pair_func(Z)
  Lr_max = max(Lr)
  if(Lr_max > th_pair){
    cpt_pair = which.max(Lr)} else {cpt_pair = NA}
  return(cbind(cpt_pair, Lr_max))
}
end.time =  Sys.time()
end.time-start.time

mean(Eval_pair[,"Lr_max"] > th_pair) #power
sum(abs(Eval_pair[,"cpt_pair"] - t0)<=1, na.rm=T)/sim_size #accuracy










##### Motivation - issues in assumption of independence between segments
phi<-function(h,c,nu)
{
  h[h==0]<-1e-10
  num1<-(2^(1-nu))/gamma(nu)
  num2<-((h^(1/2))*c)^nu
  num3<-besselK(x=((h^(1/2))*c), nu=nu)
  return(num1*num2*num3)
}



st_cov = function( h,t, sigma  , c  , nu, a , alpha , beta ){
  d=2 # 2 dimensional spatial coordinates
  psi_fun = psi(abs(t)^2,a , alpha, beta)
  phi_fun = phi((h^2/ psi_fun),c, nu)
  return((sigma^2*phi_fun) / (psi_fun^(d/2)))
  
}
stcovmat_function = function(tl,sigma,c, nu, a, alpha, beta){
  st_cov_mat=  matrix(NA,ncol = sl*tl,nrow=sl*tl)
  tseq = 1:tl
  tmat<-outer(tseq,tseq,function(t1,t2) (abs(t1-t2)))
  for(i in 1:tl)
  {
    for(j in i:tl)
    {
      tmp1<-st_cov(h = uniq.dist, t= tmat[i,j], sigma , c , nu , a, alpha, beta)
      temp[mat.index[,-3]]<-tmp1[mat.index[,3]]
      st_cov_mat[((i-1)*sl+1):((i-1)*sl+sl),((j-1)*sl+1):((j-1)*sl+sl)]<-st_cov_mat[((j-1)*sl+1):((j-1)*sl+sl),((i-1)*sl+1):((i-1)*sl+sl)]<-temp
    }
  }
  return(st_cov_mat)
}

stcovmat_lag = function(tlag,sigma, c, nu, a, alpha, beta){
  st_cov_mat=  matrix(NA,ncol = sl,nrow=sl)
  
  tmp1<-st_cov(h = uniq.dist, t= tlag, sigma , c , nu , a, alpha, beta)
  temp[mat.index[,-3]]<-tmp1[mat.index[,3]]
  st_cov_mat<-temp
  
  return(st_cov_mat)
}

## data is simulated for piecewise stationary segment (no change)
set.seed(21)
set.seed(22)

mysim1<-rmvnorm(n=1,mean = c(rep(0, times=sl*tl/2)), 
                sigma=stcovmat_function(tl/2, c=5, sigma=1, nu=0.5, a= mya, alpha=myalpha, beta=mybeta),method="chol")
mysim2<-rmvnorm(n=1,mean = c(rep(0, times=sl*tl/2)), 
                sigma=stcovmat_function(tl/2, c=5, sigma=1, nu=0.5, a= mya, alpha=myalpha, beta=mybeta),method="chol")

Z = c(mysim1, mysim2) #Z




logL_markov_ind = function(par, Z,sl,tl){
  sigma= 1; c = par[1]; nu = 0.5; a =mya; alpha=myalpha; beta=mybeta
  
  Sigma11 = stcovmat_lag(tlag = 0, sigma,c,nu,a,alpha,beta)
  Sigma22 = stcovmat_lag(tlag = 0, sigma,c,nu,a,alpha,beta)
  Sigma12 = Sigma21 = stcovmat_lag(tlag = 1, sigma,c,nu,a,alpha,beta)
  mu=0
  mu1= mu2 = mu*rep(1 , times=sl)
  
  logL=0
  for (i in tl:1) {
    
    if(i==1){
      Z1 = Z[1:sl]
      chol_st<-chol(Sigma11)
      inv_st = chol2inv(chol_st)
      loglik = -0.5*2*sum(log(diag(chol_st))) -0.5*t(Z1-mu1)%*%inv_st%*%(Z1-mu1)
      logL = logL + loglik
    } else{
      Z1 = Z[(sl*(i-1)+1):(sl*i)]; Z2 = Z[(sl*(i-2)+1):(sl*(i-1))]
      mu_cond = mu1+ Sigma12%*%solve(Sigma22)%*%(Z2- mu2)
      Sigma_cond = Sigma11- Sigma12%*%solve(Sigma22)%*%Sigma21
      chol_st<-chol(Sigma_cond)
      inv_st = chol2inv(chol_st)
      
      loglik = -0.5*2*sum(log(diag(chol_st))) -0.5*t(Z1-mu_cond)%*%inv_st%*%(Z1-mu_cond)
      logL = logL + loglik 
    }
  }
  
  return(-logL)
}


optimfull = optim(par=c(5), logL_markov_ind,Z=Z,sl=sl,tl=tl)
loglikfull = -optimfull$value

lr_ts = vector(length=tl-1)
for (k in 1:(tl-1)) {
  Z1 = Z[1:(sl*k)]
  optim1 = optim(par=c(5),  logL_markov_ind,Z=Z1,sl=sl,tl=k)
  loglik1 = -optim1$value
  
  Z2 = Z[-c(1:(sl*k))]; X2 = Z[-c(1:(sl*k))]
  optim2 =optim(par=c(5), logL_markov_ind,Z=Z2,sl=sl,tl=tl-k)
  loglik2 = -optim2$value
  
  
  lr_ts[k]  = 2*(loglik1+loglik2-loglikfull) 
}


plot(lr_ts, type="o", ylab = "LR", xlab= "t", main="LRT assuming independent segments (No changepoint)")

which.max(lr_ts)

lr_mar = LR_func(Z, tl)
plot(lr_mar, type="o", ylab = "LR", xlab= "t", main="LRT using proposed method (No changepoint)", ylim=c(0,10))




## comparison of LRT in case of dependence and a changepoint

