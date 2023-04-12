##### Ireland wind data #################
### Source: The Irish Meterological Service (https://www.met.ie/climate/available-data) ########
## wind speed in knots for 22 synoptic stations ###

stationdetails = read.csv("/home/agarwalg/Ireland Wind Data/StationDetails.csv", header=T)
head(stationdetails)
stationdetails$name


## reading files for all the stations - dates and windspeed only
temp = list.files(path="/home/agarwalg/Ireland Wind Data", pattern="*csv",full.names = T); temp = temp[-23]
myfiles = lapply(temp, read.csv, header=T, skip = 24, 
                 colClasses = c("character", rep("NULL", 9), "numeric",rep("NULL", 13) ))
stationcodes = as.numeric(substr(temp, 37, 40))
stations = stationdetails[stationdetails$station.name %in% stationcodes, ]
sort(stations$station.name)

## converting to date format (returns year-month-date)
wdsplist = list(); myfiles_hist = list(); wdsplist_hist = list()
for (i in 1:length(myfiles)) {
  myfiles[[i]]$date = as.Date(myfiles[[i]]$date, format = "%d-%b-%Y")
  myfiles_hist[[i]] = myfiles[[i]][ myfiles[[i]]$date >= "2015-01-01" & myfiles[[i]]$date <= "2019-12-31",]
  wdsplist_hist[[i]] = myfiles_hist[[i]]$wdsp
  
  myfiles[[i]] = myfiles[[i]][ myfiles[[i]]$date >= "2020-01-01" & myfiles[[i]]$date <= "2021-12-31",]
  wdsplist[[i]] = myfiles[[i]]$wdsp
}



library(rlist)
wind.data = list.cbind(wdsplist)
colnames(wind.data) = c(stationcodes)

wdsplist_hist[[14]][1824:1826] = NA
wind.histdata = list.cbind(wdsplist_hist)
colnames(wind.histdata) = c(stationcodes)
wind.histdata = data.frame(date = myfiles_hist[[1]]$date,  wind.histdata, check.names=F)
head(wind.histdata)

#reaarange data according to station data
wind.data = wind.data[,match(stations$station.name,colnames(wind.data))]
wind.data = data.frame(date = myfiles[[1]]$date,  wind.data, check.names=F)
head(wind.data)

head(stations)

## sum(is.na(wind.data)) 3 missing observations
library(ggmap)
load("/home/agarwalg/wind data/map_ire.RData")

wind.plot = data.frame(x =stations$longitude, y = stations$latitude, wind = apply(0.5148 *wind.data[,-1],2, mean, na.rm=T) )
dim(wind.plot)
library(viridis)
ggmap(map_ire) + geom_point(data=wind.plot, aes(x,y, col = wind), pch=16, size=4) +scale_color_viridis("Wind speed (m/s)", option="C")+xlab("Longitude") + ylab("Latitude")

## converting wind speed to m/s and taking a square root transformation
windsqrtms = windsqrt = sqrt(0.5148 * as.matrix(wind.data[,-1])) 

date_at = c(which(wind.data$date %in% seq(min(wind.data$date), max(wind.data$date), by="quarter")),length(wind.data$date))
plot(wind.data$date,apply(windsqrt, 1, mean), type="l", xaxt = "n",
     ylab="Average Wind Speed (sqrt)", xlab="date")
axis(1, at = c(seq(min(wind.data$date), max(wind.data$date), by="quarter"), max(wind.data$date)),
     labels = format(wind.data$date, "%d-%b-%y")[date_at])

#axis.Date(1, x= wind.data$date, format(wind.data$date, "%b-%y"))

## check average wind speed at locations
apply(as.matrix(windsqrt), 2, mean, na.rm=T)
quantile(apply(as.matrix(windsqrt), 2, mean, na.rm=T))




## plot time series of two locations lowest (Ballyhaise) and highest (Malin Head)
head(windsqrt)
plot(wind.data$date,windsqrt[,"675"], type="l", ylim=c(0.5,4.5), ylab = "Wind speed (sqrt)", 
     main="Wind speeds at two locations", xaxt="n", xlab="date")
lines(wind.data$date, windsqrt[,"1575"], lty=2)
axis(1, at = c(seq(min(wind.data$date), max(wind.data$date), by="quarter"), max(wind.data$date)),
     labels = format(wind.data$date, "%d-%b-%y")[date_at])

#plot(wind.data$date,0.5148*wind.data[,"675"], type="l", ylim=c(0,22), ylab = "Wind speed (m/s)", 
#     main="Wind speeds at two locations", xaxt="n", xlab="date")
#lines(wind.data$date, 0.5148*wind.data[,"275"], lty=2)
#axis(1, at = c(seq(min(wind.data$date), max(wind.data$date), by="quarter"), max(wind.data$date)),
#     labels = format(wind.data$date, "%d-%b-%y")[date_at])

## clear evidence of spatial and temporal dependencies


par(mfrow=c(2,1))
hist(as.matrix(0.5148 *wind.data[,-1]), xlab = "", main="Histogram of wind speed (before transformation)")
hist(windsqrt, xlab ="" , main="Histogram of wind speed (after sqrt transformation)")
par(mfrow=c(1,1))

#for (i in  1:22) plot(windsqrt[,i], type="l", main = i)



#### seasonality
Jday = 1:366

windsqrt = windsqrt - mean(windsqrt, na.rm = T)
wind.data$jday = as.numeric(format(wind.data$date, '%j'))
#daymeans = sapply(split(windsqrt, wind.data$jday), mean)

#plot(daymeans ~ Jday, xlab="day of the year")
#lines(lowess(daymeans ~ Jday, f = 0.14), lwd=2)


## seasonality using historical data
windsqrt.hist =  sqrt(0.5148 * as.matrix(wind.histdata[,-1])) 
windsqrt.hist =  windsqrt.hist - mean(windsqrt.hist, na.rm = T)
wind.histdata$jday = as.numeric(format(wind.histdata$date, '%j'))

daymeans = sapply(split(windsqrt.hist, wind.histdata$jday), mean, na.rm=T)

plot(daymeans ~ Jday, xlab="day of the year", main="Seasonality effect")
lines(lowess(daymeans ~ Jday, f= 0.2, delta = 0.01*range(Jday)), lwd=2)
lines(lowess(c(daymeans,daymeans,daymeans) ~ c(Jday-366,Jday,Jday+366), f= 0.066)$y[366+Jday], lwd=2, col="black")



0.2*366

73/(366*3)
### removing seasonality
seasonalwind = lowess(c(daymeans,daymeans,daymeans) ~ c(Jday-366,Jday,Jday+366), f= 0.066)$y[366+Jday]
velocity = apply(windsqrt, 2, function(x) { x - seasonalwind })

plot(apply(velocity, 1, mean), type="l")

# fig 3, but not really yet...
dists = spDists(cbind(stations$longitude, stations$latitude), longlat=TRUE)
corv = cor(velocity, use = "complete.obs")
sel = !(as.vector(dists) == 0)
plot(as.vector(corv[sel]) ~ as.vector(dists[sel]),
     xlim = c(0,500), ylim = c(.4, 1), xlab = "distance (km.)", pch=3,
     ylab = "correlation") 

data_h_temp  = data.frame(h=as.vector(dists[sel]), cor =as.vector(corv[sel]))
C_h = function(h, sigma, c) { sigma*exp(-c*h)}
fit.sp = nls(cor ~ C_h(h, sigma, c), data=data_h_temp, start=list(sigma=0.9, c=0.00134))
## sigma^2 = 0.9606, c = 0.00113

xdiscr = 1:500
# add correlation model:
lines(xdiscr, .9606 * exp(- .00113 * xdiscr))



uniq.dist<-sort(unique(c(dists)))
##### Unique indexing #####
mat.index<-NULL
for(i in 1:length(uniq.dist)){
  mat.index[[i]]<-cbind(which(dists==uniq.dist[i],arr.ind = T),i)
}
mat.index<-do.call(rbind,mat.index)


### empirical temporal correlations
emp_temp = apply(velocity, 2, acf, plot= F, na.action = na.pass)
emp_lagall = sapply(emp_temp, function(x) x$acf[2:25])
#plot(apply(emp_lagall,1,mean), ylim= c(-0.1,1), type="n")
#abline(b=apply(emp_lagall,1,mean))

emp_lag = sapply(emp_temp, function(x) x$acf[2:4])
apply(emp_lag,1,mean)
## empirical correlations of 0.5151, 0.2330, 0.1593, at lags 1, 2, 3.

par(mfrow=c(2,3))
set.seed(23)
for(i in sample(1:22, 6)){myacf = acf(velocity[,i],plot=F,  na.action = na.pass)
plot(myacf[1:15], ylim=c(0,1), main = stations$name[i])}
par(mfrow=c(1,1))
dim(velocity)


data_emp_temp  = data.frame(u=1:3, c =apply(emp_lag,1,mean) )
C_u = function(u, a, alpha) { (a*u^(2*alpha)+1)^(-1)}
fit.tempcor = nls(c ~ C_u(u, a, alpha), data=data_emp_temp, start=list(a=0.9, alpha=0.7))
mya = coef(fit.tempcor)[1]
myalpha =  coef(fit.tempcor)[2] 






### log likelihood using markov property
my.matern<-function(h,c,sigma,nu)
{
  h[h==0]<-1e-10
  num1<-(sigma^2)*(2^(1-nu))/gamma(nu)
  num2<-(h*c)^nu
  num3<-besselK(x=(h*c), nu=nu)
  return(num1*num2*num3)
}

psi<-function(t,a,alpha,beta)
{
  temp<-(a*(t^alpha)+1)^beta
  return(temp)
}

sp_cov_mat = function(t1,t2, sl,tl, sigma, c1, c2, nu,a,alpha,beta,delta, k ){
  
  sp_mat<-matrix(NA,ncol = sl,nrow = sl)
  
  tseq = 1:tl
  c_vec = ifelse(tseq>k, c2,c1)
  nu_vec = ifelse(tseq>k, nu,nu)
  c_mean = mean(c_vec)
  
  term1 = gamma((nu_vec[t1]+nu_vec[t2])/2)/ sqrt(gamma(nu_vec[t1])*gamma(nu_vec[t2]))
  term2 = 1/sqrt((c_vec[t1]^2)*(c_vec[t2]^2))
  term3 = fc_t1t2 = 1/(psi(t=(t1-t2)^2,a=a,alpha = alpha,beta=beta)/(c_mean^2) +
                         ((1/(c_vec[t1]^2))+(1/(c_vec[t2]^2)))/2 - psi(t=0,a=a,alpha = alpha,beta=beta)/(c_mean^2))
  term4 = 1/(psi(t=(t1-t2)^2,a=a,alpha = alpha,beta=delta)) #pure temporal cor
  
  tmp1<-term1*term2*term3*my.matern(h=uniq.dist,
                                    c=fc_t1t2^(1/2),
                                    nu=(nu_vec[t1]+nu_vec[t2])/2,
                                    sigma=1)
  sp_mat[mat.index[,-3]]<-tmp1[mat.index[,3]]
  return(sp_mat)
}



logL_markov = function(par, Z,sl,tl, k, hypo){
  
  mydelta=0
  sigma= 0.98; delta = 0; nu=0.5; a= 0.958; alpha =0.82; beta = 0.61
  if(hypo==0){
    c = par[1]
    c1=c; c2=c
  } else if (hypo==1){
    c1= par[1]; c2=par[2]
  }
  tseq = 1:tl
  museg1= museg2=0
  mu = ifelse(tseq > k,museg2, museg1)
  

  
  logL=0
  for (i in tl:2) {
    
    if(i==1){
      mu1  = mu[i]*rep(1 , times=sl)
      Sigma11 = sp_cov_mat(t1 = i, t2=i, sl= sl,tl=tl, sigma,c1,c2,nu,a,alpha,beta,mydelta, k)
      
      Z1 = Z[1:sl]
      chol_st<-chol(Sigma11)
      inv_st = chol2inv(chol_st)
      loglik = -0.5*2*sum(log(diag(chol_st))) -0.5*t(Z1-mu1)%*%inv_st%*%(Z1-mu1)
      logL = logL + loglik
    } else{
      mu1  = mu[i]*rep(1 , times=sl)
      mu2 = mu[i-1]*rep(1 , times=sl)
      Sigma11 = sp_cov_mat(t1 = i, t2=i, sl=sl,tl=tl, sigma,c1,c2,nu, a,alpha,beta,mydelta,k)
      Sigma22 = sp_cov_mat(t1 = i-1, t2=i-1,sl=sl,tl=tl, sigma,c1,c2,nu,a,alpha,beta,mydelta, k)
      Sigma12 = Sigma21 =  sp_cov_mat(t1 = i, t2=i-1,sl=sl,tl=tl, sigma,c1,c2,nu,a,alpha,beta,mydelta, k)
      
      Z1 = Z[(sl*(i-1)+1):(sl*i)]; Z2 = Z[(sl*(i-2)+1):(sl*(i-1))]
      missZ1 = which(is.na(Z1)); missZ2 = which(is.na(Z2))
      
      if(length(missZ1)==0){missZ1= length(Z1)+1}; if(length(missZ2)==0){missZ2= length(Z2)+1} # don't delete any observation if there are no missing obervations
      ### omit the missing observations and the corresponding cov entries
      Z1 = Z1[-missZ1]; Z2 = Z2[-missZ2]; mu1 = mu1[-missZ1]; mu2= mu2[-missZ2]
      Sigma11 = Sigma11[-missZ1, -missZ1]; Sigma22 = Sigma22[-missZ2, -missZ2];Sigma12 = Sigma12[-missZ1, -missZ2]; Sigma21= Sigma21[-missZ2, -missZ1]
      
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


logL_markov = function(par, Z,sl,tl, k, hypo){
  
  mydelta=0
  sigma= 0.98; delta = 0; nu=0.5; a= 0.958; alpha =0.82; beta = 0.61
  if(hypo==0){
    museg1= museg2 = par[1]; c = par[2]
    c1=c; c2=c
  } else if (hypo==1){
    museg1 = par[1]; museg2 = par[2];c1= par[3]; c2=par[4]
  }
  tseq = 1:tl
  mu = ifelse(tseq > k,museg2, museg1)
  
  logL=0
  for (i in tl:2) {
    
    if(i==1){
      mu1  = mu[i]*rep(1 , times=sl)
      Sigma11 = sp_cov_mat(t1 = i, t2=i, sl= sl,tl=tl, sigma,c1,c2,nu,a,alpha,beta,mydelta, k)
      
      Z1 = Z[1:sl]
      chol_st<-chol(Sigma11)
      inv_st = chol2inv(chol_st)
      loglik = -0.5*2*sum(log(diag(chol_st))) -0.5*t(Z1-mu1)%*%inv_st%*%(Z1-mu1)
      logL = logL + loglik
    } else{
      mu1  = mu[i]*rep(1 , times=sl)
      mu2 = mu[i-1]*rep(1 , times=sl)
      Sigma11 = sp_cov_mat(t1 = i, t2=i, sl=sl,tl=tl, sigma,c1,c2,nu, a,alpha,beta,mydelta,k)
      Sigma22 = sp_cov_mat(t1 = i-1, t2=i-1,sl=sl,tl=tl, sigma,c1,c2,nu,a,alpha,beta,mydelta, k)
      Sigma12 = Sigma21 =  sp_cov_mat(t1 = i, t2=i-1,sl=sl,tl=tl, sigma,c1,c2,nu,a,alpha,beta,mydelta, k)
      
      Z1 = Z[(sl*(i-1)+1):(sl*i)]; Z2 = Z[(sl*(i-2)+1):(sl*(i-1))]
      missZ1 = which(is.na(Z1)); missZ2 = which(is.na(Z2))
      
      if(length(missZ1)==0){missZ1= length(Z1)+1}; if(length(missZ2)==0){missZ2= length(Z2)+1} # don't delete any observation if there are no missing obervations
      ### omit the missing observations and the corresponding cov entries
      Z1 = Z1[-missZ1]; Z2 = Z2[-missZ2]; mu1 = mu1[-missZ1]; mu2= mu2[-missZ2]
      Sigma11 = Sigma11[-missZ1, -missZ1]; Sigma22 = Sigma22[-missZ2, -missZ2];Sigma12 = Sigma12[-missZ1, -missZ2]; Sigma21= Sigma21[-missZ2, -missZ1]
      
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



dim(velocity)
Sn = sl=   dim(velocity)[2]; Tn = tl = dim(velocity)[1] 
##covariance change
## likelihood ratio 
Z = c(velocity)
# a, sigma, nu, c

LR_func = function(Z, tl) {
  optimfull = optim(par=c(2, 1/100),logL_markov, Z= Z ,sl=sl,tl=tl, k=Tn, hypo=0,
                   method = "L-BFGS-B",lower = c(1,1/400), upper=c(3, 1/50) )
  optimfull = optim(par=c(0, 1/100),logL_markov, Z= Z ,sl=sl,tl=tl, k=Tn, hypo=0 )
  
  nlm(logL_markov, p =c(0,1/300),  Z= Z ,sl=sl,tl=tl, k=Tn, hypo=0)
  
  loglikfull = -optimfull$value
  
  lr_ts = vector(length=tl-1)
  for (k in 2:(tl-2)) {
    
    optimalt = optim(par=c(0,0,1/300, 1/300),logL_markov,Z= Z ,sl=sl,tl=tl, k=k, hypo=1)
    
     optimalt = optim(par=c(2,2,1/300, 1/300),logL_markov,Z= Z ,sl=sl,tl=tl, k=k, hypo=1,
                    method = "L-BFGS-B", lower = c(0,0,1/400,1/400, upper=c(4,4,1/50, 1/50)))
    loglikalt = -optimalt$value
    
    lr_ts[k]  = 2*(loglikalt-loglikfull) 
  }
  return(lr_ts)
}


lr_ts = LR_func(Z = c(velocity), tl)
plot(lr_ts, ylim=c(0,100),type="o", xlab="t", ylab = "LR", main="LRT (changepoint in spatial process)")
which.max(lr_ts)

library(doParallel)
registerDoParallel(cores=20)
start.time = Sys.time()
LR_func = function(Z, tl) {
  
optimfull = optim(par=c(0,1/100),logL_markov, Z= Z ,sl=sl,tl=tl, k=Tn, hypo=0)
loglikfull = -optimfull$value
lr_ts = rep(0, tl-1)
lr_ts = foreach(k=2:(tl-2), .combine = c) %dopar% {
  optimalt = optim(par=c(0,0,1/100, 1/100),logL_markov,Z= Z ,sl=sl,tl=tl, k=k, hypo=1)
  loglikalt = -optimalt$value
  2*(loglikalt-loglikfull)
}
return(lr_ts)
}
end.time = Sys.time()
end.time-start.time

#lr_ts = c(0,lr_ts,0,0)
plot(wind.data$date, lr_ts,type="o", xlab="t", ylab = "LR", main="LRT (changepoint in spatial process)", xaxt="n")
axis(1, at = c(seq(min(wind.data$date), max(wind.data$date), by="quarter"), max(wind.data$date)),
     labels = format(wind.data$date, "%d-%b-%y")[date_at])
abline(v= wind.data$date[572], col="red")
which.max(lr_ts)
max(lr_ts)
plot(lr_ts)
wind.data$date[572]
## 25 July 21

#### finding the threshold ######
#### Monte carlo simuation from a Gaussian process with no change
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

st_covmat_func =  function(c1,c2,nu, a, alpha, beta, delta, sigma, t0){
  tseq = 1:tl
  c_vec = ifelse(tseq>t0, c2,c1)
  nu_vec = ifelse(tseq>t0, nu,nu)
  
  myf_c<-f_cvec(t=tseq,c_vec = c_vec,a,alpha, beta)
  myzeta<-zeta(t=tseq,nu_vec = nu_vec,c_vec = c_vec,a=a,alpha=alpha, beta=beta, delta=delta)
  mynus<-outer(nu_vec,nu_vec,function(a,b) (a+b)/2)
  
  spcov<-matrix(NA,ncol = sl*tl,nrow=sl*tl)
  temp<-matrix(NA,ncol = sl,nrow = sl)
  
  spcov<-matrix(NA,ncol = sl*tl,nrow=sl*tl)
  for(i in 1:tl)
  {
    for(j in i:tl)
    {
      tmp1<-myzeta[i,j]*my.matern(h=uniq.dist,
                                  c=myf_c[i,j],
                                  nu=mynus[i,j],
                                  sigma=sigma)
      temp[mat.index[,-3]]<-tmp1[mat.index[,3]]
      spcov[((i-1)*sl+1):((i-1)*sl+sl),((j-1)*sl+1):((j-1)*sl+sl)]<-spcov[((j-1)*sl+1):((j-1)*sl+sl),((i-1)*sl+1):((i-1)*sl+sl)]<-temp
    }
  }
  return(spcov)
}

library(doParallel)
registerDoParallel(cores=22)

sim_size = 20
mydelta=0
sigma= 0.98; delta = 0;  nu=0.5; a= 0.901; alpha =0.75; beta = 0.61
tl = 730
start.time =  Sys.time()
LR_null = foreach(s = 1:sim_size, .combine =cbind) %dopar% {
  spcov = st_covmat_func(c1=1/140,c2=1/140,nu=0.5, a=0.901, alpha=0.75, beta = 0.61, delta= 0, sigma = 1)
  mysim<-rmvnorm(n=1,mean = c(rep(0, times=sl*tl/2), rep(0, times=sl*tl/2)), sigma=spcov,method="chol")
  Z<-c(mysim)
  
  Lr = (LR_func(Z, tl))
  return(Lr)
}
end.time =  Sys.time()
end.time-start.time

th_wind = quantile(apply(LR_null, 2, max), 0.95)
th_wind = 15.2797


##################    detecting multiple changes  #########################
bin.seg.func = function(data_vec, Z_alpha, sl, tl, LR_func){
  
  n = tl
  
  Cpts = vector(); Seg_list = list(c(1,n))
  
  while (length(Seg_list)!=0) {
    
    Seg = Seg_list[[1]]  
    
    s = Seg[1]; t = Seg[2]
    data_seg = data_vec[((s-1)*sl+1):(t*sl)]
    lr_ts = LR_func(data_seg, tl = t-s+1)
    if(s==1 & t==n) {lr = lr_ts}
    
    if(max(lr_ts) >= Z_alpha){
      Seg_list = Seg_list[-1]
      tau_hat = which(lr_ts==max(lr_ts))
      r = tau_hat +s -1
      Cpts = c(Cpts, r)
      if(r!=s) {Seg_list = append(Seg_list, list(c(s,r))) }
      if(r!=(t-1)) {Seg_list = append(Seg_list, list(c(r+1,t))) }
    } else {Seg_list = Seg_list[-1] }
    #  }
  }
  
  return(list(Cpts, lr))
}

th_wind=40
start.time = Sys.time()
result = bin.seg.func(c(velocity), th_wind,sl, tl, LR_func)
end.time = Sys.time()
end.time-start.time


result[[1]]

plot(result[[2]])
which.max(result[[2]])








## plotting the pattern change
####### Transforming Longitude Latitude to Universal Transverse Mercator (UTM) projections ####
library(rgdal)
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## assigning a CRS (Coordinate Reference System)
  
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep=''))) ## Transforming the CRS
  return(data.frame(x= res$X, y= res$Y))
  #return(as.data.frame(res))
}

no_change = apply(windsqrt, 2, mean, na.rm=T)
pre_change = apply(windsqrt[1:572,], 2, mean, na.rm=T)
post_change = apply(windsqrt[573:731,], 2, mean, na.rm =T)
post_change-pre_change

change_data = data.frame(x=stations$easting, y=stations$northing, pre_change, post_change, change = post_change-pre_change, no_change)
#coordinates(change_data) = ~x+y

####Define the grid for interpolation
##Define the grid extent:
library(maps)
ireland = map_data("world", region = "ireland")

x.range <- range(ireland$long)
y.range <- range(ireland$lat)
grd <- as.matrix(expand.grid(x=seq(from=x.range[1], to=x.range[2], length.out = 100),
                             y=seq(from=y.range[1], to=y.range[2], length.out = 100) ))



ireland_latlon_all=NULL
for (i in unique(ireland$group)) {
  ireland_group = subset(ireland, group==i)
  ireland_latlon = cbind(ireland_group$long, ireland_group$lat)
  ireland_latlon_all = rbind(ireland_latlon_all, ireland_latlon, c(NA,NA))
}

library(mgcv)
pred_grd = grd[in.out(ireland_latlon_all,grd),] 
dim(pred_grd)
plot(pred_grd)

LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## assigning a CRS (Coordinate Reference System)
  
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep=''))) ## Transforming the CRS
  return(data.frame(x= res$X, y= res$Y))
  #return(as.data.frame(res))
}


pred_grd_km <- LongLatToUTM(pred_grd[,1],pred_grd[,2],zone = 29)/1000
stations_km = LongLatToUTM(stations$longitude, stations$latitude,zone = 29)/1000

## ordinary kriging with exponential covariance

library(mvtnorm)
C22 = my.matern(spDists(as.matrix(stations_km)), c= 0.01047690,sigma = 0.98, nu =0.5 )

crosslocs<-spDists(as.matrix(pred_grd_km),as.matrix(stations_km))
C12 = my.matern(h= crosslocs, c= 0.01047690,sigma = 0.98, nu =0.5 )

mu_hat = 0.00167659
mu_allloc = rep(mu_hat, 11)

krig.wind = rep(mu_hat, times = dim(pred_grd_km)[1]) + C12%*%solve(C22)%*%(change_data$no_change-mu_allloc )
krig.pre = rep(mu_hat, times = dim(pred_grd_km)[1]) + C12%*%solve(C22)%*%(change_data$pre_change-mu_allloc )
krig.post = rep(mu_hat, times = dim(pred_grd_km)[1]) + C12%*%solve(C22)%*%(change_data$post_change-mu_allloc )
krig.change = rep(mu_hat, times = dim(pred_grd_km)[1]) + C12%*%solve(C22)%*%(change_data$change-mu_allloc )
pred_results = data.frame(pred_grd,krig.wind, krig.pre, krig.post, krig.change)

## Prediction map (average wind speeds)
p = ggmap(map_ire)
p + ggtitle(" Average Wind Speeds")+
  geom_point(data = pred_results, aes(x, y, colour = krig.wind),size=2, alpha=0.4, shape=15) +
  scale_color_viridis("Wind speed")+ labs(x= "Longitude (degrees)", y="Latitude (degrees)")


## Prediction map (average wind speeds (before change))
p + ggtitle(" Average Wind Speeds (before change)")+
  geom_point(data = pred_results, aes(x, y, colour = krig.pre),size=2, alpha=0.4, shape=15) +
  scale_color_viridis("Wind speed")+ labs(x= "Longitude (degrees)", y="Latitude (degrees)")

## Prediction map (average wind speeds (before change))
p + ggtitle(" Average Wind Speeds (after change)")+
  geom_point(data = pred_results, aes(x, y, colour = krig.post),size=2, alpha=0.4, shape=15) +
  scale_color_viridis("Wind speed")+ labs(x= "Longitude (degrees)", y="Latitude (degrees)")

## Prediction map (average wind speeds ( change))
p + ggtitle("Change in wind speed pattern")+
  geom_point(data = pred_results, aes(x, y, colour = krig.change),size=2, alpha=0.4, shape=15) +
  scale_color_viridis("change")+ labs(x= "Longitude (degrees)", y="Latitude (degrees)") +
  geom_point(data=stations, aes(x=longitude,y=latitude), pch=16, size=2) +  
 geom_text(data=stations,aes(longitude,y=latitude,label=name), hjust=0, vjust=0)

pred_results$krig.change
stations$station.name
summary(krig.change)
########## prediction results #######################
Z = c(velocity)
miss_Z = which(is.na(Z))


set.seed(197)
test_index = sort(sample(1:length(Z), floor(0.20*length(Z)))) ##  random selection
test_index = (sl*571):(sl*573) ## around changepoint
test_index = c(((366:tl)-1)*sl+3, ((366:tl)-1)*sl+13, ((366:tl)-1)*sl+19)  ## removing the values from 3rd station (Shannon airport), 13 (athernry), 19 mt dillon



train_Z = Z[-test_index]
test_Z = Z[test_index]

miss_train_Z = which(is.na(train_Z))
train_Z = train_Z[-miss_train_Z]

miss_test_Z = which(is.na(test_Z))
#test_Z = test_Z[-miss_test_Z]

############# no changepoint  ############
mysigma= 0.98; mydelta = 0;  mya= 0.901; myalpha =0.75; mybeta = 0.61
Cz = st_covmat_func(c1=0.01047690, c2=0.01047690, nu=0.5, a=mya, alpha=myalpha, beta = mybeta, delta= mydelta, sigma = mysigma, t0=571)
Cz = Cz[-miss_Z, -miss_Z]
loc.all = cbind(rep(stations_km$x, tl), rep(stations_km$y, tl))
test_loc = loc.all[test_index,]

C22 = Cz[-test_index, -test_index]
C12 = Cz[ test_index, -test_index]
C11 = Cz[test_index, test_index]

dim(C22); dim(C12)

mu_hat = 0.00167659
mu_all = rep(mu_hat, length = sl*tl)
mu_all = mu_all[-miss_Z]
mu_all = mu_all[-test_index]

chol_st_nochange<-chol(C22)
inv_st_nochange = chol2inv(chol_st_nochange)
rm(chol_st_nochange)

krig.test.nochange = rep(mu_hat, times = dim(test_loc)[1]) + C12%*%inv_st_nochange%*%(train_Z-mu_all )
krig.var.nochange = diag(C11 - C12%*%inv_st_nochange%*%t(C12))

rm(Cz); rm(C22); rm(inv_st_nochange)
pred.nochange = data.frame(test_Z, krig.test.nochange, krig.var.nochange)
## rmse
sqrt(mean((pred.nochange$krig.test - pred.nochange$test_Z)^2, na.rm=T))
## crps 
library(scoringRules)
crps_test_nochange = apply(as.matrix(pred.nochange), 1, function(x) crps(x[1], family = "norm", mean=x[2], sd = sqrt(x[3])))
mean(crps_test_nochange)

## logS
logs_test_nochange = apply(as.matrix(pred.nochange), 1, function(x) logs(x[1], family = "norm", mean=x[2], sd = sqrt(x[3])))
mean(logs_test_nochange)



#############  changepoint  ############
mysigma= 0.98; mydelta = 0;  mya= 0.901; myalpha =0.75; mybeta = 0.61
Cz = st_covmat_func(c1=0.010862112, c2=0.008968451, nu=0.5, a=mya, alpha=myalpha, beta = mybeta, delta= mydelta, sigma = mysigma, t0=571)
Cz = Cz[-miss_Z, -miss_Z]

C22 = Cz[-test_index, -test_index]
C12 = Cz[ test_index, -test_index]
C11 = Cz[test_index, test_index]

dim(C22); dim(C12)

mu1 = 0.058838115; mu2 = -0.196997978
mu_all = c(rep(mu1, 572*sl), rep(mu2, (tl-572)*sl))
mu_all = mu_all[-miss_Z]
mu_all = mu_all[-test_index]

chol_st_change<-chol(C22)
inv_st_change = chol2inv(chol_st_change)
rm(chol_st_change)

mu_test = ifelse(test_index < sl*572,mu1, mu2 )

krig.test = mu_test + C12%*%inv_st_change%*%(train_Z-mu_all )
krig.var = diag(C11 - C12%*%inv_st_change%*%t(C12))

rm(Cz); rm(C22); rm(inv_st_change)
pred.change = data.frame(test_Z, krig.test, krig.var)

## rmse
sqrt(mean((pred.change$krig.test - pred.change$test_Z)^2))

## crps 
library(scoringRules)
crps_test = apply(as.matrix(pred.change), 1, function(x) crps(x[1], family = "norm", mean=x[2], sd = sqrt(x[3])))
mean(crps_test)

## logS
logs_test = apply(as.matrix(pred.change), 1, function(x) logs(x[1], family = "norm", mean=x[2], sd = sqrt(x[3])))
mean(logs_test)



which.max(crps_test_nochange-crps_test)
mat = matrix(1:(sl*tl), nrow = tl, ncol=sl, byrow=T)
which(mat==test_index[2538], arr.ind = T)
length(test_index)




############## Implementing naive optimistic search  #######################
Z = c(velocity)
optimfull = optim(par=c(0,1/100),logL_markov, Z= Z ,sl=sl,tl=tl, k=Tn, hypo=0)
loglikfull = -optimfull$value

LR = function(k, data_vec) {
  
  optimalt = optim(par=c(0,0,1/100, 1/100),logL_markov,Z= data_vec ,sl=sl,tl=tl, k=k, hypo=1)
  loglikalt = -optimalt$value
  
  lr_ts = 2*(loglikalt-loglikfull) 
  
  return(lr_ts)
}

count=1
nOS = function(l, s, r, data_vec){
  
  lr_vec = vector()
  if((r-l) <= 5){
    for (i in (l+1):(r-1)) { lr_vec[i] = LR(i, data_vec); count = count+1 }
    s_hat =  which.max(lr_vec)
    return(c(s_hat, max(lr_vec, na.rm = T), count))
  }
  if((r-s) > (s-l)){
    count =count +1
    w = ceiling(r - (r-s)*step_size)
    if(LR(w, data_vec) >= LR(s, data_vec)){
      nOS(s,w, r, data_vec)} else {nOS(l,s, w, data_vec)}
  } else {
    count = count+1
    w = floor(l+(s-l)*step_size)
    if(LR(w, data_vec) >= LR(s, data_vec)){
      nOS(l,w, s, data_vec)} else {nOS(w,s, r, data_vec)}
  }
}



start.time = Sys.time()
L=1; R=730; step_size = 0.3
result = nOS(l= L, s = floor((L+ step_size*R)/ (1+step_size)),  r=R, data_vec = c(velocity))
end.time = Sys.time()
end.time- start.time



### optimistic binary segmentation

opt.bin.seg.func = function(data_vec, Z_alpha, sl, tl){
  
  n = tl
  
  Cpts = vector(); Seg_list = list(c(1,n))
  
  while (length(Seg_list)!=0) {
    
    Seg = Seg_list[[1]]  
    
    
    s = Seg[1]; t = Seg[2]
    data_seg = data_vec[((s-1)*sl+1):(t*sl)]
    
    R = t-s+1
    result  = nOS(l= L, s = floor((L+ step_size*R)/ (1+step_size)),  r=R, data_seg)
    #lr_ts = LR_func(data_seg, tl = t-s+1)
    
    if(result[2] >= Z_alpha){
      Seg_list = Seg_list[-1]
      tau_hat = result[1]
      r = tau_hat +s -1
      Cpts = c(Cpts, r)
      if(r!=s) {Seg_list = append(Seg_list, list(c(s,r))) }
      if(r!=(t-1)) {Seg_list = append(Seg_list, list(c(r+1,t))) }
    } else {Seg_list = Seg_list[-1] }
  }
  
  return(Cpts)
}

































### modeling and removing seasonality
Sn = sl = dim(stations)[1];  Tn = tl = dim(wind.data)[1]

wind_vec = c(windsqrt)


t_cov = rep(1:Tn, each= Sn)
cos_cov = cos(2*pi*t_cov/365)
sin_cov = sin(2*pi*t_cov/365)

lm.fit = lm(wind_vec~ cos_cov + sin_cov)
summary(lm.fit)

seasonality = lm.fit$coefficients[2]*cos_cov+lm.fit$coefficients[3]*sin_cov
seasonality_mat = matrix(seasonality, nrow = Tn, ncol = Sn, byrow=T)
head(seasonality_mat)

plot(seasonality_mat[,2], type="l")



