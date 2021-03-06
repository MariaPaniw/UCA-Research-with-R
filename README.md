# UCA Research with R

This is a repository of researchers at the UCA to present their work with R code. 

- [ARIMA in fisheries](#arima-in-fisheries)
- [Lasso regression](#lasso-regression)
- [IPMs and 3D graphs](#ipms-and-3d-graphs)

## ARIMA in Fisheries

**Who made this?**

I am Remedios Cabrera Castro, a faculty member of the Biology Department. Our research group specializes in population dynamics of fish stocks. Currently, we are working with fisheries data on the Pacific [halibut](https://en.wikipedia.org/wiki/Pacific_halibut) (*Hippoglossus stenolepis*, Schmidt, 1904). Specifically, we are working with [Autoregressive integrated moving average (ARIMA) models](http://people.duke.edu/~rnau/411arim.htm).


**What does this code do?**
The code below demonstrates the use of ARIMA models for fisheries using a "catch per unit effort" (CPUE) time series. This code, along with much more information on its implementation (in Spanish) is available at: [http://rodin.uca.es/xmlui/handle/10498/17877](http://rodin.uca.es/xmlui/handle/10498/17877) 


Example
========================================================


```r
#Carga de librería. 
library(tseries) # install package if necessary 

setwd("G:/Teaching/IntroR/")

#Load data (average CPUE)

X = read.csv("MM31998-2010.csv") # you can find it in the figures folder
```


```r
#Representación gráfica. 

#Definición de la serie temporal con datos propios. 

X<-ts(X,start=1998,frequency=13) 
plot(X, type="o",col ="red") 
```

![ARIMA1](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/ARIMA_1.png)

```r
#Gráficos de correlación. 
acf(X,lag.max=42) 
```

![ARIMA2](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/ARIMA_2.png)

```r
pacf(X,lag.max=42) 
```

![ARIMA3](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/ARIMA_3.png)

```r
Box.test(X, lag = 6, type = c("Box-Pierce", "Ljung-Box")) 
```

```

	Box-Pierce test

data:  X
X-squared = 419.5567, df = 6, p-value < 2.2e-16
```


```r
library(lmtest) 
n<-length(X)  
ti<-1:n  
dwtest(X~ti) 
```

```

	Durbin-Watson test

data:  X ~ ti
DW = 0.6392, p-value < 2.2e-16
alternative hypothesis: true autocorrelation is greater than 0
```

```r
#Elección del modelo, aditivo o multiplicativo. 
dif<-diff(X) 
sd(dif)/abs(mean(dif)) 
```

```
[1] 48.44421
```

```r
inc<-X[-1]/X[-length(X)] 
sd(inc)/abs(mean(inc)) 
```

```
[1] 0.03341487
```


```r
#Descomposición de la serie por componentes. 
plot(decompose(X,type="multiplicative")) 
```

![ARIMA4](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/ARIMA_4.png)

```r
adf.test(X) 
```

```

	Augmented Dickey-Fuller Test

data:  X
Dickey-Fuller = -5.3119, Lag order = 5, p-value = 0.01
alternative hypothesis: stationary
```

```r
kpss.test(X,null=c("Trend")) 
```

```

	KPSS Test for Trend Stationarity

data:  X
KPSS Trend = 0.0479, Truncation lag parameter = 3, p-value = 0.1
```

```r
X_te<-diff(X, lag=1) 
```


```r
#No ciclo estacional. 
X_te_es<-diff(X_te, lag=13) 
plot(decompose(X_te_es,type="multiplicative")) 
```

![ARIMA5](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/ARIMA_5.png)

```r
plot(X_te_es) 
```

![ARIMA6](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/ARIMA_6.png)

```r
#Correlaciones de la serie estacionaria. 
layout(matrix(c(1,2),1,2)) 
acf(X_te_es,lag.max=42) 
pacf(X_te_es,lag.max=42) 
```

![ARIMA7](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/ARIMA_7.png)

```r
#Realización del modelo. 
fit_1<-arima(X, order =c(0,0,1), seasonal=list(order=c(0,0,1), period=13)) 
X.pred<-predict(fit_1,n.ahead=65) 

#Comprobación de residuos. 
str(fit_1) 
```

```
List of 13
 $ coef     : Named num [1:3] 1 0.303 471.509
  ..- attr(*, "names")= chr [1:3] "ma1" "sma1" "intercept"
 $ sigma2   : num 271
 $ var.coef : num [1:3, 1:3] 6.07e-04 -1.75e-08 7.56e-08 -1.75e-08 4.48e-03 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:3] "ma1" "sma1" "intercept"
  .. ..$ : chr [1:3] "ma1" "sma1" "intercept"
 $ mask     : logi [1:3] TRUE TRUE TRUE
 $ loglik   : num -717
 $ aic      : num 1442
 $ arma     : int [1:7] 0 1 0 1 13 0 0
 $ residuals: Time-Series [1:169] from 1998 to 2011: -0.876 15.458 2.027 22.758 28.34 ...
 $ call     : language arima(x = X, order = c(0, 0, 1), seasonal = list(order = c(0, 0, 1),      period = 13))
 $ series   : chr "X"
 $ code     : int 0
 $ n.cond   : int 0
 $ model    :List of 10
  ..$ phi  : num(0) 
  ..$ theta: num [1:14] 1 0 0 0 0 ...
  ..$ Delta: num(0) 
  ..$ Z    : num [1:15] 1 0 0 0 0 0 0 0 0 0 ...
  ..$ a    : num [1:15] -55.94 -31.05 -5.45 -6.54 -6.46 ...
  ..$ P    : num [1:15, 1:15] 0.00 1.11e-16 0.00 0.00 0.00 ...
  ..$ T    : num [1:15, 1:15] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ V    : num [1:15, 1:15] 1 1 0 0 0 ...
  ..$ h    : num 0
  ..$ Pn   : num [1:15, 1:15] 1.01 1.00 -5.20e-17 -5.23e-17 5.27e-17 ...
 - attr(*, "class")= chr "Arima"
```

```r
layout(matrix(c(1,2),1,2)) 
acf(fit_1$residuals) 
pacf(fit_1$residuals) 
```

![ARIMA8](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/ARIMA_8.png)

```r
Box.test(fit_1$residuals, lag=7) 
```

```

	Box-Pierce test

data:  fit_1$residuals
X-squared = 184.0872, df = 7, p-value < 2.2e-16
```


## IPMs and 3D graphs

**Who made this?**

I am Maria Paniw, a PhD student at the Biology Department, attached to the project [BREATHAL](https://sites.google.com/site/ucaherrizalab/), headed by Fernando Ojeda.

In a nutshell, I have focused my research on stochastic population dynamics of disturbance-adapted plant species. In concrete, I work on the population dynamics of *Drosophyllum lusitanicum*, a fire-adapted carnivorous plant species.

![Drosophyllum](http://www.carnivorousplants.org/howto/GrowingGuides/Images/Drosophyllum.jpg)


**What does this code do?**
The code below will show you how to create several types of functions in R and run loops. You'll also see how to construct integral projection models to describe population dynamics and plot 3D data via heat maps.

**What exactly do I mean by 'population dynamics'?**

All populations consist of individuals that have specific life cycles from birth to death. The transitions of individuals between life-cycle components is governed by vital rates such as survival and growth. I use models to link vital-rate transitions of individuals to population-level processes. 

Here is the life cycle of *Drosophyllum*:

![Drosophyllum life cycle](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/DrosoLife.png)

If the tranistions associated with different stages are discrete, they can be translated into a matrix: 

![Matrix](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/lifeCycleMat.png)

The big assumption here is that all individuals grouped within a stage behave exactly the same. If this assumption doesn't hold and individuals may be best described by continuous transitions, you can use integral projection models to link continuous vital-rate transitions to population-level metrics. 

The code below shows you just how to do that. Much of it is based on [Ellner & Rees, Am. Nat., 2006, 410-428](http://www.jstor.org/stable/10.1086/499438)


```r
# load some necessary functions (intall if necessary)

library(plyr)
library(fields)
library(MASS)
library(Cairo)
```

# Parameters for vital rates

The first thing to do is to get the parameters for vital-rate functions.

I obtained the parameters from Bayesian models. 

```r
params=data.frame(a0.surv=-3.3,bc.surv=0.72,a1.surv.one=1.9,a1.surv.two=1.15,
                  a1.surv.three=-1.86,a1.surv.four=-1.2,
                  a0.gr=2.31,bc.gr=0.67,a1.gr.one=1.3,a1.gr.two=-0.31,
                  a1.gr.three=-1.01,a1.gr.four=0.01,
                  bcTSF.gr.one=-0.11,bcTSF.gr.two=-0.008,bcTSF.gr.three=0.16,
                  bcTSF.gr.four=-0.04,a0.sds=3.95,a1.sds.one=-0.68,
                  a1.sds.two=0.41,a1.sds.three=0.17,a1.sds.four=0.1,
                  a0.fl=-9.2,bc.fl=1.59,a1.fl.two=-4.79,
                  a1.fl.three=3.98,a1.fl.four=0.81,
                  bcTSF.fl.two=0.55,bcTSF.fl.three=-0.59,bcTSF.fl.four=0.04,
                  a0.fs=-3.62,bc.fs=0.66,a1.fs.two=-0.47,
                  a1.fs.three=0.33,a1.fs.four=0.14,
                  a0.fps=-0.07,bc.fps=0.24,a1.fps.two=0.2,
                  a1.fps.three=-0.04,a1.fps.four=-0.17,
                  a0.goCont=-2.55,a1.goCont.burned=0.66,a1.goCont.unburned=-0.66,
                  a0.staySB=1.06,a1.staySB.burned=-0.67,a1.staySB.unburned=0.67,
                  a0.outSB=1.06,a1.outSB.burned=-0.67,a1.outSB.unburned=0.67)

# replace some characters in the column names
colnames(params) = gsub("four",">three",colnames(params))
attach(params)
```


# Creating vital-rate functions: above-ground (continuous)

The second thing to do is to create functions that decribe the vital rates of  *Drosophyllum*. All vital rates that decribe the fates of above-ground individuals are functions of size *z* and covariate *cov* - time since fire (TSF).


```r
# SURVIVAL:

S.fun <- function(z,cov) {
 
    mu.surv=exp(a0.surv+get(paste("a1.surv.",cov,sep=""))+bc.surv*z)
  
    return(mu.surv/(1+mu.surv))
}
```



```r
# GROWTH 

GR.fun <- function(z,zz,cov){

    growth.mu=(a0.gr+get(paste("a1.gr.",cov,sep=""))+
                 (bc.gr+get(paste("bcTSF.gr.",cov,sep="")))*z)
    

    var.res=0.3986543 # we assume constant residual variance 
    # Density distribution function of the normal distribution
    gr1 = sqrt(2*pi*var.res)
    gr2 = ((zz-growth.mu)^2)/(2*var.res)
  
  return(exp(-gr2)/gr1)
}
```



```r
## SEEDLING SIZES (same approach as in growth function)

SDS.fun <- function(z,zz,cov){
  
  sds.mu=(a0.sds+get(paste("a1.sds.",cov,sep="")))
  

  var.res=0.3329295

  # Density distribution function of the normal distribution
  sds1 = sqrt(2*pi*var.res)
  sds2 = ((zz-sds.mu)^2)/(2*var.res)
  
  return(exp(-sds2)/sds1)
  
}
```



```r
# PROBABILITY OF FLOWERING 

FL.fun <- function(z,cov) {
  
  mu.fl=exp(a0.fl+get(paste("a1.fl.",cov,sep=""))+
              (bc.fl+get(paste("bcTSF.fl.",cov,sep="")))*z)
 
  return(mu.fl/(1+mu.fl))
}

# NUMBER OF FLOWERING STALKS 

FS.fun <- function(z,cov) {
  
  mu.fs=exp(a0.fs+get(paste("a1.fs.",cov,sep=""))+bc.fs*z)
  
  return(mu.fs)
}
```



```r
# NUMBER OF FLOWERS PER STALK

FPS.fun <- function(z,cov) {
  
  mu.fps=exp(a0.fps+get(paste("a1.fps.",cov,sep=""))+bc.fps*z)
  
  return(mu.fps)
}
```


# Creating vital-rate functions: below-ground (discrete) 

Now, we describe the fates of seeds, which can either:

- germinate the growing season following maturation (*goCont*)
- go into and survive in the seed bank (*staySB*)
- germinate out of the seed bank (*outSB*)

All vital rates that decribe the fates of seeds are functions of post-fire habitat status: *burned* vs. *unburned*.


```r
# IMMEDIATE GERMINATION (goCont):

goCont.fun <- function(pfs) {

      mu.goCont = exp(a0.goCont+get(paste("a1.goCont.",pfs,sep="")))


  return(mu.goCont/(1+mu.goCont))
}
```


```r
# PERCENTAGE STAYING IN THE SEED BANK (staySB):

staySB.fun <- function(pfs) {

      mu.staySB <- exp(a0.staySB+get(paste("a1.staySB.",pfs,sep="")))
  
  return(mu.staySB/(1+mu.staySB))
}


# PERCENTAGE GERMINATING OUT OF THE SEED BANK (outSB):

outSB.fun <- function(pfs) {

      mu.outSB <- exp(a0.outSB+get(paste("a1.outSB.",pfs,sep="")))

  return(mu.outSB/(1+mu.outSB))
}
```

Some of the vital rates may need additional corrections, often because they were overestimated. When these corrections are constants and not models (functions), you can just add them directly.


```r
# data frame with the correction factors for each TSF
# since we are building IPMs for each TSF

corr=data.frame(fGerm=c(NA,NA,0.06,0.06,0.03), #for goCont
                mort=c(0.16,NA,0.2,0.2,0.2), # for seed mortality
                TSF=c("zero","one","two","three",">three"))
```

# Function to discretize the integral of the projection model



```r
#Discretization occurs via minimum and maximum sizes
minsize=0 # minimum size
maxsize=9.6 # maximum size

IPMkernel<-function(n) { # n defines the size of the kernel  

  # the midpoints for the integration
  b <- minsize+c(0:n)*(maxsize-minsize)/n # interval that each cell of the matrix covers 
  h <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint

  # The zero (right after fire) matrix
  if(cov=="zero"){
    outSB <- 0.81*(1-mort) 
    staySB <- 0.1 # constant assumed
    goCont=0.01 # so that we don´t devide by zero! (see below)
    goSB=0
    goSoil <- 1 # to avoid deviding by 0
    # survival (S), growth (G), and fecundity (FecALL) are all zero
    S <-matrix(0,n,n) 
    G <- matrix(0,n,n)
    FecALL=matrix(0,n,n)
    # the relevant non-0 transition is the size of seedlings
    R <- (t(outer(h,h,SDS.fun,cov="one")))
    
    # TSF 1
    }else if(cov=="one"){
      
    outSB <- outSB.fun(pfs="burned")
    staySB <-0.05
    goCont=0.01 # so that we don´t devide by zero! (see below)
    goSB=0
    goSoil <- 1
    S <- diag(S.fun(h,cov)) # Survival Matrix 
    G <- t(outer(h,h,GR.fun,cov)) # Growth Matrix
    #Recruits distribution
    R <- (t(outer(h,h,SDS.fun,cov)))
    FecALL=matrix(0,n,n)# no seeds produced, therefore 0 fecundity
    
    # scale G and R below so columns sum to 1
    G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    
  # TSF 2 up   
  }else{
    outSB=outSB.fun(pfs="unburned")
    staySB=staySB.fun(pfs)
    goCont=goCont.fun(pfs=pfs.goCont)*fGerm
    goSB=1-(goCont/fGerm)
    goSoil = 1-mort
    S <- diag(S.fun(h,cov)) # Survival Matrix 
    G <- t(outer(h,h,GR.fun,cov)) # Growth Matrix
    
    #Recruits distribution
    R <- (t(outer(h,h,SDS.fun,cov)))
  
    #Probability of flowering
    Fec01 = (diag(FL.fun(h,cov)))
    
    #Number of flowering stalks 
    Fec02 = (diag(FS.fun(h,cov)))
    
    #Number of flowers per stalk
    
    Fec03= (diag(FPS.fun(h,cov)))
    
    #Number of seeds per flower that survive to become offspring 
    Fec04 = (diag(rep(9.8,n)))
    
    
    FecALL= Fec01*Fec02*Fec03*Fec04*goCont*goSoil 
    
    # scale D and G so columns sum to 1
    G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        
  }
 # scale R so columns sum to 1
  R <- R/matrix(as.vector(apply(R,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  
  return(list(S=S,G=G,FecALL=FecALL,R=R,meshpts=h,goCont=goCont,goSoil=goSoil,goSB=goSB,outSB=outSB,staySB=staySB,fGerm=fGerm))
  
}
```

Loop through all the TSF to create IPMs for each
================


```r
TSF.name=c("zero","one","two","three",">three")

#  To link pfs  to TSF
pfs.name=c("NA","NA","burned","unburned","unburned")
pfs.name.goCont=c("NA","NA","unburned","unburned","unburned")

bins=50 # how large is the kernel? 
discr=1 # hpw many discrete stages
matP=array(0,c(bins+discr,bins+discr,5))# save S + G matrices
matF=array(0,c(bins+discr,bins+discr,5)) # save fecundity matrices
matPF=array(0,c(bins+discr,bins+discr,5))

# I have 5 TSF:
for(t in 1:5){
  cov=TSF.name[t]
  fGerm=corr[corr$TSF==cov,"fGerm"]
  mort=corr[corr$TSF==cov,"mort"]
  pfs=pfs.name[t]
  pfs.goCont=pfs.name.goCont[t]
     
  M=IPMkernel(bins)
  
  if(t==5) M$outSB <- M$outSB * 0.17 # drastically descrease germination from the seed bank in TSF >3 to approximate germination observed in the field
  
  ### Create P and F matrices including discrete stages
    
  Pmat.cont <- M$G%*%M$S
  Pmat.discr = c(M$staySB,M$outSB*M$R[,2])
  Pmat = cbind(Pmat.discr,rbind(rep(0,length(M$meshpts)),Pmat.cont))
  
  matP[,,t] = Pmat
    
  Fmat.cont <- M$R%*%M$FecALL
  Fmat.discr =rep(0,length(M$meshpts)+1)
  Fmat=cbind(Fmat.discr,rbind(diag(M$FecALL)*M$goSB/(M$goCont),
                              Fmat.cont))
  matF[,,t] = Fmat
    
  mat <-Pmat+Fmat
  matPF[,,t] = mat
    
}
```

# Plot the resulting transitions



```r
# Plot the P kernels
maxSize=9.6
y=M$meshpts


fire=c("0","1","2","3",">3")

set.panel()

par(mfrow=c(1,5),mar=c(1.3,1.3,1,1),oma=c(3,3.5,2,4)) 

set.panel( 1,5) # 1X5 matrix of plots

# now draw all your plots using image command
count=0

  for(f in 1:length(fire)){
    count=count+1

    #Change values to make plots without the white colors
    conv=function(x){ifelse(x>0.1,x<-0.1,ifelse(-0.0006>x,x<-0,x))}    
    V=apply(matP[,,f],c(1,2),conv)
 
    if(count==1){
      image(c(-1,y),c(-0.67,y),t(V), xlab=expression(paste("size ", italic(t))),ylab=expression(paste("size ", italic(t), " + 1")),cex.axis=2,cex.lab=2.7,zlim=c(0,0.1),col=tim.colors(900), xpd="n")
      
      abline(h=-0.25,col="white",lwd=1.8)
      abline(v=-0.43,col="white",lwd=1.8)
      
      title( main=paste("TSF: ",fire[f]," years", sep=""), line =1.1,cex.main=3,xpd="n")
    } else {
      
      image(c(-1,y),c(-0.67,y),t(V), xlab=expression(paste("size ", italic(t))),ylab="",yaxt="n",cex.axis=2,cex.lab=2.7,zlim=c(0,0.1),col=tim.colors(900),xpd="n")
      abline(h=-0.25,col="white",lwd=1.8)
      abline(v=-0.43,col="white",lwd=1.8)
      title( main=paste(fire[f]," years", sep=""), line =1.1,cex.main=3,xpd="n")
    }
    
  }

par(oma=c(0,0,0,2.5))# reset margin to be much smaller.

image.plot( legend.only=TRUE,zlim=c(0,0.1),col=tim.colors(900),legend.size=1.5,xpd="n") 

set.panel()
```

# Growth and Survival

![P kernels](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/Pmat.png)

You can use the same code above, replace *matP* by *matF*, play around with the `conv()` function (upper limits should be > 1), and get the transitions for fecundity.

# Fecundity

![F kernels](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/Fmat.png)


# Stochastic simulations of population growth

You can use the matPF array, which holds IPMs for each time-since-fire (TSF) state, to simulate stochastic population dynamics where the stochastic part comes from fires, i.e., a Markov chain environment in which the states are TSF and transitions are based on probability of fire! 


```r
### Add simulations 

# Simulation function (just one of many possibilities)
# This function takes a environmental transition matrix (trans) and the number of simulation (N years)
# and returns a character vector of the sequence of environments in N years 

simula <- function(trans,N) {
  # This function tells you sample the rownames (i.e., environment at t+1), each having a probability = column name that you got at previous iteration(environment at t) 
  state.at.N <- function(char,trans) {
    sample(rownames(trans),1,prob=trans[,char])
  }
  
  sim <- character(N)
  sim[1] <- "1" # start with 0 matrix and populate 1st character with "1"
  
  # Populate characters 2-N with the row names samples obtained from (state.at.N)
  for (b in 2:N) {
    sim[b] <- state.at.N(sim[b-1],trans)
  }
  
  sim
}
### Prepare parameters


num.of.sim = 10 # number of simulations to run (here few to speed process up);
# Define frequency of fire (= 1/fire return interval) 

freqarray=c(1/10,1/30,1/50,1/100) # in the ms, we used 16 different frequencies 

# Define simulation time
tr=500 #the discard time
ts=2000 

# run each simulation for ts+tr years  
trun=ts+tr

# The different fire environments 
env=c("zero","one","two","three",">three")
env.n=length(env)


lambda.s=array(0,c(length(freqarray),num.of.sim)) # save final stochstic lambda

# Set some values needed to for further calculations 
    
    n2=(dim(matPF)[1])^2
    n=dim(matPF)[1]
    
     for(si in 1:num.of.sim){ # loop over simulations
      
      for (f in 1:length(freqarray)) {
        #create environmental transition matrices
        
        fire=freqarray[f]; 
        
        P=matrix(c(fire,1-fire,rep(0,3),
                   fire,0,1-fire,rep(0,2),
                   fire,0,0,1-fire,0,
                   fire,0,0,0,1-fire,
                   fire,0,0,0,1-fire),nrow=5,ncol=5,byrow=F)
        colnames(P) <- c("1","2","3", "4","5")
        row.names(P) <- c("1","2","3", "4","5")
        
        
        ## Simulations of TSF states
        simulate=as.numeric(simula(P,trun+1))

        # Calculate stochastic lambda
        ########################################################
        
        states <- simulate
        growth <- array(0,ts)   
        
        # Initialize population vectors (n0)
        vec1=c(1000,rep(1,n-1)) # start with 1000 seeds 
        vec1 <- vec1/sum(vec1)
        vec1 <- t(vec1) 
        
        # ITERATION TO CALCULATE LAMBDA FOR EACH TIME STEP
        for (i  in 1:trun){
          i2 <- states[i]
          mat1 <-  matPF[,,i2]
          vec1 <- mat1%*%as.numeric(vec1)
          growth1 <- sum(vec1) # population growth at one time step 
          vec1 <- vec1/growth1
          if( i > tr){ # after the burn-in, save the growth rate for each time step      
            i1 <- i - tr
            growth[i1] <- growth1
          }
          
        }
        
        a1=sum(log(growth[1:ts]))
        
        lambda.s[f,si]= a1  
        
      }
      
      
    }
    
lambda.s/ts # stochastic growth rate decreases with increasing fire return
```

```
            [,1]         [,2]         [,3]         [,4]         [,5]
[1,]  0.19403244  0.194449845  0.184284959  0.175270054  0.166475084
[2,]  0.03281634  0.041596917  0.069778985  0.046907929  0.092693754
[3,]  0.01897316  0.002842226  0.004510146  0.003573959 -0.002749331
[4,] -0.01381451 -0.032761771 -0.001564432 -0.017575502 -0.032322083
            [,6]        [,7]         [,8]         [,9]       [,10]
[1,]  0.14945419  0.17299825  0.153782749  0.139537006  0.15444790
[2,]  0.06209270  0.05104792  0.035423477  0.030231530  0.07482913
[3,]  0.01533621  0.02664775  0.012921860  0.009208009  0.03176619
[4,] -0.04393360 -0.02290194 -0.006742557 -0.023728928 -0.03458702
```


## Lasso regression

**Who made this?** 

I am Jorge del Rosario Fernandez Santos, a postdoc at the Department of Physical Education, UCA.

**What does this code do?**
The code below will show you how to fit a linear model using the Lasso regression which is a method that accomplishes both shrinkage and variable selection. For that purpose, we are going to use the [glmnet package](http://www.jstatsoft.org/article/view/v033i01) and the Bayesian version proposed by [Lykuo & Ntzoufras (2013)](http://link.springer.com/article/10.1007%2Fs11222-012-9316-x). 


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
**Variable selection**

The first step is to produce correlated normal distributed data. We are going to use the method proposed by [Hardin et al. (2013)](http://arxiv.org/abs/1106.5834) for generating a correlation matrix and the [Cholesky descomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition) to generate multivariate random correlated data.

```{r eval=FALSE}
# Generate the correlation matrix
set.seed(7)
noise.iden <- noisecor(diag(10), epsilon = .7, eidim=2)

# Use Cholesky descomposition
U = t(chol(noise.iden))
numvars = dim(U)[1]
numobs = 100
random.normal = matrix(rnorm(numvars*numobs,0,1), nrow=numvars, ncol=numobs);
X = U %*% random.normal
newX = t(X)
raw = as.data.frame(newX)
names(raw) = c('resp', 'pred1', 'pred2', 'pred3','pred4', 'pred5', 'pred6','pred7', 'pred8', 'pred9')

# Prepare the matrix Y with the response variable and the matrix X with 9 independent variables 
X <- as.matrix(raw[,2:10])
Y <- as.matrix(raw[,1])
```

### PACKAGE GLMNET 

As state in the [glmnet vignette](http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html), _glmnet is a package that fits a generalized linear model via penalized maximum likelihood_. To perform a Lasso regression we have to set the elastic-net penalty (Î±) = 1. The **glmnet** function is going to return different models for the researcher to choose from. 

```{r eval=FALSE}
library(glmnet)
fit <- glmnet(X,Y)
# And a summary of the glmnet path
print(fit)
# We can visualize the coefficients using the plot function
plot(fit, label = TRUE)
```
![fit](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/fit.png)

If we prefer that the software select one of them by cross-validation we have to use the **cv.glmnet** function. 
```{r eval=FALSE}
cvfit = cv.glmnet(X, Y)
# We can plot the cross-validation curve
plot(cvfit)
```
![cvfit](https://raw.githubusercontent.com/MariaPaniw/UCA-Research-with-R/master/Figures/cvfit.png)
```{r eval=FALSE}
# We can get the coefficients for the model with:
# the value of Î» that gives minimum mean cross-validated error
cvfit$lambda.min
# or the alue of Î» that gives the most regularized 
# model such that error is within one standard error of the minimum
cvfit$lambda.1se
# To get the coefficients
coef(cvfit, s = "lambda.min")
```

### BAYESIAN LASSO
[Lykuo & Ntzoufras (2013)](http://link.springer.com/article/10.1007%2Fs11222-012-9316-x) proposed a Bayesian implementation of the lasso regression focusing on the appropiate specification of the shrinkage parameter Î» through Bayes factors that evaluate the inclusion of each covariate in the model. The authors create a function _blvs_ which performs a Bayesian version of the Lasso by imposing the Double-Exponential prior distribution on the linear regression coefficients. For technical details please look at the reference.

```{r eval=FALSE}
# First we have to get the benchmark correlation for a  given sample size
# in our case n=100
rho = bf.rho(100, 3)
# Then the value of the shrinkage parameter
lambda = bf.lambda(rho, 100, 1)
# and finally the blvs function
bLasso = blvs(X, Y, lambda=lambda, alpha2=0.0001, gamma2=10000, nburn=1000, ndraw=20000)
```

