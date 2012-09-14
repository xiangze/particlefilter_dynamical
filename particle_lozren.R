source("lotz.f.noise.R")
source("rossler.noize.R")
source("resample.R")

dt <- 0.01
#dt <- 0.001


N <- 100
interval <- 100
lag <- 1

rv <- 0.0 #original noise var
xv <- 0.01  #system noise var
yv <- 0.03   #obesebation noise var
nv <- 0.001 # dispersion of PF

psize <- c(-1000,1000)
preshow <- T
show <- T

type <- "ros"
rparams <- c(4,5.7,14,18)
times <- 30000
#type <- "loz"
#rparams <- c(4,16,28,331)

#for(N in c(100,300)){
for(N in c(110)){
  for(ro in rparams){
    
    if(type=="loz"){
      lotz.params$ro <- ro
      lz <- function(x,v){x+dt*lotz.noise(0,x,v,lotz.params)}
      pname <-paste("loz","N",N,ro,"int",interval,".png",sep="_")
    }else{
      rossler.params$c <- ro
      lz <- function(x,v){x+dt*rossler.noise(0,x,v,rossler.params)}
      pname <-paste("ros","N",N,ro,"int",interval,".png",sep="_")
    }
    
#make original(hidden) time series (x) with system noise
#observed sequence :y
    y <- array(dim=c(times/interval,3))
    xo <- array(dim=c(1,3))
    xo[1,] <- c(0.2,0.2,0.1)

for (t in 2:times){  
  xo <- rbind(xo,lz(xo[t-1,],rv))
  if(t%%(interval)==0){
    y[t/interval,] <-xo[t,]+rnorm(3,0,yv);
  }
}

    xmin <- min(xo)
    xmax <- max(xo)  
    trange <- c(xmin,xmax)

    
if(preshow){
##   png(paste("org_phase",pname,sep="_"))
##   plot(xo[,2],xo[,3], type='l',xlim=pxsize,ylim=pysize, axes=F,xlab="",ylab="")
##   par(new=T)  
##   plot(y[,2],y[,3], type='l',xlim=pxsize,ylim=pysize,col="red",xlab="x",ylab="y")
##   dev.off()

  png(paste("org_time",pname,sep="_"))
  plot(seq(1,times),xo[,2], ylim=trange,type='l', axes=F,xlab="",ylab="")
  par(new=T)
  plot(seq(1,times,by=interval),y[,2],ylim=trange,type='b',col="red",xlab="x",ylab="y")
  dev.off()
  
}

#prediction by particles 
    x <- array(dim=c(times,N,3))
    x[1,,] <- matrix(runif(N,-0.2,0.2))

    for (t in 2:times){

      for(n in 1:N){
        x[t,n,] <- lz(x[t-1,n,],xv)
      }
    
  if(t%%(interval)==0 && t/interval>lag){
    weight <- rep(0,N) 
    for(n in 1:N){
      for(i in 0:(lag-1)){
        weight[n] <- weight[n]+dnorm(y[(t/interval-i),],x[(t-i*interval),n,],yv,log=T)
      }
      
    }
#    print(weight)
#    browser()

    wm <- max(weight)
    weight <- sapply(weight,function(x){exp(x-wm)})
#    print(wm)
    
    costs <- weight/sum(weight)*N
    x[t,,] <- resample(x[t,,],costs)
    x[t,,] <- x[t,,]+matrix(rnorm(length(x[t,,]),0,nv),nc=ncol(x[t,,]))
  }
}

if(show){
  
##   png(paste("particle_phase",pname,sep="_"))
##   for (i in 1:N){
##     plot(x[,i,2],x[,i,3], type='l',xlim=psize,ylim=psize, axes=F,xlab="",ylab="")
##     par(new=T)
##   }
##   plot(y[,2],y[,3], type='l',xlim=psize,ylim=psize,col="red",ann=T)
##   dev.off()

  png(paste("particle_time",pname,sep="_"))
  for (i in 1:N){
    plot(seq(1,times),x[,i,2],ylim=trange ,type='l', axes=F,xlab="",ylab="")
    par(new=T)
  }
  plot(seq(1,times,by=interval),y[,2],ylim=trange ,type='b',col="red",ann=T)
  dev.off()
  
  xm <- array(dim=times/interval)
  xv <- array(dim=times/interval)
  for(ta in seq(1,times,by=interval)){
    xm[ta/interval] <- mean(x[ta,,2])
    xv[ta/interval] <- var(x[ta,,2])
  }

  png(paste("particle_time_dif",pname,sep="_"))
  plot(seq(1,times/interval),xm-y[,2], type='l')
  dev.off()

  png(paste("particle_time_var",pname,sep="_"))
  plot(seq(1,times/interval),xv, type='l',col="green")
  dev.off()

}
   
  }
}


