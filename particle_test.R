#library(odesolve)
#source("lotz.f.noise.R")

dt <- 0.01
N <- 30
times <- 2000
interval <- 20

lag <- 1

rv <- 0.01 #original noise
xv <- 0.1  #system noise var
yv <- 3.   #obesebation noise var
nv <- 0.01 # of PF

psize <- c(-40,40)
preshow <- T
show <- T

lotz.params <- list(dl = 10,
               bt = 8/3,
               ro = 28)

lotz.noise <- function(t, y,v,p=lotz.params){
  y <- y+rnorm(3,v)
  dl <- p[['dl']]
  ro <- p[['ro']]
  bt <- p[['bt']]

       c(dl * (y[2] - y[1]),
           y[1] * (ro - y[3]) - y[2],
           y[1] * y[2] - bt * y[3] )
}


resample<- function(inits,costs){
   inits.n <- c()
    l <- length(costs)

   for (j in 1:l){
      nn <- costs[j]

      if(nn>=1)
        for(kk in 1:nn){
          inits.n <- rbind(inits.n,inits[j,])
        }

      if((nn%%1)>runif(1,0,1)){
        inits.n <- rbind(inits.n,inits[j,])
      }
    }
   inits.n <- rbind(inits.n,inits.n)
   inits.n <-inits.n[1:l,]
}

lz <- function(x,v){x+dt*lotz.noise(0,x,v)}
#f <- function(x){x+dt*lotz.noise(0,x,xv)}

#make original(hidden) time series (x) with system noise
#observed sequence :y
y <- array(dim=c(times/interval,3))
xo <- array(dim=c(1,3))
xo[1,] <- c(0.2,0.2,0.1)

for (t in 2:times){  
  xo <- rbind(xo,lz(xo[t-1,],rv))
  if(t%%(interval)==0){
    y[t/interval,] <-xo[t,]+rnorm(3,yv);
  }
}

if(preshow){
  png("org_phase.png")
  plot(y[,2],y[,3], type='l',xlim=psize,ylim=psize,col="red",xlab="x",ylab="y")
  par(new=T)
  plot(xo[,2],xo[,3], type='l',xlim=psize,ylim=psize, axes=F,xlab="",ylab="")
  dev.off()

  png("org_time.png")
  plot(seq(1,times),xo[,2], type='l', axes=F,xlab="",ylab="")
  par(new=T)
  plot(seq(1,times,by=interval),y[,2],type='l',col="red",xlab="x",ylab="y")
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
        weight[n] <- weight[n]-dnorm(y[(t/interval-i),],x[(t-i*interval),n,],yv,log=T)
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
  png("particles_phase.png")
  for (i in 1:N){
    plot(x[,i,2],x[,i,3], type='l',xlim=psize,ylim=psize, axes=F,xlab="",ylab="")
    par(new=T)
  }
  plot(y[,2],y[,3], type='l',xlim=psize,ylim=psize,col="red",ann=T)
  dev.off()

  png("particles_time.png")
  for (i in 1:N){
    plot(seq(1,times),x[,i,2], type='l', axes=F,xlab="",ylab="")
    par(new=T)
  }
  plot(seq(1,times,by=interval),y[,2],type='l',col="red",ann=T)
  dev.off()
}
