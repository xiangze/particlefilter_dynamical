library(odesolve)
source("lotz.f.noise.R")

N <- 10
times <- 1000
interval <- 10
lag <- 3

rv <- 0.01
xv <- 0.1 #system noise var
yv <- 3. #obesebation noise var

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

resample.o  <- function(inits,costs){
    inits.o <- inits
    k <- 1
    rest <- 0
    for (j in 1:length(costs)){
      nn <- (costs[j]+rest)
      if(nn>=1){
      for(kk in 1:nn){
        inits[k,] <- inits.o[j,]
        k <- min(k+1,length(costs))
      }
      if(debug)print(k)
      rest <- (costs[j]+rest)-kk      
    }
    }
    inits
}

#make original(hidden) time series (x) with system noise
#observed sequence :y

#xo <- array(dim=c(times,3))
xo <- c(2,3,1)
for (t in 2:times){  
  xo <- rbind(xo,lotz.noise(t,xo[t-1,],rv))
  if(t%%(interval)==0){
    y <-rbind(y,xo[t,])+rnorm(3,yv);
  }
}

#y<-(select)
png("org.png")
#par(=T)
lines(xo[,2],xo[,3], type='l')
lines(y[,2],y[,3], type='l')
dev.off()

#y <- lsoda(y0, times, f, par)
#out <- ode(init[n,-1], Time, LotVmod, Pars.loop)[,-1]

#prediction by particles 
x <- array(dim=c(times,N,3))
x[1,,] <- matrix(runif(N,-1,1))
#apply (inits,1,std.lyap.run)

f <- function(x){lotz.noise(0,x,xv)}
for (t in 2:times){
  x[t,,] <- apply(x[t-1,,],1,f)

  weight <- (y[(t-lag):t,]-x[(t-lag):t,])^2
  costs <- weight/sum(weight)*N
  x[t,,] <- resample(x[t,,],costs)
  x[t,,] <- x[t,,]+matrix(rnorm(length(x[t,,]),0,var),nc=ncol(x[t,,]))
}
