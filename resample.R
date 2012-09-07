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
