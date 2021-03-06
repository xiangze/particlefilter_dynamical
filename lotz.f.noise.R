lotz.params <- list(dl = 10,
               bt = 8/3,
               ro = 28)

lotz.noise <- function(t, y,v,p=lotz.params){
  y <- y+rnorm(3,0,v)
  dl <- p[['dl']]
  ro <- p[['ro']]
  bt <- p[['bt']]

#  list(
       c(dl * (y[2] - y[1]),
           y[1] * (ro - y[3]) - y[2],
           y[1] * y[2] - bt * y[3] )
#)
}

