rossler.params <- list(a = 0.1,
               b = 0.1,
               c = 18)

rossler.noise <- function(t, y,v,p=lotz.params){
  y <- y+rnorm(3,0,v)
  a <- p[['a']]
  b <- p[['b']]
  c <- p[['c']]

  c(-y[2]-y[3],
    y[1] +a*y[2],
    b+y[3] *(y[1]-c)
    )
}
