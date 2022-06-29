jk.logit <- function(x,...)
{
  y = log10( x / (1 - x) )
  return(y)
}
  