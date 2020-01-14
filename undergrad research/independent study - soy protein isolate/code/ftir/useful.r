linreg.compute <-
function (xi, yi) {
  mx = mean(xi)  
  my = mean(yi)
  Nx = length(xi)
  sxx = sum( (xi - mx)^2 )
  syy = sum( (yi - my)^2 )
  sxy = sum( (xi - mx)*(yi - my) )
  m = sxy / sxx
  b = my - m*mx
  sr = sqrt( (syy - m^2*sxx) / (Nx-2) )
  sm = sqrt( sr^2 / sxx )
  sb = sr * sqrt( 1 / (Nx - sum( xi )^2 / sum( xi^2 ) ) )
  
  daf = data.frame(mx=NA, my=NA, Nx=NA, sxx=NA, syy=NA, sxy=NA, m=NA, b=NA, sr=NA, sm=NA, sb=NA)
  daf$mx = mx
  daf$my = my
  daf$Nx = Nx
  daf$sxx = sxx
  daf$syy = syy
  daf$sxy = sxy
  daf$m = m
  daf$b = b
  daf$sr = sr
  daf$sm = sm
  daf$sb = sb
  
  return( daf )
}

err.add <- function (x, y) {
    myret = NULL
    myret[1] = x[1] + y[1]
    myret[2] = sqrt( x[2]^2 + y[2]^2 )
    return(myret)
}

err.sub <- function (x, y) {
    myret = NULL
    myret[1] = x[1] - y[1]
    myret[2] = sqrt( x[2]^2 + y[2]^2 )
    return(myret)
}

err.mul <- function (x,y) {
    myret = NULL
    myret[1] = x[1] * y[1]
    myret[2] = sqrt( (y[1]*x[2])^2 + (x[1]*y[2])^2 ) 
    return(myret)
}

err.mul2 <- function (x, y) {
    myret = NULL
    myret[1] = x[1] * y[1]
    myret[2] = sqrt( (x[2]/x[1])^2 + (y[2]/y[1])^2 ) * (x[1]*y[1])
    return(myret)
}

err.div <- function (x, y) {
    myret = NULL
    myret[1] = x[1] / y[1]
    myret[2] = sqrt( ((1/y[1])*x[2])^2 + ((-x[1]/y[1]^2)*y[2])^2 )
    return(myret)
}

err.div2 <- function (x, y) {
    myret = NULL
    myret[1] = x[1] / y[1]
    myret[2] = sqrt( (x[2]/x[1])^2 + (y[2]/y[1])^2 ) * (x[1]/y[1])
    return(myret)
}

linreg.solve <- function(r, v)
{
  v <- c(v, 0)
  t <- err.sub(v, c(r$b, r$sb))
  t <- err.div(t, c(r$m, r$sm))
  return ( t )
}

err.add3 <-
function (ex, ey)
    return( sqrt( ex^2 + ey^2 ) )
    
err.sub3 <-
function (ex, ey)
    return( sqrt( ex^2 + ey^2 ) )

err.mul3 <-
function (x, ex, y, ey)
    return( sqrt( (ex/x)^2 + (ey/y)^2 ) * (x*y) )

err.div3 <-
function (x, ex, y, ey)
    return( sqrt( (ex/x)^2 + (ey/y)^2 ) * (x/y) )

find.localmax <- function(xs,ys,wid,th) {
  peaks = array(dim=c(0,2))
  if (length(xs) != length(ys))
  {
    print("xs is not equal in length to ys")
    return(NA)
  }
  if (1+wid > length(xs)-wid)
  {
    print("width parameter too large")
    return(NA)
  }
  for ( i in (1+wid):(length(xs)-wid) )
  {
    tmax <- max( c(ys[(i-wid):(i-1)], ys[(i+1):(i+wid)]) )
    tmin <- min( c(ys[(i-wid):(i-1)], ys[(i+1):(i+wid)]) )
    if ( (ys[i] > tmax) && (ys[i] > tmin+th) )
    {
      peaks <- rbind(peaks, c(xs[i], ys[i]))
    }
  }
  peaks
}

df.avg <- function(d, s, ind) {
  t = array(rep(0,s))
  for (i in ind)
  {
    t <- t + d[[i]]
  }
  t = t / 3
  return( t )
}

interpol <- function(cx, lx, hx, ly, hy)
{
  ly + (cx-lx)/(hx-lx) * (hy-ly)
}

nearest <- function(vec, val)
{
  for (i in 1:(length(vec)))
  {
    if (vec[i] > val)
    {
      return(i-1)
    }
  }
}

deriv <- function(x, y)
{
  dydx <- vector()
  newx <- vector()
  for (i in 1:(length(x)-1))
  {
    newx <- c(newx, (x[i] + x[i+1]) / 2)
    dydx <- c(dydx, (y[i+1] - y[i]) / (x[i+1] - x[i])) 
  }
  array(data=c(newx,dydx), dim=c(length(x)-1, 2))
}

# http://users.fmg.uva.nl/rgrasman/rpages/2005/09/error-bars-in-plots.html
superpose.eb <- 
function (x, y, ebl, ebu = ebl, length = 0.08, ...) 
    arrows(x, y + ebu, x, y - ebl, angle = 90, code = 3, 
    length = length, ...)

###

