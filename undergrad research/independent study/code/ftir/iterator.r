
#example:
# > state <- c(0,0)
# > sizes <- c(2,2)
# > iterator.increment(state, sizes)
# 0 1
# > iterator.increment(state, sizes)
# 0 2
# > iterator.increment(state, sizes)
# 1 0

iterator.increment <- function(state, sizes)
{
  if (length(state) != length(sizes))
  {
    print("unmatched state/size vectors")
    return(NA)
  }
  did_inc <- FALSE
  i <- 1
  while( i <= length(sizes) && !did_inc )
  {
    if (state[i] == sizes[i])
    {
      if (i == length(state))
      {
        state[1:length(state)] <- 1 
      }
      else
      {
        state[i] <- 1
      }
    }
    else
    {
      state[i] <- state[i] + 1
      did_inc <- TRUE 
    }
    i <- i + 1
  }
  state  
}