build.range <- function(range, f, idxs=c(1,length(range)))
{
  ranges <- matrix(nrow=0, ncol=2)
  start_range <- 0
  end_range <- 0
  for (idx in idxs[1]:idxs[2])
  {
    result <- f(range, idx)
    if (result == TRUE)
    {
      if (start_range == 0)
      {
        start_range <- idx
        end_range <- idx
      }
      else
      {
        end_range <- idx
      }
    }
    if (result == FALSE)
    {
      if (start_range != 0 && start_range != end_range)
      {
        ranges <- rbind(ranges, c(start_range, end_range) )
        start_range <- 0
        end_range <- 0
      }
      
    }
  }
  ranges
}

var.export <- function(varlist)
{
  o <- list()
  for (var in varlist)
  {
    o[[var]] <- get( var, envir=parent.frame(), inherits=FALSE )
  }
  o
}