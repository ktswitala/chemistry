
search.peaks <- function(spectrum, v_lo, v_hi)
{
  index_lo <- which(abs(spectrum[,1]-v_lo)==min(abs(spectrum[,1]-v_lo)))
  index_hi <- which(abs(spectrum[,1]-v_hi)==min(abs(spectrum[,1]-v_hi)))
  if (index_lo > index_hi)
  {
    temp <- index_lo
    index_lo <- index_hi
    index_hi <- temp
  }
  spec_min <- min( spectrum[,2] )
  spec_max <- max( spectrum[,2] )
  th_pos <- spec_max - spec_min
  th_step <- (spec_max - spec_min) / 100

  w_pos <- (index_hi - index_lo) / 4
  w_step <- (index_hi - index_lo) / 400
  
  for (i in 1:100)
  {
    w_pos_actual <- round(w_pos)
    if (w_pos_actual >= 1)
    {
      do_w <- c(do_w, w_pos_actual)
    }
    w_pos <- w_pos - w_step
    th_pos <- th_pos - th_step
  }

  peak_results <- matrix(nrow=100, ncol=100)
  for (i in 1:100)
  {
    for (j in 1:100)
    {
      peaks <- find.localmax(spectrum[index_lo:index_hi,1], spectrum[index_lo:index_hi,2], w_pos, th_pos)
      peak_results[i,j] <- length( peaks )
    }
    th_pos <- th_pos - th_step
  }
}