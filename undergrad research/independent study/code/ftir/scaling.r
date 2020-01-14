
convert.rangetoindex <- function(column, ranges)
{
  lo_index <- which(abs(column-ranges[1])==min(abs(column-ranges[1])))
  hi_index <- which(abs(column-ranges[2])==min(abs(column-ranges[2])))
  if (lo_index > hi_index)
  {
    hi_index:lo_index
  }
  else
  {
    lo_index:hi_index
  }
}

find.common.wavenumbers <- function(spec_large, spec_small)
{
  #assumes both spectrums are decreasing
  start_index <- 1
  large_indices <- NULL
  small_indices <- NULL
  for (i in 1:length(spec_small))
  {
    for (j in start_index:length(spec_large))
    {
      if (abs(spec_small[i] - spec_large[j]) < 0.01)
      {
        small_indices <- c(small_indices, i)
        large_indices <- c(large_indices, j)
      }
    }
  }
  print (length(small_indices) / length(spec_small) )
}

match.spectra.basic <- function(spectrum1, spectrum2, indexes, scale_steps=100)
{
  spectrum_ratio <- (spectrum1[indexes,2] / spectrum2[indexes,2])
  spectrum_ratio <- spectrum_ratio[ !is.infinite( spectrum_ratio ) ]
  scale_mean <- mean( spectrum_ratio )
  scale_stdev <- sd( spectrum_ratio )

  scale_step <- (scale_stdev * 4) / scale_steps
  scale_pos <- scale_mean - scale_stdev * 2

  scale_amounts <- vector(length=scale_steps)
  step_distances <- vector(length=scale_steps)
  for( k in 1:scale_steps )
  {
    scale_pos <- scale_pos + scale_step
    scaled_spectrum <- spectrum2[indexes,2] * scale_pos
    scale_amounts[k] <- scale_pos
    step_distances[k] <- sum( abs( scaled_spectrum - spectrum1[indexes,2] ) )    
  }
  spectrum2[,2] <- spectrum2[,2] * scale_amounts[ which.min( step_distances ) ]
  rets <- list()
  rets$spec <- spectrum2
  rets$dist <- step_distances[ which.min( step_distances ) ]
  rets
}

match.spectra <- function(spectra, ranges, scale_steps=50)
{
  spectra_count <- ncol(spectra)-1
  if (spectra_count == 1)
  {
    print("only one spectra in match.spectra")
    return(NA)
  }

  indexes <- vector()
  for(i in 1:nrow(ranges))
  {
    lo_index <- which(abs(spectra[,1]-ranges[i,1])==min(abs(spectra[,1]-ranges[i,1])))
    hi_index <- which(abs(spectra[,1]-ranges[i,2])==min(abs(spectra[,1]-ranges[i,2])))
    if (lo_index > hi_index)
    {
      indexes <- c(indexes, hi_index:lo_index)
    }
    else
    {
      indexes <- c(indexes, lo_index:hi_index)
    }
  }
  # calculate paired distances
  step_distances <- matrix(nrow=scale_steps, ncol=(spectra_count-1) )
  scale_amounts <- matrix(nrow=scale_steps, ncol=(spectra_count-1) )
  dist_index <- 1

  spectra_ref <- spectra[indexes,2]
  for( back_spectra_i in 3:(spectra_count+1) )
  {
    spectrum <- spectra[indexes,back_spectra_i]
    spectrum_ratio <- spectra_ref / spectrum
    scale_mean <- mean( spectrum_ratio )
    scale_stdev <- sd( spectrum_ratio )

    scale_step <- (scale_stdev * 4) / scale_steps
    scale_pos <- scale_mean - scale_stdev * 2

    for( k in 1:scale_steps )
    {
      scale_pos <- scale_pos + scale_step
      scaled_spectrum <- spectrum * scale_pos
      scale_amounts[k,dist_index] <- scale_pos
      step_distances[k,dist_index] <- sum( (scaled_spectrum - spectra_ref )^2 )    
    }
    dist_index <- dist_index + 1
  }

  # calculate global minimum
  state <- c( rep(1, spectra_count-1) )
  sizes <- c( rep(scale_steps, spectra_count-1) )

  min_dist <- NA
  min_state <- NA
  for( i in 1:prod(sizes) )
  {
    dist <- 0

    for (k in 1:length(state))
    {
      dist <- dist + step_distances[ state[k], k ]
    }
    for( front_spectra_i in 3:(spectra_count) )
    {
      for( back_spectra_i in (front_spectra_i+1):(spectra_count+1) )
      { 
        scale_1 <- scale_amounts[ state[front_spectra_i-2], front_spectra_i-2 ]
        scale_2 <- scale_amounts[ state[back_spectra_i-2], back_spectra_i-2 ] 
        dist <- dist + sum( abs( spectra[indexes,front_spectra_i]*scale_1 - spectra[indexes,back_spectra_i]*scale_2 ) )
      }
    }
    if (is.na(min_dist) || dist < min_dist) 
    {
      min_dist <- dist
      min_state <- state
    }
    state <- iterator.increment( state, sizes )
  }
  scaled_spectra <- spectra
  for ( i in 3:(spectra_count+1) )
  {
    scaled_spectra[,i] <- scaled_spectra[,i] * scale_amounts[ min_state[i-2], i-2 ]
  }
  scaled_spectra
}