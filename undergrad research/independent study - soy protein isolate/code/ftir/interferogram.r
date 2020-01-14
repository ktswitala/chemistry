
source(paste(ftir_path,"misc.r",sep=""))

calc.lorentzian <- function(v, v0, aL, scaling=1)
{
  scaling * (1/pi) * (aL / ( (v - v0)^2 + aL^2))
 }

add.lorentzian <- function(spec, v0, aL, scaling=1)
{
  for( i in 1:length(spec[,1]) )
  {
    spec[,2] <- spec[,2] + calc.lorentzian(spec[,1], v0, aL, scaling=scaling)
  }
  spec
}

spectrum.magnitude <- function(spectrum)
{
  abs( spectrum )
}

spectrum.phase <- function(spectrum)
{
  atan( Im(spectrum) / Re(spectrum) )
}

#finds the smooth parts of the phase spectrum
phase.good.filter <- function(phase_spectrum)
{
  f <- c(-1,-1,4,-1,-1)
  temp <- filter(phase_spectrum, f)
  temp[which(is.na(temp))] <- FALSE
  temp
}

phase.good <- function(phase_spectrum, th)
{
  phase.good.filter(phase_spectrum) < th
}

# discards half of the data
fft.discard <- function(spectrum)
{
  save_freq <- floor( length(spectrum)/2 )
  new_spectrum <- matrix(ncol=1, nrow=save_freq )
  new_spectrum[,1] <- spectrum[2:(1+save_freq)] 
  new_spectrum
}

# shifts data to the right with a cyclic boundary
ig.shift <- function(ig, n)
{
  new_ig <- rep(0, length(ig) )
  src_bound <- length(ig)-n
  new_ig[1:n] <- ig[(src_bound+1):length(ig)]
  new_ig[(n+1):length(ig)] <- ig[1:src_bound]
  new_ig
}

# inserts m interpolation points into the wavenumbers
ig.zerofill.wn <- function(wavenumbers, m)
{
  interp_ct <- m*(length(wavenumbers)-1)
  new_size <- interp_ct + length(wavenumbers)
  new_wn <- rep(0, new_size)
  idx <- 1
  for ( wn_idx in 1:(length(wavenumbers)-1) )
  {
    new_wn[idx] <- wavenumbers[wn_idx]
    idx <- idx + 1
    incr <- (wavenumbers[wn_idx+1] - wavenumbers[wn_idx]) / (m+1)
    for ( i in 1:(m) )
    {
      new_wn[idx] <- new_wn[idx-1] + incr
      idx <- idx + 1
    }
  }
  new_wn[new_size] <- wavenumbers[length(wavenumbers)]
  new_wn
}

ig.zerofill <- function(ig, m)
{
  new_size <- (2^m) * length(ig)
  new_ig <- rep(0, new_size )
  new_ig[1:length(ig)] <- ig[1:length(ig)]
  new_ig
}

ig.opus.extract <- function(ifg)
{
  # discard time data
  ifg <- ifg[,2]

  # forward and backward scans
  mid <- (length(ifg)/2)
  fw <- ifg[1:mid]
  bw <- ifg[(mid+1):length(ifg)]
  var.export( c("fw", "bw", "mid") )
}

ig.opus.fft <- function(o_ifg)
{
  fft_pts <- floor(length(o_ifg$fw) / 2)
  
  wvn_increment <- lwn / (16384/res)
  
  fw_fft <- fft(o_ifg$fw)
  bw_fft <- fft(o_ifg$bw)
  fw_fft_trunc <- matrix(nrow=fft_pts, ncol=2)
  bw_fft_trunc <- matrix(nrow=fft_pts, ncol=2)
  fw_fft_trunc[,1] <- (1:fft_pts)*wvn_increment
  fw_fft_trunc[,2] <- fw_fft[2:(1+fft_pts)]
  bw_fft_trunc[,1] <- (1:fft_pts)*wvn_increment
  bw_fft_trunc[,2] <- bw_fft[2:(1+fft_pts)]
  
  c(o_ifg, var.export( c("fw_fft", "bw_fft", "fw_fft_trunc", "bw_fft_trunc") ) )
}

ig.broken <- function(ifg)
{
  #find centerburst
  ifg_centerburst_idx <- which.max( abs(ifg) )

  # transform to frequency domain and discard the last half + a0
  true_spectrum <- fft( ifg )
  spectrum <- fft.discard( true_spectrum )

  # calculate phase spectrum
  magnitude_spectrum <- spectrum.magnitude( spectrum )
  phase_spectrum <- spectrum.phase( spectrum )

  # find the reasonable parts of the phase  
  good_phase <- phase.good( phase_spectrum, 0.1 )

  good_phase_ranges <- build.range(good_phase, function(range,idx) { return(range[idx] == TRUE && range[idx+1] == TRUE) }, c(1,length(good_phase)-1) )
  spec_range <- good_phase_ranges[ which(good_phase_ranges[,2] - good_phase_ranges[,1] > 1000), ]
  if (length(spec_range) != 2)
  {
    print("Error: found more than one acceptable range")
  }
  var.export( ls() )
}

ig.ratio <- function(o_fg, o_bg)
{
  o_fg$ratio_spectrum <- -log10(o_fg$magnitude_spectrum / o_bg$magnitude_spectrum)
  o_fg
}

split.freqs <- function(sz)
{
  l <- list()
  if (sz %% 2 == 1)
  {
    l$positive <- 2:(floor(sz/2)+1) 
    l$negative <- (floor(sz/2)+2):sz
  } else 
  {
    l$positive <- 2:(floor(sz/2)+1) 
    l$negative <- (floor(sz/2)+1):sz
  }
  l
}
