
source("F:\\s4\\research\\code\\paths.r")
source(paste(ftir_path,"interferogram.r",sep=""))
source(paste(ftir_path,"scaling.r",sep=""))

deconv <- function(spec, big_T, sigma_start, sigma_incr, new_size_scale=1, xlims=c(50,150), title1="", title2="")
{
  plot(spec, type="l", xlim=xlims, xlab="frequency", ylab="intensity")
  title(main=title1)
  savePlot("C:\\ftir\\deconvolution\\spec-0.0.png")

  # FT
  spec_fft_x <- 1:length(spec[,1])
  spec_fft_y <- fft(spec[,2]) / length(spec[,2])
  plot(spec_fft_x,spec_fft_y,type="l")
  savePlot("C:\\ftir\\deconvolution\\spec_fft-0.0.png")

  for (i in 1:30)
  {
    freq_index <- split.freqs(length(spec_fft_x))
    positive <- freq_index$positive
    negative <- freq_index$negative

    # deconvolution
    scaling <- (i*sigma_incr) + sigma_start
    #big_T <- (.301675*K) / scaling
    boxcar <- function(t)
    {
      if (t > big_T)
      {
        return(0)
      }
      else
      {
        return(1)   
      }
    }
    triangle <- function(t)
    {
      if (t > big_T)
      {
        return(0)
      }
      else
      {
        return (1 - (t/big_T))
      }
    }

    t <- seq(0,1,by=(1/(length(positive)-1)))
    t_filter <- exp(2*pi*(scaling)*t) * vapply( t, triangle, 0 )^2
    #t_filter <- exp(2*pi*(scaling)*t) * (1-(t/big_T)^2)^2 * vapply( t, boxcar, 0 )
    #t_filter <- exp(2*pi*(scaling)*t) * vapply( t, boxcar, 0 )

    spec_fft2_y <- spec_fft_y
    spec_fft2_y[positive] <- spec_fft_y[positive] * t_filter
    spec_fft2_y[negative] <- Conj(rev(spec_fft2_y[positive]))
    plot(spec_fft_x,Re(spec_fft2_y),type="l")
    savePlot(paste("C:\\ftir\\deconvolution\\spec_fft-",format(scaling,nsmall=1),".png",sep=""))

    # interpolation
    wvn_new <- spec[,1]
    sz <- length(spec[,2])
    sz_i <- new_size_scale
    while (sz_i > 1)
    {
      wvn_new <- ig.zerofill.wn( wvn_new, 1 )
      sz <- sz + (sz-1)
      sz_i <- sz_i-1
    }
    new_size <- sz
    zerofilled <- rep(0, new_size)
    zerofilled[1] <- spec_fft2_y[1]
    zerofilled[positive] <- spec_fft2_y[positive]
    zerofilled[negative+(new_size-length(spec[,1]))] <- spec_fft2_y[negative]

    # FT -1
    spec2 <- fft(zerofilled, inverse=TRUE)
    plot(wvn_new, Re(spec2), type="l", xlim=xlims, xlab="frequency", ylab="intensity")
    title(main=title2)
    savePlot(paste("C:\\ftir\\deconvolution\\spec-", format(scaling,nsmall=1),".png",sep=""))
  }
}

spec_f <- read.table("C:\\ftir\\opus-standardized\\april 16 spi hydrate.0.ab");
indices <- convert.rangetoindex( spec_f[,1], c(1600,1690) )
spec <- cbind( spec_f[indices,1], spec_f[indices,2] )
#deconv(spec, 0.45, 0, 0.2, xlims=c(1600,1690), new_size_scale=4, title1="SPI hydrated - Original", title2="SPI hydrated - Deconvoluted")

spec_f <- read.table("C:\\ftir\\opus-standardized\\april 8 spi.3.ab");
indices <- convert.rangetoindex( spec_f[,1], c(1600,1690) )
spec <- cbind( spec_f[indices,1], spec_f[indices,2] )
deconv(spec, 0.42, 0, 0.2, xlims=c(1600,1690), new_size_scale=4, title1="SPI - Original", title2="SPI - Deconvoluted")

spec <- matrix(nrow=188, ncol=2)
spec[,1] <- 1:length(spec[,1])
spec[,2] <- 0
spec <- add.lorentzian(spec, 100, 4)
spec <- add.lorentzian(spec, 106, 4, scaling=0.75)
#deconv(spec, 0.8, 0, 0.1, new_size_scale=4, title1="Test peaks - Original", title2="Test peaks - Deconvoluted")
