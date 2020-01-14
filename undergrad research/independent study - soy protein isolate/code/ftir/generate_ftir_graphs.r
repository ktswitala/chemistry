
source("F:\\s4\\research\\code\\paths.r")
source(paste(ftir_path,"scaling.r",sep=""))

pressed_cool <- list.files("C:\\ftir\\sorted-samples\\pressed-cool\\", full.names=TRUE)
pressed_hot <- list.files("C:\\ftir\\sorted-samples\\pressed-hot\\", full.names=TRUE)
unpressed_cool <- list.files("C:\\ftir\\sorted-samples\\unpressed-cool\\", full.names=TRUE)
unpressed_hot <- list.files("C:\\ftir\\sorted-samples\\unpressed-hot\\", full.names=TRUE)
unpressed_hydrate <- list.files("C:\\ftir\\sorted-samples\\unpressed-hydrate\\", full.names=TRUE)
compare1 <- list.files("C:\\ftir\\sorted-samples\\compare1\\", full.names=TRUE)
compare2 <- list.files("C:\\ftir\\sorted-samples\\compare2\\", full.names=TRUE)
compare3 <- list.files("C:\\ftir\\sorted-samples\\compare3\\", full.names=TRUE)

process <- function(flist,picname, titlebar=NULL, wvns=c(1480,1690), legendx=1480 )
{
  spectra_ref <- read.table( flist[1] )

  sz <- length(flist)
  rb <- rainbow(sz)
  par(new=FALSE)
  i <- 1
  
  indices <- convert.rangetoindex( spectra_ref[,1], wvns )
  y_hi <- max(spectra_ref[indices,2]) * 1.1
  y_low <- min(spectra_ref[indices,2]) * 0.9
  for (f in flist)
  {

    dat <- read.table(f)
    spectra_scaled <- match.spectra.basic(spectra_ref, dat, indices)
    
    plot(spectra_scaled$spec, type="l", col=rb[i], xlim=wvns, ylim=c(y_low, y_hi), xlab="wavenumber (cm^-1)", ylab="absorbance" )
    par(new=TRUE)
    i <- i + 1
  }
  if (!is.null(titlebar))
  {
    title(main=titlebar)
  }
  legend(legendx, y_hi, basename(flist), lty=rep(1,length(flist)), col=rb)
  savePlot(paste("C:\\ftir\\sorted-samples\\",picname,".png",sep="") )
  par(new=FALSE)
}

process(pressed_cool,"presscool")
process(pressed_hot,"presshot")
process(unpressed_cool,"upresscool")
process(unpressed_hot,"upresshot")
process(unpressed_hydrate,"hydrate")
process(compare1, "compare1")
process(compare2, "compare2", titlebar="FTIR-ATR spectrum - Amide I and II bands")
process(compare3, "compare3", titlebar="FTIR-ATR spectrum - forms of soy protein", wvns=c(900,1690), legendx=1200 )
