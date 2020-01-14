
read.fluor <- function(fn, name) 
{
  df <- list()
  df$data <- read.table(fn)
  df$name <- name
  df
}

find.em.max <- function(dat)
{
  names <- NULL
  wvn <- NULL
  for( d in dat )
  {
    names <- c(names, d$name)
    idx <- convert.rangetoindex( d$data[,1],c(300,400) )
    wvn <- c(wvn, d$data[ (idx[1]-1) + which.max(d$data[idx,2]) , 1 ] )
  }
  data.frame(names=names, wvn=wvn)
}

make.graph <- function(filenames, legends, saveto, scalings=NULL, xlimits=c(250,800), titlebar=NULL)
{
  rb <- rainbow(length(filenames))
  par(new=FALSE)
  for(i in 1:length(filenames))
  {
    spec <- read.table(filenames[i], header=FALSE)
    if (!is.null(scalings))
    {
      spec[,2] <- spec[,2] / scalings[i]
    }
    plot(spec, xlim=xlimits, ylim=c(0,100), xlab="wavelength(nm)", ylab="intensity", type="l", col=rb[i])
    par(new=TRUE)
  }
  legend(400, 100, legends, lty=rep(1,length(filenames)), col=rb)
  if (!is.null(titlebar))
  {
    title(main=titlebar)
  }
  savePlot(saveto)
  par(new=FALSE)
}

setwd("F:\\s4\\research\\data\\fluor")

spi1 <- paste(".\\428\\",c("spi_1700","spi_1300","spi_1_00","spi_6_00"),".ASC",sep="")
make.graph( spi1, 
            c("1mg","3mg","10mg","20mg"), "C:\\spi.png")

spi2 <- paste(".\\422\\",c("spi2_300","spi2_400","spi2_500"),".ASC",sep="")
make.graph( spi2, 
            c("66mg","88mg","111mg"), "C:\\spi2.png")

make.graph( c(spi1, spi2), c("1mg","3mg","10mg","20mg","66mg","88mg","111mg"), "C:\\spi3.png", titlebar="Excitation spectrum of SPI",
            scalings=c(1,1,1,1,2.77,2.77,2.77))

make.graph( paste(".\\422\\",c("130_8_00","130_9_00"),".ASC",sep=""),
            c("62mg","84mg"), "C:\\spi130.png")
            
make.graph( paste(".\\422\\",c("hy2_1_00","hy2_5_00"),".ASC",sep=""),
            c("65mg","122mg"), "C:\\hy2.png")

make.graph( paste(".\\422\\",c("spi2_300","130_6_00","hy2_1_00"),".ASC",sep=""),
            c("spi","spi130","hydrate"), "C:\\wt60mg.png", scalings=c(1,0.5,1))

make.graph( paste(".\\422\\",c("spi2_500","130_6_00","hy2_5_00"),".ASC",sep=""),
            c("spi 111mg","spi130 60mg","hydrate 120mg"), "C:\\wt60-120mg.png", scalings=c(1,0.5,1))

make.graph( paste(".\\422\\",c("spi2_300","130_6_00","hy2_1_00","ph_1__00","pc_1__00"),".ASC",sep=""),
            c("spi","spi heated","spi hydrate","pressed hot","pressed cool"), "C:\\all.png", scalings=c(1,0.5,1,2,2), xlimits=c(300,600), titlebar="Excitation spectrum @ ~60mg protein")
            
make.graph( paste(".\\428\\",c("spi_1700", "1302_100"),".ASC",sep=""),
            c("spi", "spi heated"), "C:\\low_conc.png", titlebar="Excitation spectrum @ ~1mg protein")
           