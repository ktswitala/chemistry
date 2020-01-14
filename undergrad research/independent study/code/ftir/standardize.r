
source(paste(ftir_path,"interferogram.r",sep=""))

files <- list.files("C:\\ftir\\opus-out\\")

for (f in files)
{
  spec <- read.table(paste("C:\\ftir\\opus-out\\",f,sep=""))
  res <- (16384/15801.02) * (spec[1,1]-spec[2,1])

  if (abs(round(res)-res) > 0.01)
  {
    print(paste(res, " ",files[1],sep=""))
  }
}

do.stuff <- function(files)
{
  for (f in files)
  {
    spec <- read.table(paste("C:\\ftir\\opus-out\\",f,sep=""))
    res <- (16384/15801.02) * (spec[1,1]-spec[2,1])
    if (res >= 0.75)
    {
      doublings <- log2(round(res) / 0.5)
      wvn_new <- spec[,1]
      sz <- length(spec[,2])
      for (i in 1:doublings)
      {
        wvn_new <- ig.zerofill.wn( wvn_new, 1 )
        sz <- sz + (sz-1)
      }
      zeros <- sz - length(spec[,2])
      old_freqs <- split.freqs(length(spec[,2]))
    
      ifg <- fft(spec[,2]) / length(spec[,2])
      new_ifg <- c( ifg[1], ifg[old_freqs$positive], rep(0, zeros), ifg[old_freqs$negative])
    
      spec_new <- Re(fft(new_ifg,inverse=TRUE))
      
      if (length(wvn_new) == length(spec_new))
      {
        spec2 <- cbind(wvn_new,spec_new)
        write.table(spec2, paste("C:\\ftir\\opus-standardized\\",f,sep=""), row.names=FALSE, col.names=FALSE)
      }
      else
      {
        print(length(wvn_new))
        print(length(spec_new))
      }
      print(f)
    }
  }
}

do.stuff(files)