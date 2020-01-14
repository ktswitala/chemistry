
dev.setcount <- function(c)
{
  while (length( dev.list() ) > c )
  {
    dev.off()
  }
  while (length( dev.list() ) < c )
  {
    dev.new()
  }
  dev.set( dev.list()[1] )
}