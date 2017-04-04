## ACTUALLY ROATING THIS MATRIX IS NOT TRIVIAL ....
plotRotatedMatrix <- function(x,plot.new = TRUE,col.palette = heat.colors(10),...)
{
  ## clear lower triangle
  #x[lower.tri(x)] <- NA
  
  ## calculate diag
  nr <- nrow(x)
  nc <- ncol(x)
  d1 <- sqrt(nr^2 + nc^2)
  d2 <- 0.5 * d1
  
  ## empty plot area
  if(plot.new == TRUE)
  {
    plot(0, type="n", xlim=c(0, d1), ylim=c(-nrow(x)/2,nrow(x)/2),
         xlab="", ylab="", asp = 1, xaxs = "i",yaxs = "i",
         axes = F)
  }
  ## plot matrix and rotate 45
  x <- matrix(col.palette[cut(x,length(col.palette))],nc=ncol(x),byrow=T)
  plot.raster <- as.raster(x)
  rasterImage(plot.raster,xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
              interpolate=FALSE, angle=45, xpd = T)
  
}

