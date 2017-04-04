## FROM GAVIN
## https://gavinband.wordpress.com/2015/04/27/computing-ld-in-r/

compute.ld <- function( haps, focus.i = 1:nrow(haps), variant.names = rownames( haps ) ) {
  # haps is a 0-1 matrix with L SNPs (in rows) and N haplotypes (in columns).
  # Since the values are 0, 1, we rely on the fact that a*b=1 iff a and b are 1.
  # assume focus.i contains l values in range 1..L
  L = nrow( haps )
  l = length( focus.i )
  # p11 = lxL matrix.  i,jth entry is probability of 11 haplotype for ith focal SNP against jth SNP.
  focus.hap = haps[ focus.i, , drop = FALSE ]
  p11 <- ( focus.hap %*% t( haps )) / ncol( haps )
  # p1. = lxL matrix.  ith row is filled with the frequency of ith focal SNP.
  p1. <- matrix( rep( rowSums( focus.hap ) / ncol( haps ), L ), length( focus.i ), L, byrow = FALSE )
  # p.1 = lxL matrix.  jth column is filled with the frequency of jth SNP.
  frequency = rowSums( haps ) / ncol( haps )
  p.1 <- matrix( rep( frequency, length( focus.i ) ), length( focus.i ), L, byrow = TRUE )
  
  # Compute D
  D <- p11 - p1. * p.1
  
  # Compute D'
  denominator = pmin( p1.*(1-p.1), (1-p1.)*p.1 )
  wNeg = (D < 0)
  denominator[ wNeg ] = pmin( p1.*p.1, (1-p1.)*(1-p.1) )[wNeg]
  denominator[ denominator == 0 ] = NA
  Dprime = D / denominator 
  
  # Compute correlation, this result should agree with cor( t(haps ))
  R = D / sqrt( p1. * ( 1 - p1. ) * p.1 * ( 1 - p.1 ))
  R[ frequency[ focus.i ] == 0 | frequency[ focus.i ] == 1, frequency == 0 | frequency == 1 ] = NA
  if( !is.null( variant.names )) {
    rownames( D ) = rownames( Dprime ) = rownames( R ) = variant.names[ focus.i ]
    colnames( D ) = colnames( Dprime ) = colnames( R ) = variant.names
    names( frequency ) = variant.names
  }
  
  return( list( D = D, Dprime = Dprime, frequency = frequency, R = R ) ) ;
}