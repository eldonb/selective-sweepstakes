#!/usr/bin/R
library(R.utils)
## using sed to replace text :
## sed -i 's/old-text/new-text/g' file
n <- 136
vn <- log( (1:(n-1))/n ) -  log( 1 -  ((1:(n-1))/n ) )
x <- scan('path-to-data', nlines = 1)
x <- x[-c(1, n + 1)]
pdf('prefix.pdf')
plot( vn, log(x) - log(1 - x), col='black', pch=1, ylim=c(-10,3), las=1, bty='l', xlab='Derived allele frequency (logit)', ylab='Normalised site frequency (logit)', cex.lab=1.3 )
lines(vn, scan('./ds_omega_c6'), lwd=2, col='red')
points( vn,  log( x/sum(x)) - log( 1 -  (x/sum(x)) ), pch = 1, lwd=2)
                                        x <- scan('./wfsfs_treeseq_resout')
##
x <- scan('path-to-slim-sfs')
d <- data.frame( 'vy'= log(x[x > 0] ) - log(1- x[x > 0] ), 'vx' = vn[x > 0], na.action=getOption('na.action') )
f <- loess( vy ~ vx, data=d,  na.action=getOption("na.action") )
points(vn, log(x) - log(1-x), col=colours()[45], pch=6, lwd=2)
lines( vn,  predict( f,  data.frame( vx=vn) ), col=colours()[45], lwd=3, lty=5)
##
##
legend(-4,2, legend=c('south coast',  'random sweepstakes + selection', 'neg: mean = -0.05, shape = 2', 'pos: h = 1, s = 1'), pch=c(1,6,6,6), cex=1.4, bty='n', ncol=1, col=c('black',   colors()[45],colors()[45],  colors()[45] ))
legend( 3,1, legend=c('DS', 'loess'), lwd=3, lty=c(1,2), col=c('red',  colors()[45]) )
graphics.off()
