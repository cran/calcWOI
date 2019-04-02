\name{calcWOI-package}
\alias{calcWOI-package}
\alias{calcWOI}
\docType{package}
\title{
Calculates the Wavelet-Based Organization Index
}
\description{
Calculates the original wavelet-based organization index, the modified wavelet-based organization index and the local wavelet-based organization index of an arbitrary 2D array.
}
\details{
This package provides all functions and tools to calculate the wavelet-based organization index (Brune et al., 2018). The central function within the calcWOI package is WOI. This function calculates the original WOI, the modified WOI and the local WOI. All these indexes based on wavelet transforms (DaubExPhase4) done with the function wavtra, where we use parts of the cddews function within the package LS2W to create the function wavtra and constant data like the inverse correction matrix AI and the centre of mass matrix.

Other functions in calcWOI like flatten, blowup and buildperiodic are used to ensure that the incoming field is quadratic and periodic.
}
\author{
Sebastian Brune, Sebastian Buschow, Florian Kapp, Petra Friederichs. Maintainer: Sebastian Brune <sbrune@uni-bonn.de>
}
\references{
Eckley, I.A., Nason, G.P. and Treloar, R.L. (2010) Locally stationary wavelet fields with application to the modelling and analysis of image texture. Journal of the Royal Statistical Society (Series C), 59, 595-616.

Eckley, I.A. and Nason, G.P. (2011). LS2W: Implementing the Locally Stationary 2D Wavelet Process Approach in R, Journal of Statistical Software, 43(3), 1-23.

Brune, S., Kapp, F., & Friederichs, P. (2018). A wavelet-based analysis of convective organization in ICON large-eddy simulations. Quarterly Journal of the Royal Meteorological Society, 144(717), 2812-2829.
}
\seealso{
\code{\link[LS2W:LS2W-package]{LS2W}}
}
\examples{
# Calculate WOI for a random precipitation
# field using an 350x300 array with
# random positive numbers 
x <- array(5 + rnorm(350*300), dim = c (250,300))
s <- c(1,2)
l <- c(3,4)
thres <- 0.1
flat <- 25
WOIres <- WOI (x = x, s = s, l = l, thres = thres, flat = flat,
verbose = TRUE)
# original WOI (Brune et al., 2018)
WOIorig <- WOIres$WOIorig
print(paste("Original WOI:", WOIorig))
# modified WOI
WOI <- WOIres$WOI
print(paste("Modified WOI:", WOI))
# local WOI
LWOI <- WOIres$LWOI
image(LWOI)
}
