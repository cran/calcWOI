\name{calcWOI-package}
\alias{calcWOI-package}
\alias{calcWOI}
\docType{package}
\title{
Calculates the Wavelet-Based Organization Index
}
\description{
Calculates the original wavelet-based organization index, the modified wavelet-based organization index and the local wavelet-based organization index of an arbitrary 2D array. Since version 1.0.3 the function LW is added, which calculates the local wavelet-based organization index with help of Dual-tree wavelets.
}
\details{
This package provides all functions and tools to calculate the wavelet-based organization index (Brune et al., 2018). The central function within the calcWOI package is WOI. This function calculates the original WOI, the modified WOI and the local WOI (Brune et al., 2020a). All these indexes based on wavelet transforms (DaubExPhase4) done with the function wavtra, where we use parts of the cddews function within the package LS2W to create the function wavtra and constant data like the inverse correction matrix AI and the centre of mass matrix.

Other functions in calcWOI like flatten, blowup and buildperiodic are used to ensure that the incoming field is quadratic and periodic.

The function LW allows the calculation of the revised LWOI, namely LW, of an arbitrary array. Details on the new LW will be published in Brune et al (2020b).
}
\author{
Sebastian Brune, Sebastian Buschow, Florian Kapp, Petra Friederichs. Maintainer: Sebastian Brune <sbrune@uni-bonn.de>
}
\references{
Eckley, I.A., Nason, G.P. and Treloar, R.L. (2010) Locally stationary wavelet fields with application to the modelling and analysis of image texture. Journal of the Royal Statistical Society (Series C), 59, 595-616.

Eckley, I.A. and Nason, G.P. (2011). LS2W: Implementing the Locally Stationary 2D Wavelet Process Approach in R, Journal of Statistical Software, 43(3), 1-23.

Brune, S., Kapp, F., & Friederichs, P. (2018). A wavelet-based analysis of convective organization in ICON large-eddy simulations. Quarterly Journal of the Royal Meteorological Society, 144(717), 2812-2829.

Brune, S., Buschow, S., & Friederichs, P. (2020a). Observations and high-resolution simulations of convective precipitation organization over the tropical Atlantic. Quarterly Journal of the Royal Meteorological Society.

Brune, S., Buschow, S., & Friederichs, P. (2020b). The Local Wavelet-based Organization Index - Quantification, Localization and Classification of Convective Organization from Radar and Satellite Data. In preparation.
}
\examples{
# Calculate WOI, modified WOI and LWOI for a random precipitation
# field using an 230x200 array with
# random positive numbers 
x <- array(5 + rnorm(230*200), dim = c(230, 200))
s <- c(1,2)
l <- c(3,4)
thres <- 0.1
flat <- 5
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
par(mfrow = c(2, 2))
image(LWOI, main = "LWOI")

# Calculate the three components of the 
# local wavelet-based organization index with Dual-tree wavelets
LW <- LW(x, thres = thres, Nx = 2^ceiling(log2(max(dim(x)))), 
  Ny = 2^ceiling(log2(max(dim(x)))), boundaries = "pad")
image(LW$LWsc, main = "LWsc", zlim = 0:1)
image(LW$LWin, main = "LWin", zlim = 0:1)
image(LW$LWai, main = "LWai", zlim = 0:1)
}
