\name{LW}
\alias{LW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates the wavelet-based organization index with the dualtree wavelet transform
%%  ~~function to do ... ~~
}
\description{This function calculates the locally wavelet-based organization index (LW) as defined in Brune et al. (2020) based on the dualtree complex wavelet transform.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LW(x, thres = 0.1, Nx = 2^ceiling(log2(max(dim(x)))), 
  Ny = 2^ceiling(log2(max(dim(x)))), boundaries = "pad", verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{2D numeric array with dimensions larger than 16 x 16. dim(x)[1] has not to be equal to dim(x)[1], but NA and/or NaN are not allowed.
%%     ~~Describe \code{x} here~~
}
  \item{thres}{0 or a positive number. Threshold for rain rate. Default is 0.1 mm/h, because we calculate LWOI only for grid points, where rain rate is >= 0.1 mm/h. For brightness temperatures we use 245 K.
%%     ~~Describe \code{thres} here~~
}
  \item{Nx}{The number of grid points in x direction of the array, which is put into the dualtree wavelet transform. The default is the log2 of maximum dimension of the ingoing array. Nx must be equal or larger than the first dimension of x.
%%     ~~Describe \code{Nx} here~~
}
  \item{Ny}{The number of grid points in y direction of the array, which is put into the dualtree wavelet transform. The default is the log2 of maximum dimension of the ingoing array. Nx must be equal or larger than the second dimension of x.
%%     ~~Describe \code{Ny} here~~
}
  \item{boundaries}{Handles the boundary conditions, either "pad", "mirror" or "periodic".
%%     ~~Describe \code{boundaries} here~~
}
  \item{verbose}{Default FALSE. Set TRUE for print statements.
%%     ~~Describe \code{verbose} here~~ 
}
}
\value{This functions returns a list with follwing elements:
\item{LWsc}{LWOI scale calculated with dualtrees. The array is masked with thres. The dimension of LWsc is equal to the dimension of x, if the boundaries are only padded with zeros.}
\item{LWin}{LWOI intensity calculated with dualtrees. The array is masked with thres. The dimension of LWsc is equal to the dimension of x, if the boundaries are only padded with zeros.}
\item{LWai}{LWOI anisotropy calculated with dualtrees. The array is masked with thres. The dimension of LWsc is equal to the dimension of x, if the boundaries are only padded with zeros.}
\item{LWuu}{LWOI u component calculated with dualtrees. The array is masked with thres. The dimension of LWsc is equal to the dimension of x, if the boundaries are only padded with zeros.}
\item{LWvv}{LWOI v component with dualtrees. The array is masked with thres. The dimension of LWsc is equal to the dimension of x, if the boundaries are only padded with zeros.}
\item{angle}{Angle of spectrum. The array is masked with thres. The dimension of angle is equal to the dimension of x, if the boundaries are only padded with zeros. Angle ranges between 0 degree and 180 degree.}
\item{thres}{Threshold of LW calculation.}
\item{mask}{Mask defined by thres.}
\item{x}{Ingoing array.}
\item{ts}{Computation time in seconds.}
}
\references{Brune, S., Kapp, F., & Friederichs, P. (2018). A wavelet-based analysis of convective organization in ICON large-eddy simulations. Quarterly Journal of the Royal Meteorological Society, 144(717), 2812-2829.
Brune, S., Buschow, S., & Friederichs, P. (2020). Observations and high-resolution simulations of convective precipitation organization over the tropical Atlantic. Quarterly Journal of the Royal Meteorological Society.
Brune, S., Buschow, S., & Friederichs, P. (2020). The Local Wavelet-based Organization Index - Quantification, Localization and Classification of Convective Organization from Radar and Satellite Data.
%% ~put references to the literature/web site here ~
}
\author{Sebastian Brune
%%  ~~who you are~~
}
\note{This function calculates the locally wavelet-based organization on the basis of dualtree wavelet spectra. To calculate the WOI, modified WOI or LWOI as presented in Brune et al (2020), use the function WOI.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{The input array must be numeric without NA or NaN and maximal of size 1024 x 1024.}

\examples{
# Random array of dim 230 x 200
x <- array((rnorm(230*200)), dim = c (230, 200))
thres <- 0.1
LWres <- LW(x, thres = thres, Nx = 2^ceiling(log2(max(dim(x)))), 
  Ny = 2^ceiling(log2(max(dim(x)))), boundaries = "pad")
print(paste("The LWOI calculation took", LWres$ts, "seconds."))
# plot data
par(mfrow = c(3, 2))
# original data
image(LWres$x, main = "Original data", zlim = 0:1)
# scale component
image(LWres$LWsc, main = "Scale", zlim = 0:1)
# intensity component
image(LWres$LWin, main = "Intensity", zlim = 0:1)
# anisotropy component
image(LWres$LWai, main = "Anisotropy", zlim = 0:1)
# u component
image(LWres$LWuu, main = "u direction", zlim = 0:1)
# v component
image(LWres$LWvv, main = "v direction", zlim = 0:1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
