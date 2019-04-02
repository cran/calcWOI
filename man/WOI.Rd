\name{WOI}
\alias{WOI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates WOI, modified WOI and LWOI
%%  ~~function to do ... ~~
}
\description{This function calculates the wavelet-based organization index (WOI) as defined in Brune et al. (2018), a modified version of WOI and the local WOI using DaubExPhase4 wavelet.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
WOI(x = x, s = c(1, 3), l = c(4, 7), thres = 0.1, flat = 25,
    verbose = FALSE, periodic = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{2D numeric array with dimensions larger than 16 x 16. dim(x)[1] has not to be equal to dim(x)[1], but NA and/or NaN are not allowed.
%%     ~~Describe \code{x} here~~
}
  \item{s}{Vector (length 2) of smallest small convective scale and largest small convective scale. Default: s = c(1, 3).
%%     ~~Describe \code{s} here~~
}
  \item{l}{Vector (length 2) of smallest large convective scale and largest large convective scale. Default: l = c(4, 7).
%%     ~~Describe \code{l} here~~
}
  \item{thres}{0 or a positive number. Threshold for rain rate. Default is 0.1 mm/h, because we calculate LWOI only for grid points, where rain rate is >= 0.1 mm/h.
%%     ~~Describe \code{thres} here~~
}
  \item{flat}{The number of grid points at each side, which should be smoothed lineraly. Default is 25. For quadratic arrays with dimension 2^n x 2^n boundaries are not smoothed and flat is a dummy variable.
%%     ~~Describe \code{flat} here~~
}
  \item{verbose}{If TRUE, the function prints progress statements and calculation time. Default is FALSE.
%%     ~~Describe \code{verbose} here~~
}
  \item{periodic}{If TRUE, the field x is already of size 2^n x 2^n and has perdiodic boundaries. Default is FALSE.
%%     ~~Describe \code{periodic} here~~
}
}
\value{This functions returns a list with follwing elements:
\item{WOI1orig}{Original WOI1 (but calculated with DaubExPhase4).}
\item{WOI2orig}{Original WOI2 (but calculated with DaubExPhase4).}
\item{WOI3orig}{Original WOI3 (but calculated with DaubExPhase4).}
\item{WOIorig}{Original WOI (but calculated with DaubExPhase4).}
\item{WOI1}{Modified WOI1.}
\item{WOI2}{Modified WOI2.}
\item{WOI3}{Modified WOI3.}
\item{WOI}{Modified WOI.}
\item{LWOI1}{2D array of LWOI1.}
\item{LWOI2}{2D array of LWOI2.}
\item{LWOI3}{2D array of LWOI3.}
\item{LWOI}{2D array of LWOI.}
\item{s}{Smallest and largest small convective scale.}
\item{l}{Smallest and largest large convective scale.}
\item{flat}{The number of grid points at each side, which are smoothed.}
\item{quad}{TRUE or FALSE, if the ingoing array is quadratic and 2^n x 2^n.}
\item{thres}{Threshold of LWOI calculation.}
\item{RR}{Ingoing array.}
\item{ts}{Computation time in seconds.}
}
\references{Brune, S., Kapp, F., & Friederichs, P. (2018). A wavelet-based analysis of convective organization in ICON large-eddy simulations. Quarterly Journal of the Royal Meteorological Society, 144(717), 2812-2829.
%% ~put references to the literature/web site here ~
}
\author{Sebastian Brune
%%  ~~who you are~~
}
\note{This function preprocesses the ingoing array. If x is quadratic with dimension 2^n, the function generates periodic boundaries by mirroring. If x is not 2^n x 2^n, the boundaries are smoothed regarding flat and 0 are added. In this case, the dimensions of the LWOI arrays are 2*n smaller than dim(x).
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{This function calculates WOI/LWOI only for arrays up to 2048 x 2048. Minimum size is 16 x 16.}

\examples{
# Random array of dim 350 x 300
x <- array((rnorm(350*300)), dim = c (350, 300))
s <- c(1, 2)
l <- c(3, 4)
thres <- 0.1
flat <- 25
WOIres <- WOI(x = x, s = s, l = l, thres = thres, flat = flat, 
	    verbose = TRUE)
# original WOI (Brune et al., 2018)
WOIorig <- WOIres$WOIorig
print(paste("Original WOI:", WOIorig))
# modified WOI
WOI <- WOIres$WOI 
print(paste("Modified WOI:", WOI))
# local WOI
LWOI <- WOIres$LWOI
par(mfrow = c(1, 2))
image(WOIres$RR, main = "Rain")
image(LWOI, main = "LWOI")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
