\name{wavtra}
\alias{wavtra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs the wavelet transform
%%  ~~function to do ... ~~
}
\description{This function uses parts the wavelet transform of Eckely et al. (2010). We use the DaubExPhase4 wavelet for all calculations. The inverse A matrix is loaded from constants.rda. The resulting value of each transform is written to the centre of mass of the spectrum, which is also saved in constants.rda due to computation time.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wavtra(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{2D array of dimension 2^n x 2^n with n = 4, 5, ...  or 12. Periodic boundaries are assumed.
%%     ~~Describe \code{x} here~~
}
}
\value{Returns a 3D array with 2^n x 2^n x 3*n wavelet coefficients. The third dimension includes the wavelet coefficients of North-South scales 1-n, East-West scales 1-n and Diagonal scales 1-n.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Eckley, I.A., Nason, G.P. and Treloar, R.L. (2010) Locally stationary wavelet fields with application to the modelling and analysis of image texture. Journal of the Royal Statistical Society (Series C), 59, 595-616.

Eckley, I.A. and Nason, G.P. (2011). LS2W: Implementing the Locally Stationary 2D Wavelet Process Approach in R, Journal of Statistical Software, 43(3), 1-23.
%% ~put references to the literature/web site here ~
}
\author{Sebastian Brune, Sebastian Buschow
%%  ~~who you are~~
}
\examples{
x <- array(1:(2^12), dim = c(2^6, 2^6))
print(dim(x))
}
