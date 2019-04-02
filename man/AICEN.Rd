\name{AICEN}
\alias{AICEN}
\docType{data}
\title{File with Inverse A matrix and centre of mass matrix
%%   ~~ data name/kind ... ~~
}
\description{This file includes the lists AIx and CENx for x = 16, 32, 64, 128, 256, 512, 1024, 2048 or 4096. AIx is the inverse matrix calculated by Eckley et al. (2010). CENx includes the centre of mass for DaubExPhase4, calculated with help of the LS2W package by Eckley and Nason (2011).
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(AICEN)}
\format{List of 18 elements (9 AIx and 9 CENx). CENx include lists for x and y direction.
}
\author{Sebastian Buschow, Sebastian Brune
}
\references{
Eckley, I.A., Nason, G.P. and Treloar, R.L. (2010) Locally stationary wavelet fields with application to the modelling and analysis of image texture. Journal of the Royal Statistical Society (Series C), 59, 595-616.

Eckley, I.A. and Nason, G.P. (2011). LS2W: Implementing the Locally Stationary 2D Wavelet Process Approach in R, Journal of Statistical Software, 43(3), 1-23.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(AICEN)
image(AICEN$AI256)
str(AICEN)
}
\keyword{datasets}
