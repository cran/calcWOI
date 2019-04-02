\name{flatten}
\alias{flatten}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reduces the boundary gradients
%%  ~~function to do ... ~~
}
\description{This function smoothes the boundaries with a linear filter.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
flatten(x, filter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{2D array, which boundaries should be smoothed.
%%     ~~Describe \code{x} here~~
}
  \item{filter}{The smoothing vector with increasing elements from 0...1. The length of the vector corresponds to the number of smoothed points at each side.
%%     ~~Describe \code{filter} here~~
}
}
\value{Returns the incoming field with smoothed boundaries.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Florian Kapp
%%  ~~who you are~~
}
\examples{
# Smooth outer 25 grid points
x <- array(10, dim = c(100, 200))
xflat <- flatten(x = x, filter = seq(0, 1, , 25))
par(mfrow = c(1, 2))
image(x, main = "Original")
image(xflat, main = "Smoothed Bound")
}
