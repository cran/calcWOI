\name{blowup}
\alias{blowup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adds zeros around 2D array
%%  ~~function to do ... ~~
}
\description{This function adds zeros around a 2D array to get a M x M field.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
blowup(x, M, number = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{2D array, which should be blown up to a M x M array.
%%     ~~Describe \code{x} here~~
}
  \item{M}{Dimension of new array including zeros at boundaries. M should be larger than max(dim(x)).
%%     ~~Describe \code{M} here~~
}
  \item{number}{Number, that should be added to the incoming array. Default is 0.
%%     ~~Describe \code{number} here~~
}
}
\value{Returns the M x M array with the original field in the centre.
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
# Add zeros around a 3x3 matrix to generate a 8x8 matrix.
x <- matrix(1:9, nrow = 3)
new <- blowup(x = x, M = 8, number = 0)
new
}
