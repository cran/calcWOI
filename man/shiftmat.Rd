\name{shiftmat}
\alias{shiftmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Shifts the elements of a matrix
%%  ~~function to do ... ~~
}
\description{This function shifts the elements of an array to the right and the top. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
shiftmat(x, dx = 0, dy = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{2D array.
%%     ~~Describe \code{x} here~~
}
  \item{dx}{Integer number. Number of grid points to shift the array to the north. Should be smaller than dim(x). Default is 0.
%%     ~~Describe \code{dx} here~~
}
  \item{dy}{Integer number. Number of grid points to shift the array to the west. Should be smaller than dim(x). Default is 0.
%%     ~~Describe \code{dy} here~~
}
}
\value{Returns an array with shifted elements.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Sebastian Brune
%%  ~~who you are~~
}
\examples{
# shift the matrix dx = 1 and dy = 2 grid points
x <- array(1:48, dim = c(6, 8))
xshift <- shiftmat(x = x, dx = 1, dy = 2) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
