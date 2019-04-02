\name{buildperiodic}
\alias{buildperiodic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Builds periodic array by mirroring at side and top
%%  ~~function to do ... ~~
}
\description{This function generates periodic boundaries by mirroring at side and top. The resulting array is four times larger than the incoming array.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
buildperiodic(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{2D quadratic array.
%%     ~~Describe \code{x} here~~
}
}
\value{Returns a 2D array, which is four times larger than the incoming array.
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
x <- matrix(1:12, nrow = 4)
out <- buildperiodic(x = x)
out
}
