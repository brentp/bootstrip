\name{bootstrip.iter}
\alias{bootstrip.iter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrip.iter(mat, mod, fn, iterations = c(200, 5000, 1e+05), p_samples = 0.5, mc.cores = 12, smooth.sd = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
}
  \item{mod}{
%%     ~~Describe \code{mod} here~~
}
  \item{fn}{
%%     ~~Describe \code{fn} here~~
}
  \item{iterations}{
%%     ~~Describe \code{iterations} here~~
}
  \item{p_samples}{
%%     ~~Describe \code{p_samples} here~~
}
  \item{mc.cores}{
%%     ~~Describe \code{mc.cores} here~~
}
  \item{smooth.sd}{
%%     ~~Describe \code{smooth.sd} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mat, mod, fn, iterations = c(200, 5000, 1e+05), p_samples = 0.5, 
    mc.cores = 12, smooth.sd = 0) 
{
    lims = NA
    subset = rep(TRUE, nrow(mat))
    for (it in iterations) {
        message(paste("sampling", sum(subset), "rows", it, "times"))
        beta = bootstrip(mat[subset, ], mod, fn, it, p_samples, 
            mc.cores)
        if (any(is.na(lims))) {
            lims = bootstrip.limits(beta)
        }
        else {
            lims[subset, ] = bootstrip.limits(beta)
        }
        lims[subset, "samples"] = it
        subset = (lims$pvalue <= (2/it))
        if (sum(subset) < 2) {
            break
        }
    }
    lims$beta.orig = fn(lmFit(mat, mod))
    lims
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
