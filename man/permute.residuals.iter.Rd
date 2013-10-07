\name{permute.residuals.iter}
\alias{permute.residuals.iter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
permute.residuals.iter(mat, mod, mod0, iterations = c(200, 5000, 1e+05, 2e+06), p_samples = 1, mc.cores = 12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
}
  \item{mod}{
%%     ~~Describe \code{mod} here~~
}
  \item{mod0}{
%%     ~~Describe \code{mod0} here~~
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
function (mat, mod, mod0, iterations = c(200, 5000, 1e+05, 2e+06), 
    p_samples = 1, mc.cores = 12) 
{
    stopifnot(nrow(mod) == ncol(mat))
    subset = rep(TRUE, nrow(mat))
    df = NA
    for (it in iterations) {
        message(paste("sampling", sum(subset), "rows", it, "times"))
        if (any(is.na(df))) {
            df = permute.residuals(mat, mod, mod0, it, p_samples, 
                mc.cores)
        }
        else {
            df[subset, ] = permute.residuals(mat[subset, ], mod, 
                mod0, it, p_samples, mc.cores)
        }
        df[subset, "samples"] = it
        subset = (df$pvalue <= (2/it))
        if (sum(subset) == 0) {
            break
        }
    }
    df
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
