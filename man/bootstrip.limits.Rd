\name{bootstrip.limits}
\alias{bootstrip.limits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrip.limits(beta, probs = c(0.025, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{probs}{
%%     ~~Describe \code{probs} here~~
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
function (beta, probs = c(0.025, 0.975)) 
{
    pvalue = apply(beta, 1, function(row) {
        val = ecdf(row)(0)
        max(min(val, 1 - val) * 2, 1/length(row))
    })
    df = data.frame(pvalue)
    df$beta.mean = apply(beta, 1, mean)
    df$beta.median = apply(beta, 1, median)
    for (p in probs) {
        df[, paste0("beta.pct_", gsub(".", "p", p, fixed = T))] = apply(beta, 
            1, function(row) {
                quantile(row, probs = p)
            })
    }
    df
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
