\name{bootstrip}
\alias{bootstrip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrip(mat, mod, fn, iterations = 100, p_samples = 0.5, mc.cores = 12)
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
function (mat, mod, fn, iterations = 100, p_samples = 0.5, mc.cores = 12) 
{
    stopifnot(nrow(mod) == ncol(mat))
    ta = mclapply(1:iterations, function(core_i) {
        idx = sample.int(ncol(mat), p_samples * ncol(mat), replace = TRUE)
        msub = mat[, idx, drop = TRUE]
        fit = lmFit(msub, mod[idx, , drop = TRUE])
        fn(fit)
    }, mc.cores = mc.cores)
    beta = matrix(NA, nrow(mat), ncol = iterations)
    for (i in 1:iterations) {
        beta[, i] = ta[[i]]
    }
    colnames(beta) = paste0("sample_", 1:iterations)
    rownames(beta) = rownames(mat)
    beta
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
