\name{permute.residuals}
\alias{permute.residuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
permute.residuals(mat, mod, mod0, iterations = 100, p_samples = 1, mc.cores = 12)
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
function (mat, mod, mod0, iterations = 100, p_samples = 1, mc.cores = 12) 
{
    stopifnot(nrow(mod) == ncol(mat))
    reduced_lm = lmFit(mat, mod0)
    reduced_residuals = residuals(reduced_lm, mat)
    reduced_fitted = fitted(reduced_lm)
    fit = lmFit(mat, mod)
    size = p_samples * nrow(mod)
    coef.name = setdiff(colnames(mod), colnames(mod0))
    beta.orig = coefficients(fit)[, coef.name]
    rm(reduced_lm, fit)
    gc()
    nc = ncol(reduced_residuals)
    beta.list = mclapply(1:iterations, function(ix) {
        if (p_samples < 1) {
            sub_ids = sample.int(nc, size = size)
        }
        else {
            sub_ids = 1:nc
        }
        mat_sim = reduced_fitted[, sub_ids] + reduced_residuals[, 
            sample(sub_ids)]
        coefficients(lmFit(mat_sim, mod[sub_ids, ]))[, coef.name]
    }, mc.cores = mc.cores)
    beta = matrix(NA, nrow(mat), ncol = iterations)
    for (i in 1:iterations) {
        beta[, i] = beta.list[[i]]
    }
    df = data.frame(pvalue = unlist(lapply(1:nrow(beta), function(i) {
        row = beta[i, ]
        val = ecdf(row)(beta.orig[i])
        max(min(val, 1 - val) * 2, 1/length(row))
    })), beta.orig = beta.orig)
    rownames(df) = rownames(mat)
    df
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
