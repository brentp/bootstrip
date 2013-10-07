options(stringsAsFactors=FALSE)
library(limma)
library(parallel)

# fn returns the statistic of interest from a limma fit object.
bootstrip = function(mat, mod, fn, iterations=100, p_samples=0.5, mc.cores=12){
    stopifnot(nrow(mod) == ncol(mat))

    ta = mclapply(1:iterations, function(core_i){ 
            idx = sample.int(ncol(mat), p_samples * ncol(mat), replace=TRUE)
            msub = mat[,idx,drop=TRUE]
            fit = lmFit(msub, mod[idx,,drop=TRUE])
            fn(fit)
        }, mc.cores=mc.cores)
    beta = matrix(NA, nrow(mat), ncol=iterations)
    for(i in 1:iterations){
        beta[,i] = ta[[i]]
    }
    colnames(beta) = paste0("sample_", 1:iterations)
    rownames(beta) = rownames(mat)
    beta
}

bootstrip.iter = function(mat, mod, fn, iterations=c(200, 5000, 1e5), p_samples=0.5, mc.cores=12, smooth.sd=0){
    lims = NA
    subset = rep(TRUE, nrow(mat))
    for(it in iterations){
        message(paste("sampling", sum(subset), "rows", it, "times"))
        beta = bootstrip(mat[subset,], mod, fn, it, p_samples, mc.cores)
        if(any(is.na(lims))){
            lims = bootstrip.limits(beta)
        } else {

            lims[subset,] = bootstrip.limits(beta)
        }
        lims[subset, "samples"] = it
        # repeat on this subset. with a larger number of samples
        subset = (lims$pvalue <= (2 / it))
        if(sum(subset) < 2){ break }
    }
    lims$beta.orig = fn(lmFit(mat, mod))
    lims
}

bootstrip.limits = function(beta, probs=c(0.025, 0.975)){
    # pvalue is 2-tailed probability that the observed
    # distribution overlaps 0.
    pvalue = apply(beta, 1, function(row){
        val=ecdf(row)(0); max(min(val, 1-val) * 2, 1/length(row)) 
    })
    df = data.frame(pvalue)
    df$beta.mean = apply(beta, 1, mean)

    df$beta.median = apply(beta, 1, median)
    for(p in probs){
        df[,paste0("beta.pct_", gsub(".", "p", p, fixed=T))] = 
                      apply(beta, 1, function(row){ quantile(row, probs=p) })
    }
    df
}

permute.residuals = function(mat, mod, mod0, iterations=100, p_samples=1, mc.cores=12){
    stopifnot(nrow(mod) == ncol(mat))

    reduced_lm = lmFit(mat, mod0)
    reduced_residuals = residuals(reduced_lm, mat)
    reduced_fitted = fitted(reduced_lm)

    fit = lmFit(mat, mod)
    size = p_samples * nrow(mod)
    
    coef.name = setdiff(colnames(mod), colnames(mod0))
    beta.orig = coefficients(fit)[,coef.name]

    rm(reduced_lm, fit); gc() 
    nc = ncol(reduced_residuals)

    beta.list = mclapply(1:iterations, function(ix){
        # TODO: make sure cols from fit match mod when taking subset
        if( p_samples < 1){
            sub_ids = sample.int(nc, size=size)
        } else {
            sub_ids = 1:nc
        }
        # take the original model and fitted data in same order, but
        # shuffle residuals.
        mat_sim = reduced_fitted[, sub_ids] + reduced_residuals[,sample(sub_ids)]
        coefficients(lmFit(mat_sim, mod[sub_ids,]))[,coef.name]
    }, mc.cores=mc.cores)

    beta = matrix(NA, nrow(mat), ncol=iterations)
    for(i in 1:iterations){
        beta[,i] = beta.list[[i]]
    }

    df = data.frame(
       pvalue = unlist(lapply(1:nrow(beta), function(i){
        row = beta[i,]
        val=ecdf(row)(beta.orig[i]);
        max(min(val, 1-val) * 2, 1/length(row)) 
    })), beta.orig = beta.orig)
    rownames(df) = rownames(mat)
    df
}

permute.residuals.iter = function(mat, mod, mod0, iterations=c(200, 5000, 1e5, 2e6),
                                  p_samples=1, mc.cores=12){
    stopifnot(nrow(mod) == ncol(mat))
    subset = rep(TRUE, nrow(mat))
    df = NA
    for(it in iterations){
        message(paste("sampling", sum(subset), "rows", it, "times"))
        if(any(is.na(df))){
            df = permute.residuals(mat, mod, mod0, it, p_samples, mc.cores)
        } else {
            # TODO: only sending in the most significant to permute. will bias results.
            df[subset,] = permute.residuals(mat[subset,], mod, mod0, it, p_samples, mc.cores)
        }
        df[subset, "samples"] = it
        # repeat on this subset. with a larger number of samples
        subset = (df$pvalue <= (2 / it))
        if(sum(subset) == 0){ break }
    }
    df
}

#library(multtest)

#mult.permute = function(mat, mod, coef){
#
#    MTP(mat, Z=mod, test="lm.XvsZ", Z.incl=1:ncol(mod), Z.test=coef,
#        get.adjp=TRUE, get.cr=TRUE)
