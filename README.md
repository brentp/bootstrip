calculate p-value for a particular beta (or any measure), using a bootstrap-like measure. allows subsampling where as bootstrap is sampling with replacement. allows shuffling of residuals for linear model. the .iter version of each function will perform 200 iterations on the full dataset, then 5000 iterations on those with a simulated p-value < 2 / 200, then 1e5 iteration on probes with simulated p-value < 2 / 1e5 ... etc. So that we get higher precision at lower p-values.