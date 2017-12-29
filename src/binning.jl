#####
# Binning error analysis
#####

"""
    binning_error(X[; binsize=0, warnings=false])

Calculates statistical error (eff. standard deviation) for correlated data.
How: Binning of data and assuming statistical independence of bins
(i.e. R plateau has been reached). (Eq. 3.18 basically)

The default `binsize=0` indicates automatic binning.
"""
function binning_error(X::AbstractVector{T}; binsize=0, warnings=false) where T<:Real
    # Data: real numbers
    if binsize == 0
        binsize = floor(Int, length(X)/32)
        binsize = binsize==0?1:binsize
    end

    isinteger(length(X) / binsize) || !warnings ||
        warn("Non-integer number of bins $(length(X) / binsize). " *
             "Last bin will be smaller than all others.")
    
    bin_means = map(mean, Iterators.partition(X, binsize))
    return sqrt(1/length(bin_means) * var(bin_means))
end
function binning_error(X::AbstractVector{T}; binsize=0, warnings=false) where T<:Complex
    # Data: complex numbers
    sqrt(binning_error(real(X), binsize=binsize, warnings=warnings)^2 +
        binning_error(imag(X), binsize=binsize, warnings=warnings)^2)
end

# Data: arrays
function binning_error(X::AbstractArray{T}; binsize=0, warnings=false) where T<:Number
    ndimsX = ndims(X)
    mapslices(y->binning_error(y; binsize=binsize, warnings=warnings), X, ndimsX)[(Colon() for _ in 1:ndimsX-1)...,1]
end
function binning_error(X::AbstractVector{T}; binsize=0, warnings=false) where T<:(AbstractArray{S} where S)
    binning_error(cat(3, X...); binsize=binsize, warnings=warnings)
end


#####
# Useful plots
#####
"""
Plots the binning error coefficient function `R(binsize)` for multiple bin sizes.
Ideally, this plot shows a plateau (Fig. 3.3), indicating that the bin averages have become independent.
For correlated data `R>≈1` and `sqrt(R)` quantifies how much one would have underestimated
the one-sigma errorbar.

Returns bss and R from `R_function(X)`.
"""
function plot_binning_R(X::Vector{T}; min_nbins=32, figsize=(4,3)) where T<:Real

    try PyPlot
    catch
        warn("PyPlot not loaded. Loading it now...")
        eval(Expr(:toplevel, Expr(:using, Symbol("PyPlot"))))
    end

    bss, R = R_function(X, min_nbins=min_nbins)
    figure(figsize=figsize)
    plot(bss, R, "m.-")
    xlabel("bin size")
    ylabel("R");

    return bss, R
end

"""
Plot histogram of X with mean and statistical errorbars (including correlation effects) 
and autocorrelation function C(t).

Optional keyword `error` can be "tau_binning", "tau_integrated", or "tau_fitting", "binning",
and sets method for error estimation.
"""
function plot_error(X::Vector{T}; binsize=0, histbins=50, error="binning", figsize=(10,4), digits=3) where T<:Real

    try PyPlot
    catch
        warn("PyPlot not loaded. Loading it now...")
        eval(Expr(:toplevel, Expr(:using, Symbol("PyPlot"))))
    end

    fig, ax = subplots(1, 2, figsize=figsize)

    Xmean = mean(X)
    Xstd = std(X)

    if error == "binning"
        err = binning_error(X, binsize=binsize, warnings=false)
    elseif error == "tau_binning"
        err = tau_binning_error(X, binsize=binsize)
    elseif error == "tau_integrated"
        err = tau_integrated_error(X)
    elseif error == "tau_fitting"
        err = tau_fitting_error(X)
    end

    ax[1][:hist](X, histbins, color="gray", alpha=.5, normed=1)
    ax[1][:set_ylabel]("\$ P \$")
    ax[1][:set_xlabel]("\$ X \$")
    ax[1][:set_yticks]([])
    ax[1][:axvline](Xmean, color="black", label="\$ $(round.(Xmean, digits)) \$", linewidth=2.0)
    
    ax[1][:axvline](Xmean+err, color="r", label="\$ \\pm $(round.(err, digits)) \$ (binning σ)", linewidth=2.0)
    ax[1][:axvline](Xmean-err, color="r", linewidth=2.0)

    ax[1][:axvline](Xmean+2*err, color="r", alpha=.3, label="\$ \\pm $(round.(2*err, digits)) \$ (binning 2σ)", linewidth=2.0)
    ax[1][:axvline](Xmean-2*err, color="r", alpha=.3, linewidth=2.0)

    ax[1][:axvline](Xmean+Xstd, color="g", label="\$ \\pm $(round.(Xstd, digits)) \$ (std)", linewidth=2.0)
    ax[1][:axvline](Xmean-Xstd, color="g", linewidth=2.0)
    
    ax[1][:legend](frameon=false, loc="best")
    
    ax[2][:plot](autocor(X), "-", color="k", linewidth=2.0)
    ax[2][:set_xlabel]("\$ t \$")
    ax[2][:set_ylabel]("\$ C(t) \$")

    tight_layout()
    return ax
end

#####
# Calculation of error coefficient (function) R. (Ch. 3.4 in QMC book)
#####
"""
Groups datapoints in bins of varying size `bs`.
Returns the used binsizes `bss` and the error coefficient function `R(bss)` (Eq. 3.20) which should feature
a plateau, i.e. `R(bs_p) ~ R(bs)` for `bs >= bs_p`. (Fig. 3.3)

Optional keyword `min_nbins`. Only bin sizes used that lead to at least `min_nbins` bins.
"""
function R_function(X::Vector{T}; min_nbins=10) where T<:Real
    bss = Int[]
    N = length(X)
    Xmean = mean(X)
    Xvar = var(X)

    for bs in 1:N
        N%bs == 0 ? push!(bss, bs) : nothing
    end

    n_bins = Int.(N./bss)
    bss = bss[n_bins .>= min_nbins] # at least 500 bins for every binsize

    R = zeros(length(bss))
    for (i, bs) in enumerate(bss)

        blockmeans = vec(mean(reshape(X, (bs,n_bins[i])), 1))

        blocksigma2 = 1/(n_bins[i]-1)*sum((blockmeans - Xmean).^2)

        R[i] = bs * blocksigma2 / Xvar
    end

    return bss, R
end

"""
Groups datapoints in bins of fixed binsize and returns error coefficient R. (Eq. 3.20)
"""
function R_value(X::Vector{T}, binsize::Int) where T<:Real
    N = length(X)
    n_bins = div(N,binsize)
    lastbs = rem(N,binsize)
    lastbs == 0 || (lastbs >= binsize/2 || warn("Binsize leads to last bin having less than binsize/2 elements."))


    blockmeans = vec(mean(reshape(X[1:n_bins*binsize], (binsize,n_bins)), 1))
    if lastbs != 0 
        vcat(blockmeans, mean(X[n_bins*binsize+1:end]))
        n_bins += 1
    end

    blocksigma2 = 1/(n_bins-1)*sum((blockmeans - mean(X)).^2)
    return binsize * blocksigma2 / var(X)
end