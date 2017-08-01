"""
**ErrorAnalysis** is a collection of tools to calculate errors for Monte Carlo data.
"""
module ErrorAnalysis

using StatsBase
using LsqFit
using PyPlot


#####
# Statistical error (effective standard deviation)
#####

"""
Statistical error for correlated data from binning analysis.
"""
error_binning(X::Vector{T}; binsize=0) where T<:Real = sqrt(var(X)*(1+2*tau_binning(X, binsize=binsize)))

error_binning(X::Vector{T}; binsize=0) where T<:Complex = abs(error_binning(real(X), binsize=binsize) + error_binning(imag(X), binsize=binsize))
# error_binning(X::Vector{T}; binsize=0) where T<:Complex = error_binning(real(X), binsize=binsize) + im*error_binning(imag(X), binsize=binsize)

error_binning(X::Array{T}; binsize=0) where T<:Number = squeeze(mapslices(ts->error_binning(ts, binsize=binsize), X, ndims(X)), ndims(X))
export error_binning

Neff_binning(X::Vector{T}; binsize=0) where T<:Real = floor(Int, length(X)/(1+2*tau_binning(X, binsize=binsize)))
export Neff_binning

# function error_binning(X::Array{T}; binsize=0) where T<:Real
#     N = size(X, ndims(X)) # length of time series
#     linX = reshape(X, (:, N))
#     Nel = size(linX, 1) # number of elements/time series
#     el_shape = size(X)[1:end-1]
    
#     errs = zeros(T, Nel)
#     for i in 1:size(errs, 1)
#         errs[i] = error_binning(linX[i,:], binsize=binsize)
#     end
#     return reshape(errs, el_shape)
# end

"""
Statistical error for correlated data from integrated autocorrelation.
"""
error_integrated(X::Vector{T}) where T<:Real = sqrt(var(X)*(2*tau_integrated(X)))
error_integrated(X::Vector{T}) where T<:Complex = abs(error_integrated(real(X)) + error_integrated(imag(X)))
error_integrated(X::Array{T}) where T<:Number = squeeze(mapslices(ts->error_integrated(ts), X, ndims(X)), ndims(X))
export error_integrated

Neff_integrated(X::Vector{T}) where T<:Real = floor(Int, length(X)/(2*tau_integrated(X)))
export Neff_integrated

"""
Statistical error for correlated data from fitting `C(t) ~ exp(-t/tau)`.
"""
error_fitting(X::Vector{T}) where T<:Real = sqrt(var(X)*(1+2*tau_fitting(X)))
error_fitting(X::Vector{T}) where T<:Complex = abs(error_fitting(real(X)) + error_fitting(imag(X)))
error_fitting(X::Array{T}) where T<:Number = squeeze(mapslices(ts->error_fitting(ts), X, ndims(X)), ndims(X))
export error_fitting

Neff_fitting(X::Vector{T}) where T<:Real = floor(Int, length(X)/(1+2*tau_fitting(X)))
export Neff_fitting


#####
# Multiple ways to calculate an autocorrelation time (tau).
#####
"""
Integrated autocorrelation time. Defined as `tau_int = sum(autocor(M))+0.5`.
"""
function tau_integrated(X::Vector{T}) where T<:Real
    return sum(autocor(X))+0.5
end
export tau_integrated


"""
Autocorrelation time obtained from fitting autocorrelation C(t) to `exp(-t/tau_fitting)`.
"""
function tau_fitting(X::Vector{T}) where T<:Real
    logautocor = T[]
    Cx = autocor(X)

    firstneg = findfirst(x->x<=0, Cx)
    logCx = firstneg != 0 ? log.(Cx[1:firstneg-1]) : log.(Cx)
    length(logCx) <= 2 ? warn("Fitting log autocorrelation with less than 3 datapoints.") : nothing

    # Fitting
    model(x,p) = p[1].*x+p[2]
    pstart = [-1., 1.]
    fitres = curve_fit(model, 1.0:1.0:length(logCx), logCx, pstart)

    return -1/fitres.param[1]
end
export tau_fitting


"""
Autocorrelation time obtained from either static or dynamic binning analysis.

Optional keyword `binsize`. Decides wether we do static (fixed, finite binsize)
or dynamic (varying binsize) binning.
"""
function tau_binning(X::Vector{T}; binsize=0) where T<:Real

    if binsize == 0
        #dynamic binning with varying binsize
        bss, R = binning(X)
        length(R) > 0 || error("Binning failed. Maybe min_nbins to large?")
        Rplateau = R[end]
        for i in 1:length(bss)-2
            width = abs(maximum(R[i:i+2]) - minimum(R[i:i+2]))
            if width <= 0.05
                Rplateau = R[i]
                break;
            end
        end
        Rplateau != R[end] || warn("Couldn't find plateau. Maybe not converged?")
    else
        #static binning with fixed, single binsize
        Rplateau = binning(X, binsize)
    end

    return 1/2*(Rplateau-1)
end
export tau_binning




#####
# Static and Dynamic Binning
#####
"""
Performs DYNAMIC binning analysis, i.e. groups datapoints in bins of varying size `bs`.
Returns the used binsizes `bss` and the characteristic function `R(bss)` which should feature
a plateau, i.e. `R(bs_p) ~ R(bs)` for `bs >= bs_p`.

Optional keyword `min_nbins`. Only bin sizes used that lead to at least `min_nbins` bins.
Optional keyword `Rplot`. Creates a plot R vs bin size.
"""
function binning(X::Vector{T}; min_nbins=500) where T<:Real
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
Performs STATIC binning analysis, i.e. groups datapoints in bins of fixed binsize.
Returns characteristic coefficient R.
"""
function binning(X::Vector{T}, binsize::Int) where T<:Real
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
export binning


"""
Plots the binning factor R for multiple bin sizes.
Ideally, this plot shows a plateau.

Returns bss and R from binning(X).
"""
function binning_plot(X::Vector{T}; min_nbins=500, figsize=(6,4)) where T<:Real
    bss, R = binning(X, min_nbins=min_nbins)
    figure(figsize=figsize)
    plot(bss, R, "m.-")
    xlabel("bin size")
    ylabel("R");

    return bss, R
end
export binning_plot

"""
Plot histogram of X with mean and statistical errorbars (including correlation effects) and autocorrelation function C(t).

Optional keyword `error` can be "binning", "integrated", or "fitting".
"""
function error_plot(X::Vector{T}; binsize=0, histbins=50, error="binning", figsize=(10,4), digits=3) where T<:Real
    fig, ax = subplots(1, 2, figsize=figsize)

    Xmean = mean(X)
    Xstd = std(X)

    if error == "binning"
        err = error_binning(X, binsize=binsize)
    elseif error == "integrated"
        err = error_integrated(X)
    elseif error == "fitting"
        err = error_fitting(X)
    end

    ax[1][:hist](X, histbins, color="gray", alpha=.5, normed=1)
    ax[1][:set_ylabel]("\$ P \$")
    ax[1][:set_xlabel]("\$ X \$")
    ax[1][:set_yticks]([])
    ax[1][:axvline](Xmean, color="black", label="\$ $(round(Xmean, digits)) \$", linewidth=2.0)
    
    ax[1][:axvline](Xmean+err, color="r", label="\$ \\pm $(round(err, digits)) \$", linewidth=2.0)
    ax[1][:axvline](Xmean-err, color="r", linewidth=2.0)

    ax[1][:axvline](Xmean+Xstd, color="g", label="\$ \\pm $(round(Xstd, digits)) \$ (std)", linewidth=2.0)
    ax[1][:axvline](Xmean-Xstd, color="g", linewidth=2.0)
    
    ax[1][:legend](frameon=false, loc="best")
    
    ax[2][:plot](autocor(X), "-", color="k", linewidth=2.0)
    ax[2][:set_xlabel]("\$ t \$")
    ax[2][:set_ylabel]("\$ C(t) \$")

    tight_layout()
    return ax
end
export error_plot


include("Jackknife.jl")


end # module
