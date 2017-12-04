using LsqFit
using StatsBase

#####
# Multiple ways to estimate an autocorrelation time (tau).
#####
"""
Integrated autocorrelation time. Defined as `tau_int = sum(autocor(M))+0.5`.
"""
function tau_integrated(X::Vector{T}) where T<:Real
    return sum(autocor(X))+0.5
end

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

"""
Autocorrelation time obtained from either static or dynamic binning analysis.

Optional keyword `binsize`. Decides wether we do static (fixed, finite binsize)
or dynamic (varying binsize) binning. (Based on Ch. 3.4)
"""
function tau_binning(X::Vector{T}; binsize=0) where T<:Real

    if binsize == 0
        #dynamic binning with varying binsize
        bss, R = R_function(X)
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
        Rplateau = R_value(X, binsize)
    end

    return 0.5*(Rplateau-1)
end



#####
# Statistical error (effective standard deviation) calculated from explicit estimate of autocorrelation time
#####
## CORRECT?

"""
Statistical error for correlated data from tau_binning analysis (find plateau (default) or fixed bin size). (Eq. 3.14)
"""
tau_binning_error(X::Vector{T}; binsize=0) where T<:Real = sqrt(var(X)*(1+2*tau_binning(X, binsize=binsize))) # 1/length(X) factor "left out" to have tau_binning_error = std for uncorrelated data
tau_binning_error(X::Vector{T}; binsize=0) where T<:Complex = sqrt(tau_binning_error(real(X), binsize=binsize)^2 + tau_binning_error(imag(X), binsize=binsize)^2)
tau_binning_error(X::Array{T}; binsize=0) where T<:Number = squeeze(mapslices(ts->tau_binning_error(ts, binsize=binsize), X, ndims(X)), ndims(X))

Neff_binning(X::Vector{T}; binsize=0) where T<:Real = floor(Int, length(X)/(1+2*tau_binning(X, binsize=binsize)))

# function tau_binning_error(X::Array{T}; binsize=0) where T<:Real
#     N = size(X, ndims(X)) # length of time series
#     linX = reshape(X, (:, N))
#     Nel = size(linX, 1) # number of elements/time series
#     el_shape = size(X)[1:end-1]
    
#     errs = zeros(T, Nel)
#     for i in 1:size(errs, 1)
#         errs[i] = tau_binning_error(linX[i,:], binsize=binsize)
#     end
#     return reshape(errs, el_shape)
# end

"""
Statistical error for correlated data from integrated autocorrelation time. (Eq. 3.14)
"""
tau_integrated_error(X::Vector{T}) where T<:Real = sqrt(var(X) *(2*tau_integrated(X)))
tau_integrated_error(X::Vector{T}) where T<:Complex = sqrt(tau_integrated_error(real(X))^2 + tau_integrated_error(imag(X))^2)
tau_integrated_error(X::Array{T}) where T<:Number = squeeze(mapslices(ts->tau_integrated_error(ts), X, ndims(X)), ndims(X))


Neff_integrated(X::Vector{T}) where T<:Real = floor(Int, length(X)/(2*tau_integrated(X)))

"""
Statistical error for correlated data from fitting `C(t) ~ exp(-t/tau)`. (Eq. 3.14)
"""
tau_fitting_error(X::Vector{T}) where T<:Real = sqrt(var(X)*(1+2*tau_fitting(X)))
tau_fitting_error(X::Vector{T}) where T<:Complex = sqrt(tau_fitting_error(real(X))^2 + tau_fitting_error(imag(X))^2)
tau_fitting_error(X::Array{T}) where T<:Number = squeeze(mapslices(ts->tau_fitting_error(ts), X, ndims(X)), ndims(X))

Neff_fitting(X::Vector{T}) where T<:Real = floor(Int, length(X)/(1+2*tau_fitting(X)))