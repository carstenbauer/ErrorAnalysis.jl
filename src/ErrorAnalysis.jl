"""
**ErrorAnalysis** is a collection of tools to calculate errors for Monte Carlo data.
"""
module ErrorAnalysis

include("binning.jl")
include("Jackknife.jl")
include("autocorrtimes.jl") # old

"""
    iswithinerrorbars(a, b, δ[, print=false])

Checks whether numbers `a` and `b` are equal up to given error `δ`.
Will print `x ≈ y + k·δ` for `print=true`.

Is equivalent to `isapprox(a,b,atol=δ,rtol=zero(b))`.
"""
function iswithinerrorbars(a::T, b::S, δ::Real, print::Bool=false) where T<:Number where S<:Number
  equal = isapprox(a,b,atol=δ,rtol=zero(T))
  if print && !equal
    out = a>b ? abs(a-(b+δ))/δ : -abs(a-(b-δ))/δ
    println("x ≈ y + ",round(out,4),"·δ")
  end
  return equal
end

"""
    iswithinerrorbars(A::AbstractArray{T<:Number}, B::AbstractArray{T<:Number}, Δ::AbstractArray{<:Real}[, print=false])

Elementwise check whether `A` and `B` are equal up to given real error matrix `Δ`.
Will print `A ≈ B + K.*Δ` for `print=true`.
"""
function iswithinerrorbars(A::AbstractArray{T}, B::AbstractArray{S},
                           Δ::AbstractArray{<:Real}, print::Bool=false) where T<:Number where S<:Number
  size(A) == size(B) == size(Δ) || error("A, B and Δ must have same size.")

  R = iswithinerrorbars.(A,B,Δ,false)
  allequal = all(R)

  if print && !all(R) && T<:Real && S<:Real
    O = similar(A, promote_type(T,S))
    for i in eachindex(O)
      a = A[i]; b = B[i]; δ = Δ[i]
      O[i] = R[i] ? 0.0 : round(a>b ? abs(a-(b+δ))/δ : -abs(a-(b-δ))/δ,4)
    end
    println("A ≈ B + K.*Δ, where K is:")
    display(O)
  end

  return allequal
end


# Exports
export iswithinerrorbars

export binning_error
export jackknife_error

export plot_binning_R
export plot_error

# export tau_binning
# export tau_integrated
# export tau_fitting

# export Neff_binning
# export Neff_integrated
# export Neff_fitting

# export tau_binning_error
# export tau_integrated_error
# export tau_fitting_error


end # module