"""
**ErrorAnalysis** is a collection of tools to calculate errors for Monte Carlo data.
"""
module ErrorAnalysis

using PyPlot

include("binning.jl")
include("Jackknife.jl")
include("autocorrtimes.jl") # old


# Exports
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