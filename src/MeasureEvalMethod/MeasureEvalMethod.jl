module MeasureEvalMethod

# Import the necessary packages.
import FastGaussQuadrature
using ..InfiniteOpt

# include jl files here
include("methods.jl")
# Export functions here
export generate_measure_data, MC_sampling, Gauss_Legendre, Gauss_Hermite,
Gauss_Laguerre, infinite_transform

end # end module
