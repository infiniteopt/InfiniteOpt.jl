module MeasureEvalMethods

# Import the necessary packages.
import FastGaussQuadrature
import JuMP
import Distributions
const JuMPC = JuMP.Containers
using ..InfiniteOpt

# include jl files here
include("methods.jl")
# Export functions here
export generate_measure_data, MC_sampling, Gauss_Legendre, Gauss_Hermite,
Gauss_Laguerre, infinite_transform, support_formatting, measure_dispatch

end # end module
