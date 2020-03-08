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
export generate_measure_data, mc_sampling, gauss_legendre, gauss_hermite,
gauss_laguerre, infinite_transform, generate_supports_and_coeffs,
get_method_registry, set_method_registry

# Export constants
export default_set_types, default_methods

end # end module
