module MeasureToolbox # TODO make this an independent package

# Import the necessary packages.
import FastGaussQuadrature
import JuMP
import Distributions
const JuMPC = JuMP.Containers
using ..InfiniteOpt

# include jl files here
include("integrals.jl")
include("expectations.jl")
include("support_sums.jl")

export Automatic, UniTrapezoid, UniMCSampling, UniIndepMCSampling, Quadrature,
GaussRadau, GaussHermite, GaussLegendre, GaussRadau, GaussLobatto, FEGaussLobatto,
GaussJacobi, GaussChebyshev, GaussLaguerre, MultiMCSampling, MultiIndepMCSampling, 
uni_integral_defaults, set_uni_integral_defaults, clear_uni_integral_defaults,
integral, multi_integral_defaults, set_multi_integral_defaults, 
clear_multi_integral_defaults, expect, support_sum,
generate_integral_data, 𝔼, ∫, InternalGaussLobatto

end # end module
