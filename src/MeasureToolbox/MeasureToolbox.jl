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

end # end module
