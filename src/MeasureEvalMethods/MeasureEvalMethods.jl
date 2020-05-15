module MeasureEvalMethods

# Import the necessary packages.
import FastGaussQuadrature
import JuMP
import Distributions
const JuMPC = JuMP.Containers
using ..InfiniteOpt

# include jl files here
include("integrals.jl")

end # end module
