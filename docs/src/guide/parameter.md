# Infinite Parameters
A guide and manual to the definition and use of infinite parameters in
`InfiniteOpt`.

## Overview
Infinite parameters are what parameterize the infinite decision spaces for
infinite dimensional mathematical programs. In dynamic optimization this
corresponds to time and in stochastic optimization this to uncertain parameters
that follow a certain underlying statistical distribution. `InfiniteOpt`
considers principally two kinds of infinite parameters, ones defined over
continuous intervals and ones characterized by a distribution (others can be
added by defining a custom type). These can be used to parameterize
infinite variables, point variables, measures, and can be used directly inside
constraints.

## Basic Usage
First, we need to initialize and add infinite parameters to our `InfiniteModel`.
This can be accomplished using [`@infinite_parameter`](@ref). For example, let's
define a parameter for time in a time interval from 0 to 10:
```jldoctest basic
julia> using InfiniteOpt, JuMP

julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10])
t
```
Now `t` is a Julia variable that stores a [`ParameterRef`](@ref) which points to
where the time parameter is stored in `model`. We should also specify supports
that will later be used to discretized the time domain. This can be accomplished
with [`add_supports`](@ref):
```jldoctest basic
julia> add_supports(t, [0, 0.25, 0.75, 1])

julia> supports(t)
4-element Array{Number,1}:
 0.0
 0.25
 0.75
 1.0
```
Here only 4 supports are specified for the sake of example. Alternatively, we
could have specified initialized the parameter and added supports in just one
step using the `supports` keyword argument:
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10], supports = [0, 0.25, 0.75, 1])
t
```


We could also define a random parameter described by a statistical
distribution. This can be accomplished using [`@infinite_parameter`](@ref) in
combination with a distribution from
[`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl). For
example let's define a vector of random parameters described by a Normal
distribution:
```jldoctest basic
julia> using Distributions

julia> @infinite_parameter(model, xi[i = 1:3] in Normal())
3-element Array{ParameterRef,1}:
 xi[1]
 xi[2]
 xi[3]
```
Note that we could have used `i` as an index to assign a different distribution
to each parameter. Supports should also be added for each parameter as shown
above. Now we have infinite parameters `t` and `xi` that are ready to be used
in defining infinite variables and constraints. We also mention here that the
[`@infinite_parameter`](@ref) macro is designed to closely emulate
[`JuMP.@variable`](@ref) and thus handles arrays and keyword arguments in the
same way. This is described in more detail below.

## Infinite Sets
The domain of a given infinite parameter is described by an infinite set
inherited from [`AbstractInfiniteSet`](@ref). `InfiniteOpt` natively supports
two such sets. The first is [`IntervalSet`](@ref) which describes a continuous
interval from some lower bound up to some upper bound. Typically, this range
is inclusive of the bounds. Such sets often arise for parameters that pertain to
time and/or spatial position. For example, to define a position interval
``[-2, 2]`` we would use:
```jldoctest; setup = :(using InfiniteOpt)
julia> set = IntervalSet(-2, 2)
IntervalSet(-2.0, 2.0)
```
In principle infinite bounds are acceptable, but are not actively supported at
the moment.

The second kind of set is that of [`DistributionSet`](@ref) which contains
distribution from [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl)
that characterizes an infinite parameter. Univariate and multivariate
distributions can be used, but multivariate distributions are only valid when
defining an array of parameters with equivalent dimensions. For example, let's
make a `DistributionSet` that depends on a Beta distribution:
```jldoctest; setup = :(using InfiniteOpt, Distributions)
julia> set = DistributionSet(Beta(2,2))
DistributionSet{Beta{Float64}}(Distributions.Beta{Float64}(α=2.0, β=2.0))
```
User-defined distributions are also permissible so long as they are created in
accordance with `Distributions.jl`.

Furthermore, custom infinite sets that inherit `AbstractInfiniteSet` can also
be defined. See [Extensions](@ref) for more information.

## Parameter Definition
Defining/initializing an infinite parameter principally involves the following
steps:
- Define an [`AbstractInfiniteSet`](@ref)
- Define support points within the set to later discretize the parameter
- Specify if parameter is independent (only for multi-dimensional parameter groups)
- Construct an [`InfOptParameter`](@ref) to store this information
- Add the `InfOptParameter` object to an `InfiniteModel` and assign a name
- Create a [`ParameterRef`](@ref) that points to the parameter object

Infinite set definition is described above in [Infinite Sets](@ref). The
supports should be a vector of finite numbers that are drawn from the domain of
the infinite set. These supports will be used to transcribe the `InfiniteModel`
in preparation for it to be optimized. If desired, the supports can be specified
after the parameter is defined and an empty vector will be used to construct the
infinite parameter. We'll discuss independence a little further below. The
[`build_parameter`](@ref) function is used to construct the `InfOptParameter`.
For example, let's create a time parameter ``t \in [0, 10]`` with supports
`[0, 2, 5, 7, 10]`:
```jldoctest time_define; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> set = IntervalSet(0, 10)
IntervalSet(0.0, 10.0)

julia> t_param = build_parameter(error, set, supports = [0, 2, 5, 7, 10])
InfOptParameter{IntervalSet}(IntervalSet(0.0, 10.0), [0, 2, 5, 7, 10], false)
```  
Now we have a `InfOptParameter` that contains an `IntervalSet` and supports.
Note that the `num_params` and `independent` arguments are primarily meant to be
used by [`@infinite_parameter`](@ref) to define array parameters and specify if
they are correlated or independent.

Let's now add `t_param` to our `InfiniteModel` using [`add_parameter`](@ref)
and assign it the name of `t`:
```jldoctest time_define
julia> t_ref = add_parameter(model, t_param, "t")
t
```  

## Parameter Queries


## Parameter Modification


## Datatypes
```@index
Pages   = ["parameter.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
AbstractInfiniteSet
IntervalSet
DistributionSet
InfOptParameter
ParameterRef
```

## Methods/Macros
```@index
Pages   = ["parameter.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
@infinite_parameter
IntervalSet(::Number, ::Number)
build_parameter
add_parameter
used_by_constraint(::ParameterRef)
used_by_measure(::ParameterRef)
used_by_variable(::ParameterRef)
is_used(::ParameterRef)
JuMP.delete(::InfiniteModel, ::ParameterRef)
JuMP.is_valid(::InfiniteModel, ::ParameterRef)
JuMP.name(::ParameterRef)
JuMP.set_name(::ParameterRef, ::String)
num_parameters(::InfiniteModel)
infinite_set(::ParameterRef)
set_infinite_set(::ParameterRef, ::AbstractInfiniteSet)
JuMP.has_lower_bound(::ParameterRef)
JuMP.lower_bound(::ParameterRef)
JuMP.set_lower_bound(::ParameterRef, ::Number)
JuMP.has_upper_bound(::ParameterRef)
JuMP.upper_bound(::ParameterRef)
JuMP.set_upper_bound(::ParameterRef, ::Number)
num_supports(::ParameterRef)
has_supports(::ParameterRef)
supports(::ParameterRef)
supports(::AbstractArray{<:ParameterRef})
set_supports(::ParameterRef, ::Vector{<:Number})
add_supports(::ParameterRef, ::Union{Number, Vector{<:Number}})
delete_supports(::ParameterRef)
fill_in_supports!(::InfiniteModel, ::Int64, ::Int64)
fill_in_supports!(::ParameterRef, ::Int64, ::Int64)
generate_supports(::IntervalSet, ::Int64, ::Int64)::Vector{Float64}
generate_supports(::DistributionSet, ::Int64, ::Int64)::Vector{Float64}
group_id(::ParameterRef)
group_id(::AbstractArray{<:ParameterRef})
is_independent(::ParameterRef)
set_independent(::ParameterRef)
unset_independent(::ParameterRef)
parameter_by_name(::InfiniteModel, ::String)
all_parameters(::InfiniteModel)
```
