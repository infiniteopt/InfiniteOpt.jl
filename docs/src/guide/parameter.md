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
could have initialized the parameter and added supports in just one step using
the `supports` keyword argument:
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

### Manual Definition
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
We can also create an anonymous infinite parameter by dropping the name from
the `add_parameter` function call. For example:
```jldoctest time_define
julia> t_ref_noname = add_parameter(model, t_param)
noname
```
Now suppose we want to create an infinite parameter that is a random variable
with a given distribution. We follow the same procedure as above, except we use
distributions from `Distributions.jl` to define [`DistributionSet`](@ref). For
example, let's consider a random variable ``x \in \mathcal{N}(0,1)`` with
supports `[-0.5, 0.5]`:
```jldoctest rand_define; setup = :(using InfiniteOpt, Distributions; model = InfiniteModel())
julia> dist = Normal(0., 1.)
Normal{Float64}(μ=0.0, σ=1.0)

julia> set = DistributionSet(dist)
DistributionSet{Normal{Float64}}(Normal{Float64}(μ=0.0, σ=1.0))

julia> x_param = build_parameter(error, set, supports = [-0.5, 0.5])
InfOptParameter{DistributionSet{Normal{Float64}}}(DistributionSet{Normal{Float64}}(Normal{Float64}(μ=0.0, σ=1.0)), [-0.5, 0.5], false)
```
Again, we use [`add_parameter](@ref) to add `x_param` to the `InfiniteModel` and
assign it the name `x`:
```jldoctest rand_define
julia> x_ref = add_parameter(model, x_param, "x")
x
```

### Macro Definition
One user-friendly way of defining infinite parameters is by macro
[`@infinite_parameter`](@ref). Again, let's consider a time parameter
``t \in [0, 10]`` with supports `[0, 2, 5, 7, 10]`. Similar to
[`JuMP.@variable`](@ref), we can use comparison operators to set lower bounds
and upper bounds for the infinite parameter:
```jldoctest time_macro_define; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, 0 <= t <= 10, supports = [0, 2, 5, 7, 10])
t
```
More generally, we use `in` to define the set that an infinite parameter is
subject to. The set could be an interval set, or a distribution set. For example,
we can define the same parameter `t` as above in the following way:
```jldoctest time_macro_define_in; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10], supports = [0, 2, 5, 7, 10])
t
```
In a similar way, we can define a random infinite parameter subject to some
distribution. For example, a Gaussian infinite parameter with mean 0 and
standard deviation 1 can be defined by
```jldoctest rand_macro_define; setup = :(using InfiniteOpt, Distributions; model = InfiniteModel())
julia> dist = Normal(0., 1.)
Normal{Float64}(μ=0.0, σ=1.0)

julia> @infinite_parameter(model, x in dist, supports = [-0.5, 0.5])
x
```



things to add: containers, independent, num_supports (cases of error), base_name,
(non)anonymous (keywords, expression forms), multi-dim definition, mechanism
of macro (makes parameter, defines julia variable, register)

## Supports
For an infinite parameter, its supports are a finite set of points that the
parameter will take (or possibly take, if the parameter is random). During the
transcription stage, the supports specified will become part of the grid points
that approximate all functions parameterized by the infinite parameter.

Once an infinite parameter is defined, users can access the supports using
[`supports`](@ref) function:
```jldoctest time_macro_define
julia> supports(t_ref)
5-element Array{Int64,1}:
  0
  2
  5
  7
 10

```
We also provide functions that access other related information about the
supports. For example, [`has_supports`](@ref) checks whether a parameter has
supports, while [`num_supports`](@ref) gives the number of supports associated
with a parameter:
```jldoctest time_macro_define
julia> has_supports(t_ref)
true

julia> num_supports(t_ref)
5

```
Now suppose we want to add more supports to the `t`, which is already assigned
with some supports. We can use [`add_supports`](@ref) function to achieve this
goal:
```jldoctest time_macro_define
julia> add_supports(t_ref, [3, 8])

julia> supports(t_ref)
7-element Array{Int64,1}:
  0
  2
  5
  7
 10
  3
  8

```
At times we might want to change the supports completely. In those cases, the
function [`set_supports`](@ref) resets the supports for a certain parameter with
new supports provided:
```jldoctest time_macro_define
julia> set_supports(t_ref, [0,3,5,8,10], force = true)

julia> supports(t_ref)
5-element Array{Int64,1}:
  0
  3
  5
  8
 10

```
Note that the keyword argument [`force`] must be set as [`true`] if the
parameter has been assigned with supports.

### Automatic Support Generation During Parameter Definition
For the examples in the [Parameter Definition](@ref), we have seen how to
manually add supports to an infinite parameter. For a quick automatic
generation of support points, though, users do not have to input the support
points. Instead, the number of support points generated is supplied.

For an infinite parameter subject to an `IntervalSet`, uniformly spaced supports
including both ends are generated across the interval. For example, defining a
time parameter ``t \in [0, 10]`` with 4 supports using `build_parameter` gives
```jldoctest; setup = :(using InfiniteOpt)
julia> set = IntervalSet(0, 10)
IntervalSet(0.0, 10.0)

julia> t_param = build_parameter(error, set, num_supports = 4, sig_fig = 3)
InfOptParameter{IntervalSet}(IntervalSet(0.0, 10.0), [0.0, 3.33, 6.67, 10.0], false)
```  
Using macro definition we have
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, 0 <= t <= 10, num_supports = 4, sig_fig = 3)
t

julia> supports(t)
4-element Array{Float64,1}:
  0.0   
  3.33
  6.67
 10.0   

```
Note that the user can use the keyword argument `sig_fig` to dictate the
significant figures for the supports. The default value of `sig_fig` is 5.

For an infinite parameter that follows a univariate distribution,
supports are sampled from the underlying distribution. For example, we can
define an infinite parameter subject to a normal distribution with mean 0 and
variance 1:
```jldoctest; setup = :(using InfiniteOpt, Distributions; model = InfiniteModel(); dist = Normal(0., 1.))
julia> @infinite_parameter(model, x in dist, num_supports = 4)
x

julia> supports(x)
4-element Array{Float64,1}:
  0.51125
  0.94272
 -0.75805
  0.72906

```
For multivariate distributions, though, we require support points are provided
in the definition. However, we can use [`fill_in_supports!`](@ref) to generate
supports for parameters following multivariate distributions. See
[Automatic Support Generation For Defined Parameters](@ref) for details.

### Automatic Support Generation For Defined Parameters


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
group_id(::ParameterRef)
group_id(::AbstractArray{<:ParameterRef})
is_independent(::ParameterRef)
set_independent(::ParameterRef)
unset_independent(::ParameterRef)
parameter_by_name(::InfiniteModel, ::String)
all_parameters(::InfiniteModel)
```
