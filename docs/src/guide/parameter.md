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
Note that `add_parameter` does not register the name of the parameters into the
model that it adds to. As shown in [Macro Definition](@ref), the macro definition
does not allow for multiple parameters sharing the same name and will throw an
error if it happens.

### Macro Definition
#### One-Dimensional Parameters
One user-friendly way of defining infinite parameters is by macro
[`@infinite_parameter`](@ref). The macro executes the same process as the
manual definition (steps listed in [Parameter Definition](@ref)), but allows
the users to manipulate several features of the defined infinite parameters.
Again, let's consider a time parameter ``t \in [0, 10]`` with supports
`[0, 2, 5, 7, 10]`. Similar to [`JuMP.@variable`](@ref), we can use comparison
operators to set lower bounds and upper bounds for the infinite parameter:
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
Additional ways of defining infinite parameters are provided using keyword
arguments. For example, we can use `lower_bound` and `upper_bound` to define an
infinite parameter in an interval set:
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t, lower_bound = 0, upper_bound = 10, supports = [0, 2, 5, 7, 10])
t
```
A bit more generally, we can also use `set` to directly input the
`AbstractInfiniteSet` that the parameter is in. For example:
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t, set = IntervalSet(0, 10), supports = [0, 2, 5, 7, 10])
t
```
The parameter definition methods using keyword arguments will be useful later
when we introduce how to define anonymous parameters. See
[Anonymous Parameter Definition](@ref) for more details.

All the definitions above return a [`ParameterRef`](@ref) that refer to the
defined parameter. Note that we can also ignore the `supports` keyword argument
and the macro will define an empty array of supports for that parameter.

#### Multi-Dimensional Parameter
Using macro definition, we can also define multi-dimensional infinite parameters
in a concise way. For example, consider a position parameter `x` in a
3-dimensional space constrained in a unit cube (i.e. in the interval `[0, 1]`  
for all dimensions). This parameter can be defined in one line as follows:
```jldoctest 3d_macro; setup = :(using InfiniteOpt; model= InfiniteModel())
julia> @infinite_parameter(model, x[1:3] in [0, 1], supports = [0.3, 0.7])
3-element Array{ParameterRef,1}:
 x[1]
 x[2]
 x[3]
```
Note that we can set different supports to different dimension using an index
inside the macro similar to [`JuMP.variable`](@ref). For example,
```jldoctest; setup = :(using InfiniteOpt; model= InfiniteModel())
julia> points = [0.2 0.8; 0.3 0.7]
2×2 Array{Float64,2}:
 0.2  0.8
 0.3  0.7

julia> @infinite_parameter(model, a[i = 1:2] in [0, 1], supports = points[i,:])
2-element Array{ParameterRef,1}:
 a[1]
 a[2]

julia> supports(a[1])
2-element Array{Float64,1}:
 0.2
 0.8
```
In a similar way we can define an infinite parameter subject to a multivariate
distribution. For example, a 2-dimensional parameter `xi` subject to a
2-D normal distribution can be created as follows:
```jldoctest; setup = :(using InfiniteOpt, Distributions; model = InfiniteModel())
julia> dist = MvNormal([0., 0.], [1. 0.; 0. 2.])
FullNormal(
dim: 2
μ: [0.0, 0.0]
Σ: [1.0 0.0; 0.0 2.0]
)

julia> @infinite_parameter(model, xi[1:2] in dist)
2-element Array{ParameterRef,1}:
 xi[1]
 xi[2]
```

#### Containers for Multi-Dimensional Parameters
Note that for all the cases of multi-dimensional parameter definition above, the
macro always returns an `Array` of [`ParameterRef`](@ref). For most cases this is
true. However, we can explicitly dictate the kind of containers we want to hold
the defined parameters using the keyword `container`. For example, we use
[`SparseAxisArray`](@ref) from the [`JuMP`](@ref) package for the space
parameter `x`:
```jldoctest; setup = :(using InfiniteOpt, JuMP; model= InfiniteModel())
julia> @infinite_parameter(model, x[1:3] in [0, 1], container = SparseAxisArray)
JuMP.Containers.SparseAxisArray{ParameterRef,1,Tuple{Any}} with 3 entries:
  [3]  =  x[3]
  [2]  =  x[2]
  [1]  =  x[1]
```

#### `independent` for Multi-Dimensional Parameters
For examples up to now we did not specify the value for the keyword
`independent`, which is set as `false` by default. The keyword `independent`
applies to multi-dimensional infinite parameters and dictates whether the
supports for different dimensions are independent. Setting `independent` as
`true` would be useful if the users want to generate a grid of supports for
a multi-dimensional parameter. For example, consider the position parameter `x`
in a 3D space. Say `x` is bounded in `[0, 1]` in all three dimensions, and the
user wants to generate grid points with interval `0.5` in all three dimensions.
In this case, we can define `x` in the following way:
```jldoctest; setup = :(using InfiniteOpt, JuMP; model= InfiniteModel())
julia> pts = collect(range(0, stop = 1, length = 3))
3-element Array{Float64,1}:
 0.0
 0.5
 1.0

julia> @infinite_parameter(model, x[1:3] in [0, 1], supports = pts, independent = true)
3-element Array{ParameterRef,1}:
 x[1]
 x[2]
 x[3]
```
If `independent` is set as `false`, the transcription step will generate
JuMP variables for values of any variable parameterized by `x` at `[0.0, 0.0, 0.0]`,
`[0.5, 0.5, 0.5]` and `[1.0, 1.0, 1.0]`, a total of 3 transcribed variables.
Instead, if `independent` is set as `true`, the transcription step will obtain
a unique permutation of these supports and each transcribe parameterized
variable accordingly, leading to a total of 27 transcribed variables in this
case.

#### Anonymous Parameter Definition and `base_name`
As mentioned above, we can define anonymous parameters using keyword arguments
in the macro [`@infinite_parameter`](@ref). For instance, we can create an
anonymous position parameter in a 3D space, referred to by a list of [`ParameterRef`](@ref)
called `x`:
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> x = @infinite_parameter(m, [1:3], lower_bound = 0, upper_bound = 1)
3-element Array{ParameterRef,1}:
 noname
 noname
 noname

julia> typeof(x)
Array{ParameterRef,1}

julia> name(x[1])
""
```
This syntax creates a 1D parameter if the part `[1:3]` is neglected.

Note that this macro definition automatically assigns an empty string to the
`base_name`. We can also assign a nontrivial base name to an anonymous parameter
using the keyword argument `base_name`. For example,
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(m, [1:3], lower_bound = 0, upper_bound = 1, base_name = "x")
3-element Array{ParameterRef,1}:
 x[1]
 x[2]
 x[3]

julia> @infinite_parameter(m, [1:3], lower_bound = -1, upper_bound = 0, base_name = "x")
3-element Array{ParameterRef,1}:
 x[1]
 x[2]
 x[3]
```
We can see that anonymous parameter definition allows for multiple parameters
sharing the same base name. This is not permitted with non-anonymous parameter
definition. In fact, in anonymous parameter definition, the macro does not
register the name of the parameters in the model, so when the model checks for
repeated names it will not detect the `x`. Refer to
[Detailed Mechanism of Macro Definition](@ref) if more details are desired.

#### Detailed Mechanism of Macro Definition

things to add:
num_supports (cases of error)
mechanism of macro (makes parameter, defines julia variable, register)

## Supports
For an infinite parameter, its supports are a finite set of points that the
parameter will take (or possibly take, if the parameter is random). During the
transcription stage, the supports specified will become part of the grid points
that approximate all functions parameterized by the infinite parameter.

Once an infinite parameter is defined, users can access the supports using
[`supports`](@ref) function:
```jldoctest supports; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, 0 <= t <= 10, supports = [0, 2, 5, 7, 10])
t

julia> supports(t)
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
```jldoctest supports
julia> has_supports(t)
true

julia> num_supports(t)
5
```
Now suppose we want to add more supports to the `t`, which is already assigned
with some supports. We can use [`add_supports`](@ref) function to achieve this
goal:
```jldoctest supports
julia> add_supports(t, [3, 8])

julia> supports(t)
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
```jldoctest supports
julia> set_supports(t, [0,3,5,8,10], force = true)

julia> supports(t)
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
