# Extensions
Here we provide guidance to various ways `InfinitOpt` can be extended.

## Overview
Extendibility is one of the core ideas of `InfinitOpt` so that it can serve as a
convenient tool for those developing and implementing advanced techniques for
infinite dimensional optimization problems. Thus, `InfiniteOpt` is developed in
a modular manner to readily accommodate user-defined functionality and/or to
serve as useful base in writing a `JuMP` extension. Admittedly, this modularity
might be imperfect and comments/suggestions are welcomed to help us improve this.

## Infinite Sets
Infinite sets are used to characterize the behavior of infinite parameters and
used to govern the behavior of supports in `InfiniteOpt`. Here we walk through
how user-defined sets can be added to various degrees of functionality. A
template is provided in [`./test/extensions/infinite_set.jl`](https://github.com/pulsipher/InfiniteOpt.jl/blob/master/test/extensions/infinite_set.jl).
The extension steps employed are:
1. Define the new `struct` infinite set type (only thing required as bare minimum)
2. Extend [`InfiniteOpt.supports_in_set`](@ref) (enables error checking of supports)
3. Extend [`InfiniteOpt.generate_support_values`](@ref) (enables support generation via `num_supports` keyword arguments)
4. If a lower bound and upper bound can be reported, extend `JuMP` lower bound and upper bound methods (enables automatic bound detection in `integral`)

As an example, let's create a univariate disjoint interval set as an infinite set type.
This corresponds to the set ``[lb_1, ub_1] \cup [lb_2, ub_2]`` where
``ub_1 \leq lb_2``. First, we need to create the `DataType` with inheritance from [`InfiniteScalarSet`](@ref):
```jldoctest set_ext; output = false
using InfiniteOpt, JuMP

struct DisjointSet <: InfiniteOpt.InfiniteScalarSet
    lb1::Float64
    ub1::Float64
    lb2::Float64
    ub2::Float64
    # constructor
    function DisjointSet(lb1::Number, ub1::Number, lb2::Number, ub2::Number)
        if lb1 > ub1 || lb2 > ub2 || ub1 > lb2
            error("Invalid bounds")
        end
        return new(convert(Float64, lb1), convert(Float64, ub1),
                   convert(Float64, lb2), convert(Float64, ub2))
    end
end

# output


```
Notice that we also define the constructor function to error check and convert as
needed (this is recommended, but not required). For basic functionality this is
all we have to do to add a set in `InfiniteOpt`.

We can now define infinite parameters using this set via
[`@infinite_parameter`](@ref) both anonymously and explicitly:
```jldoctest set_ext
julia> model = InfiniteModel();

julia> t = @infinite_parameter(model, set = DisjointSet(0, 1, 3, 4), base_name = "t")
t

julia> @infinite_parameter(model, t in DisjointSet(0, 1, 3, 4))
t
```
Once defined (without further extension), these parameters can be used as normal
with the following limitations:
- Supports must be specified manually (`num_supports` is not enabled)
- Supports will not be checked if they are in the domain of the infinite set
- Set bounds cannot be queried.
- The [`DiscreteMeasureData`](@ref) or [`FunctionalDiscreteMeasureData`](@ref) must be provided explicitly to evaluate measures
However, all of these limitations except for the last one can be eliminated by 
extending a few functions as outlined above. To address the last one, we need 
to extend [`generate_integral_data`](@ref). See [`Measure Evaluation Techniques`]
for details. 

To enable support domain checking which is useful to avoid strange bugs, we will
extend [`InfiniteOpt.supports_in_set`](@ref). This returns a `Bool` to indicate
if a vector of supports are in the set's domain:
```jldoctest set_ext; output = false
function InfiniteOpt.supports_in_set(supports::Union{Number, Vector{<:Number}},
                                     set::DisjointSet)::Bool
    return all((set.lb1 .<= supports .<= set.ub1) .| (set.lb2 .<= supports .<= set.ub2))
end

# output


```
Now the checks are enabled so, the following would yield an error because the
support is not in the set domain:
```jldoctest set_ext
julia> @infinite_parameter(model, set = DisjointSet(0, 1, 3, 4), supports = 2)
ERROR: At none:1: `@infinite_parameter(model, set = DisjointSet(0, 1, 3, 4), supports = 2)`: Supports violate the set domain bounds.
```

To enable automatic support generation via the `num_supports` keyword and with
functions such as [`fill_in_supports!`](@ref), we will extend
[`InfiniteOpt.generate_support_values`](@ref):
```jldoctest set_ext; output = false
struct DisjointGrid <: InfiniteOpt.PublicLabel end

function InfiniteOpt.generate_support_values(set::DisjointSet;
                                             num_supports::Int = InfiniteOpt.DefaultNumSupports,
                                             sig_digits::Int = InfiniteOpt.DefaultSigDigits)::Tuple{Vector{<:Real}, DataType}
    length_ratio = (set.ub1 - set.lb1) / (set.ub1 - set.lb1 + set.ub2 - set.lb2)
    num_supports1 = Int64(ceil(length_ratio * num_supports))
    num_supports2 = num_supports - num_supports1
    supports1 = collect(range(set.lb1, stop = set.ub1, length = num_supports1))
    supports2 = collect(range(set.lb2, stop = set.ub2, length = num_supports2))
    return round.([supports1; supports2], sigdigits = sig_digits), DisjointGrid
end

# output


```
Now automatic support generation is enabled, for example:
```jldoctest set_ext
julia> par = @infinite_parameter(model, set = DisjointSet(0, 2, 3, 4), num_supports = 10)
noname

julia> supports(par)
10-element Array{Float64,1}:
 0.0
 0.333333333333
 0.666666666667
 1.0
 1.33333333333
 1.66666666667
 2.0
 3.0
 3.5
 4.0
```

Finally, we can extend the appropriate `JuMP` upper and lower bound functions
if desired which are:
- [`JuMP.has_lower_bound`](@ref JuMP.has_lower_bound(::AbstractInfiniteSet))
- [`JuMP.lower_bound`](@ref JuMP.lower_bound(::AbstractInfiniteSet))
- [`JuMP.set_lower_bound`](@ref JuMP.set_lower_bound(::AbstractInfiniteSet, ::Union{Real, Vector{<:Real}}))
- [`JuMP.has_upper_bound`](@ref JuMP.has_upper_bound(::AbstractInfiniteSet))
- [`JuMP.upper_bound`](@ref JuMP.upper_bound(::AbstractInfiniteSet))
- [`JuMP.set_upper_bound`](@ref JuMP.set_upper_bound(::AbstractInfiniteSet, ::Union{Real, Vector{<:Real}}))
However, if we want `has_lower_bound = false` and `has_upper_bound = false` then
no extension is needed. For our current example we won't do this since lower
and upper bounds aren't exactly clear for a disjoint interval. Please refer to
the template in `./InfiniteOpt/test/extensions/infinite_set.jl` to see how this
is done.

## Measure Evaluation Techniques
Measure evaluation methods are used to dictate how to evaluate measures. Users
may wish to apply evaluation methods other than Monte Carlo sampling and/or
Gaussian quadrature methods. To create multiple measures using the same new
evaluation methods, users may want to embed the new evaluation method under the
[`integral`](@ref) function that
does not require explicit construction of [`AbstractMeasureData`](@ref).

The basic way to do that is to write a function that creates [`AbstractMeasureData`](@ref)
object, and pass the object to the [`measure`](@ref). For instance, let's
consider defining a function that
enables the definition of a uniform grid for a univariate or multivariate
infinite parameter in [`IntervalSet`](@ref). The function, denoted `uniform_grid`,
generates uniform grid points as supports for univariate parameter and each component of
independent multivariate parameter. The univariate version of this function
can be defined as follows:

```jldoctest measure_eval; output = false, setup = :(using JuMP, InfiniteOpt)
function uniform_grid(param::GeneralVariableRef, lb::Real, ub::Real, num_supports::Int)::DiscreteMeasureData
    increment = (ub - lb) / (num_supports - 1)
    supps = [lb + (i - 1) * increment for i in 1:num_supports]
    coeffs = ones(num_supports) / num_supports * (ub - lb)
    return DiscreteMeasureData(param, coeffs, supps, lower_bound = lb, upper_bound = ub)
end

# output
uniform_grid (generic function with 1 method)

```
It is necessary to pass infinite parameter reference since the
construction of measure data object needs parameter information. Now let's
apply the new `uniform_grid` function to infinite parameters in
intervals. We consider a time parameter `t` and 2D spatial parameter `x`, and
two variables `f(t)` and `g(x)` parameterized by `t` and `x`, respectively:
```jldoctest measure_eval
julia> m = InfiniteModel();

julia> @infinite_parameter(m, t in [0, 5]);

julia> @infinite_variable(m, f(t));
```
Now we can use `uniform_grid` to construct a [`DiscreteMeasureData`](@ref) and
create a measure using the measure data, as shown below:

```jldoctest measure_eval
julia> tdata = uniform_grid(t, 0, 5, 6)
DiscreteMeasureData{GeneralVariableRef,1,Float64}(t, [0.833333, 0.833333, 0.833333, 0.833333, 0.833333, 0.833333], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], UniqueMeasure{Val{Symbol("##858")}}, InfiniteOpt.default_weight, 0.0, 5.0, false)

julia> f_meas = measure(f, tdata)
measure{t ∈ [0, 5]}[f(t)]

julia> expand(f_meas)
0.8333333333333333 f(0) + 0.8333333333333333 f(1) + 0.8333333333333333 f(2) + 0.8333333333333333 f(3) + 0.8333333333333333 f(4) + 0.8333333333333333 f(5)
```

An alternate way of extending new measure evaluation methods is to extend
[`InfiniteOpt.MeasureToolbox.generate_integral_data`](@ref). This will
allow users to use their custom measure evaluation methods in the
[`integral`](@ref) function that does not explicitly require a measure data
object. A template for how such an extension is accomplished is provided in
[`./test/extensions/measure_eval.jl`](https://github.com/pulsipher/InfiniteOpt.jl/blob/master/test/extensions/measure_eval.jl).
In general, such an extension can be created as follows:
1. Define a new empty `struct` (e.g. `my_new_fn`) that dispatches your function
2. Extend [`InfiniteOpt.MeasureToolbox.generate_integral_data`](@ref),
where `method` is of the type `my_new_fn`, and `set` needs to be a subtype
of [`AbstractInfiniteSet`](@ref) that you wish to apply the new evaluation method
to.
Note that this procedure can be used to generate new measure evaluation methods not only for existing
infinite sets, but also for user-defined infinite sets. 

For example, an extension of [`InfiniteOpt.MeasureToolbox.generate_integral_data`](@ref)
that implements uniform grid for univariate and multivariate parameters in
[`IntervalSet`](@ref) can be created as follows:

```jldoctest measure_eval; output = false
const JuMPC = JuMP.Containers
struct UnifGrid <: InfiniteOpt.MeasureToolbox.AbstractUnivariateMethod end
function InfiniteOpt.MeasureToolbox.generate_integral_data(
    pref::InfiniteOpt.GeneralVariableRef,
    lower_bound::Real,
    upper_bound::Real,
    method::Type{UnifGrid};
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function = InfiniteOpt.default_weight
    )::InfiniteOpt.AbstractMeasureData # REPLACE WITH ACTUAL ALIAS
    increment = (upper_bound - lower_bound) / (num_supports - 1)
    supports = [lower_bound + (i - 1) * increment for i in 1:num_supports]
    coeffs = ones(num_supports) / num_supports * (upper_bound - lower_bound)
    return InfiniteOpt.DiscreteMeasureData(
        pref, coeffs, supports,
        weight_function = weight_func,
        lower_bound = lower_bound, 
        upper_bound = upper_bound)
end

# output


```

Also notice that users are free to pass keyword arguments for their new
functions in addition to the required positional arguments. This might be needed
in case if the new evaluation method requires additional information not
captured in the default positional arguments. For example, the multivariate
parameter version above needs to know if the multivariate parameter is
independent in order to throw a warning when needed.

We create measure for `f` and `g` using the `uniform_grid` method
```jldoctest measure_eval
julia> f_int = integral(f, t, num_supports = 6, eval_method = UnifGrid)
integral{t ∈ [0, 5]}[f(t)]

julia> expand(f_int)
0.8333333333333333 f(0) + 0.8333333333333333 f(1) + 0.8333333333333333 f(2) + 0.8333333333333333 f(3) + 0.8333333333333333 f(4) + 0.8333333333333333 f(5)
```
Here we go! We can freely use `UnifGrid` for infinite parameters in
[`IntervalSet`](@ref) now.

## Measure Data
Measures are used to evaluate over infinite domains. Users may wish to employ
measure abstractions that cannot be readily represented with coefficients and
discretized supports, and thus may wish to extend `InfiniteOpt`'s
measure framework to accommodate other paradigms. This can be accomplished my
implementing a user-defined measure data structure that inherits from
[`AbstractMeasureData`](@ref). A template for how such an extension is
accomplished is provided in [`./test/extensions/measure_data.jl`](https://github.com/pulsipher/InfiniteOpt.jl/blob/master/test/extensions/measure_data.jl).
The extension steps employed are:
1. Define the new data struct inheriting from [`AbstractMeasureData`](@ref) (required)
2. Extend [`InfiniteOpt.parameter_refs`](@ref parameter_refs(::AbstractMeasureData)) (required)
3. Extend [`InfiniteOpt.expand_measure`](@ref) (required)
4. Extend [`InfiniteOpt.supports`](@ref supports(::AbstractMeasureData)) (required if parameter supports are employed in any way)
5. Extend [`InfiniteOpt.add_supports_to_parameters`](@ref) (required if parameter supports are employed in measure evaluation)
5. Extend [`InfiniteOpt.measure_data_in_hold_bounds`](@ref) (enables hold variable bound checking with measures)
6. Extend [`InfiniteOpt.coefficients`](@ref) (useful getter method if applicable)
7. Extend [`InfiniteOpt.weight_function`](@ref) (useful getter method if applicable)
8. Extend [`InfiniteOpt.support_label`](@ref) (needed to enable deletion if supports are added.)
8. Make simple measure constructor wrapper of [`measure`](@ref) to ease definition.

To illustrate how this process can be done, let's consider extending `InfiniteOpt`
to include measure support for assessing the variance of random expressions. The
variance of an expression ``f(x, \xi)`` where ``x \in \mathbb{R}^n`` are hold
variables and ``\xi \in \mathbb{R}^m`` are random infinite parameters is defined:
```math
\mathbb{V}[f(x, \xi)] = \mathbb{E}\left[(f(x, \xi) - \mathbb{E}[f(x, \xi)])^2 \right].
```
Note, we could just accomplish this by nested use of [`expect`](@ref), but we
implement this example to illustrate the mechanics of extension.

First, let's define our new `struct` inheriting from `AbstractMeasureData`:
```jldoctest measure_data; output = false
using InfiniteOpt, JuMP, Distributions

const JuMPC = JuMP.Containers

struct DiscreteVarianceData <: AbstractMeasureData
    parameter_refs::Union{GeneralVariableRef, Vector{GeneralVariableRef}}
    supports::Vector
    label::DataType
    # constructor
    function DiscreteVarianceData(
        parameter_refs::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}},
        supports::Vector,
        label::DataType = InfiniteOpt.generate_unique_label())
        # convert input as necessary to proper array format
        if parameter_refs isa AbstractArray
            parameter_refs = convert(Vector, parameter_refs)
            supports = [convert(Vector, arr) for arr in supports]
        end
        return new(parameter_refs, supports, label)
    end
end

# output


```

We have defined our data type, so let's extend the measure data query
methods to enable its definition. These include
[`parameter_refs`](@ref parameter_refs(::AbstractMeasureData)), 
[`supports`](@ref supports(::AbstractMeasureData)), and 
[`support_label`](@ref support_label(::AbstractMeasureData)):
```jldoctest measure_data; output = false
function InfiniteOpt.parameter_refs(data::DiscreteVarianceData)
    return data.parameter_refs
end

function InfiniteOpt.supports(data::DiscreteVarianceData)::Vector
    return data.supports
end

function InfiniteOpt.support_label(data::DiscreteVarianceData)::DataType
    return data.label
end

# output


```

We also need to extend [`InfiniteOpt.add_supports_to_parameters`](@ref)
since support points will be used for measure evaluation later:
```jldoctest measure_data; output = false
function InfiniteOpt.add_supports_to_parameters(data::DiscreteVarianceData)::Nothing
    pref = parameter_refs(data)
    supps = supports(data)
    label = support_label(data)
    add_supports(pref, supps, label = label)
    return
end

# output


```
Note that extending `supports` is not needed for abstractions that don't involve
discretization of the infinite parameter(s), such as the case for certain
outer approximation techniques. 
Our extension is now sufficiently constructed to
allow us to define out the new variance measure via
[`measure`](@ref). For example,
```jldoctest measure_data; setup = :(using Random; Random.seed!(42))
# Setup the infinite model
model = InfiniteModel()
@infinite_parameter(model, xi in Normal(), num_supports = 2) # few for simplicity
@infinite_variable(model, y(xi))
@hold_variable(model, z)

# Define out new variance measure
data = DiscreteVarianceData(xi, supports(xi))
mref = measure(2y + z, data, name = "Var")

# output
Var{xi}[2 y(xi) + z]
```
Thus, we can define measure references that employ this our new data type.

We can define variance measures now, but now let's extend
[`expand_measure`](@ref) so that they can be expanded into finite expressions:
```jldoctest measure_data; output = false
function InfiniteOpt.expand_measure(expr::JuMP.AbstractJuMPScalar,
                                    data::DiscreteVarianceData,
                                    write_model::JuMP.AbstractModel
                                    )::JuMP.AbstractJuMPScalar
    # define the expectation data
    expect_data = DiscreteMeasureData(
                      data.parameter_refs,
                      1 / length(data.supports) * ones(length(data.supports)),
                      data.supports, is_expect = true, label = data.label)
    # define the mean
    mean = measure(expr, expect_data)
    # return the expansion of the variance using the data mean
    return expand_measure((copy(expr) - mean)^2, expect_data, write_model)
end

# output


```
Notice that we reformulated our abstraction in terms of measures with
[`DiscreteMeasureData`](@ref) so that we could leverage the existing
[`expand_measure`](@ref) library. Now, new the measure type can be expanded and
moreover infinite models using this new type can be optimized. Let's try
expanding the measure we already defined:
```jldoctest measure_data
julia> expand(mref)
2 y(-0.556026876146)² + 2 z*y(-0.556026876146) - 2 y(-0.556026876146)² - 2 y(-0.44438335711)*y(-0.556026876146) - 2 z*y(-0.556026876146) + 0 z² - y(-0.556026876146)*z - y(-0.44438335711)*z + 0.5 y(-0.556026876146)² + 0.5 y(-0.556026876146)*y(-0.44438335711) + 0.5 y(-0.556026876146)*z + 0.5 y(-0.44438335711)*y(-0.556026876146) + 0.5 y(-0.44438335711)² + 0.5 y(-0.44438335711)*z + 0.5 z*y(-0.556026876146) + 0.5 z*y(-0.44438335711) + 2 y(-0.44438335711)² + 2 z*y(-0.44438335711) - 2 y(-0.556026876146)*y(-0.44438335711) - 2 y(-0.44438335711)² - 2 z*y(-0.44438335711) - y(-0.556026876146)*z - y(-0.44438335711)*z + 0.5 y(-0.556026876146)² + 0.5 y(-0.556026876146)*y(-0.44438335711) + 0.5 y(-0.556026876146)*z + 0.5 y(-0.44438335711)*y(-0.556026876146) + 0.5 y(-0.44438335711)² + 0.5 y(-0.44438335711)*z + 0.5 z*y(-0.556026876146) + 0.5 z*y(-0.44438335711)
```

Finally, as per recommendation let's make a wrapper method to make defining
variance measures more convenient:
```jldoctest measure_data; output = false
function variance(expr::Union{JuMP.GenericAffExpr, GeneralVariableRef},
                  params::Union{GeneralVariableRef, AbstractArray{GeneralVariableRef}};
                  name::String = "Var", num_supports::Int = 10,
                  use_existing::Bool = false)::GeneralVariableRef
    # get the supports
    if use_existing
        supps = supports.(params)
    else
        supps = generate_support_values(infinite_set(first(params)),
                                           num_supports = num_supports)
    end
    # make the data
    data = DiscreteVarianceData(params, supps)
    # built the measure
    return measure(expr, data, name = name)
end

# output

variance (generic function with 1 method)
```
Notice in this case that we only permit linear expressions for `expr` since
it will be squared by our new measure and we currently only support quadratic
expressions. (This could be overcome by defining place a place holder variable
for `expr`.

Now let's use our constructor to repeat the above measure example:
```jldoctest measure_data
julia> expand(variance(2y + z, xi, use_existing = true))
2 y(-0.556026876146)² + 2 z*y(-0.556026876146) - 2 y(-0.556026876146)² - 2 y(-0.44438335711)*y(-0.556026876146) - 2 z*y(-0.556026876146) + 0 z² - y(-0.556026876146)*z - y(-0.44438335711)*z + 0.5 y(-0.556026876146)² + 0.5 y(-0.556026876146)*y(-0.44438335711) + 0.5 y(-0.556026876146)*z + 0.5 y(-0.44438335711)*y(-0.556026876146) + 0.5 y(-0.44438335711)² + 0.5 y(-0.44438335711)*z + 0.5 z*y(-0.556026876146) + 0.5 z*y(-0.44438335711) + 2 y(-0.44438335711)² + 2 z*y(-0.44438335711) - 2 y(-0.556026876146)*y(-0.44438335711) - 2 y(-0.44438335711)² - 2 z*y(-0.44438335711) - y(-0.556026876146)*z - y(-0.44438335711)*z + 0.5 y(-0.556026876146)² + 0.5 y(-0.556026876146)*y(-0.44438335711) + 0.5 y(-0.556026876146)*z + 0.5 y(-0.44438335711)*y(-0.556026876146) + 0.5 y(-0.44438335711)² + 0.5 y(-0.44438335711)*z + 0.5 z*y(-0.556026876146) + 0.5 z*y(-0.44438335711)
```

We have done it! Now go and extend away!

## [Optimizer Models] (@id extend_optimizer_model)
`InfiniteOpt` provides a convenient interface and abstraction for modeling
infinite-dimensional optimization problems. By default, `InfiniteModel`s are
reformulated into a solvable `JuMP.Model` (referred to as an optimizer model)
via `TranscriptionOpt` which discretizes the model in accordance with the
infinite parameter supports. However, users may wish to employ some other
reformulation method to produce the optimizer model. This section will explain
how this can be done in `InfiniteOpt`. A template for implementing this
extension is provided in [`./test/extensions/optimizer_model.jl`](https://github.com/pulsipher/InfiniteOpt.jl/blob/master/test/extensions/optimizer_model.jl).
Our default sub-module `InfiniteOpt.TranscriptionOpt` also serves as a good
example.

A new reformulation method and its corresponding optimizer model can be
extended using the following steps:
1. Define a `mutable struct` for variable/constraint mappings and other needed info (required)
2. Define a `JuMP.Model` constructor that uses (1.) in `Model.ext[:my_ext_key]` (recommended)
3. Extend [`build_optimizer_model!`](@ref) to in accordance with the new optimizer model (required)
4. Extend [`optimizer_model_variable`](@ref) if possible (enables result queries)
5. Extend [`optimizer_model_expression`](@ref) if possible (enables result queries)
6. Extend [`optimizer_model_constraint`](@ref) if possible (enables result queries)
7. Extend [`InfiniteOpt.variable_supports`](@ref) if appropriate
8. Extend [`InfiniteOpt.expression_supports`](@ref) if appropriate
9. Extend [`InfiniteOpt.constraint_supports`](@ref) if appropriate
10. If steps 4-6 are skipped then extend the following:
    - [`InfiniteOpt.map_value`](@ref) (enables `JuMP.value`)
    - [`InfiniteOpt.map_optimizer_index`](@ref) (enables `JuMP.optimizer_index`)
    - [`InfiniteOpt.map_dual`](@ref) (enables `JuMP.dual`)
    - [`InfiniteOpt.map_shadow_price`](@ref) (enables `JuMP.shadow_price`)
    - [`InfiniteOpt.map_lp_rhs_perturbation_range`](@ref) (enables `JuMP.lp_rhs_perturbation_range`)
    - [`InfiniteOpt.map_lp_objective_perturbation_range`](@ref) (enables `JuMP.lp_objective_perturbation_range`)
11. Extend [`InfiniteOpt.add_measure_variable`](@ref) to use [`expand_measure`](@ref) without modifying the infinite model
12. Extend [`InfiniteOpt.delete_reduced_variable`](@ref) to use [`expand_measure`](@ref) without modifying the infinite model and delete unneeded reduced variables.

For the sake of example, let's suppose we want to define a reformulation method
for `InfiniteModel`s that are 2-stage stochastic programs (i.e., only
`DistributionSet`s are used, infinite variables are random 2nd stage variables,
and hold variables are 1st stage variables). In particular, let's make a simple
method that replaces the infinite parameters with their mean values, giving us
the deterministic mean-valued problem.

First, let's define the `mutable struct` that will be used to store our variable
and constraint mappings. This case it is quite simple since our deterministic
model will have a 1-to-1 mapping:
```jldoctest opt_model; output = false
using InfiniteOpt, JuMP, Distributions

mutable struct DeterministicData
    # variable and constraint mapping
    infvar_to_detvar::Dict{GeneralVariableRef, VariableRef}
    infconstr_to_detconstr::Dict{InfOptConstraintRef, ConstraintRef}
    # constructor
    function DeterministicData()
        return new(Dict{GeneralVariableRef, VariableRef}(),
                   Dict{InfOptConstraintRef, ConstraintRef}())
    end
end

# output

```

Now let's define a constructor for optimizer models that will use
`DeterministicData` and let's define a method to access that data:
```jldoctest opt_model; output = false
const DetermKey = :DetermData

function DeterministicModel(args...; kwargs...)::Model
    # initialize the JuMP Model
    model = Model(args...; kwargs...)
    model.ext[DetermKey] = DeterministicData()
    return model
end

function deterministic_data(model::Model)::DeterministicData
    haskey(model.ext, DetermKey) || error("Model is not a DeterministicModel.")
    return model.ext[DetermKey]
end

# output
deterministic_data (generic function with 1 method)

```

!!! note
    The use of an extension key such as `DetermKey` is required since it used to
    dispatch reformulation and querying methods making optimizer model
    extensions possible.

With the constructor we can now specify that a given `InfiniteModel` uses a
`DeterministicModel` instead of a `TranscriptionModel` using the `OptimizerModel`
keyword argument or via [`set_optimizer_model`](@ref):
```jldoctest opt_model; output = false
using Ipopt

# Make model using Ipopt and DeterministicModels
model = InfiniteModel(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0),
                      OptimizerModel = DeterministicModel)

# Or equivalently
model = InfiniteModel()
set_optimizer_model(model, DeterministicModel())
set_optimizer(model, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))

# output

```
Now `model` uses a `DeterministicModel` as its optimizer model! With that we can
build our `InfiniteModel` as normal, for example:
```jldoctest opt_model
@infinite_parameter(model, ξ in Uniform())
@infinite_variable(model, y[1:2](ξ) >= 0)
@hold_variable(model, x)
@objective(model, Min, x + expect(y[1] + y[2], ξ))
@constraint(model, 2y[1] - x <= 42)
@constraint(model, y[2]^2 + ξ == 2)
print(model)

# output
Min x + expect{ξ}[y[1](ξ) + y[2](ξ)]
Subject to
 y[1](ξ) ≥ 0.0, ∀ ξ ~ Uniform
 y[2](ξ) ≥ 0.0, ∀ ξ ~ Uniform
 2 y[1](ξ) - x ≤ 42.0, ∀ ξ ~ Uniform
 y[2](ξ)² + ξ = 2.0, ∀ ξ ~ Uniform
```

We have defined our `InfiniteModel`, but now we need to specify how to
reformulate it into a `DeterministicModel`. This is accomplished by extending
[`build_optimizer_model!`](@ref). This will enable the use of `optimize!`. First,
let's define an internal function `_make_expression` that will use dispatch to
convert and `InfiniteOpt` expression into a `JuMP` expression using the mappings
stored in `opt_model` in its `DeterministicData`:
```jldoctest opt_model; output = false
## Make dispatch methods for converting InfiniteOpt expressions
# GeneralVariableRef
function _make_expression(opt_model::Model, expr::GeneralVariableRef)
    return _make_expression(opt_model, expr, index(expr))
end
# IndependentParameterRef
function _make_expression(opt_model::Model, expr::GeneralVariableRef, 
                          ::IndependentParameterIndex)
    return mean(infinite_set(expr).distribution) # assuming univariate
end
# FiniteParameterRef
function _make_expression(opt_model::Model, expr::GeneralVariableRef, 
                          ::FiniteParameterIndex)
    return parameter_value(expr)
end
# DependentParameterRef
function _make_expression(opt_model::Model, expr::GeneralVariableRef, 
                          ::DependentParameterIndex)
    return mean(infinite_set(expr).distribution) # assuming valid dist.
end
# DecisionVariableRef
function _make_expression(opt_model::Model, expr::GeneralVariableRef, 
                          ::Union{InfiniteVariableIndex, HoldVariableIndex})
    return deterministic_data(opt_model).infvar_to_detvar[expr]
end
# MeasureRef --> assume is expectation
function _make_expression(opt_model::Model, expr::GeneralVariableRef,
                          ::MeasureIndex)
    return _make_expression(opt_model, measure_function(expr))
end
# AffExpr
function _make_expression(opt_model::Model, expr::GenericAffExpr)
    return @expression(opt_model, sum(c * _make_expression(opt_model, v) 
                       for (c, v) in linear_terms(expr)) + constant(expr))
end
# QuadExpr
function _make_expression(opt_model::Model, expr::GenericQuadExpr)
    return @expression(opt_model, sum(c * _make_expression(opt_model, v1) * 
                       _make_expression(opt_model, v2) for (c, v1, v2) in quad_terms(expr)) + 
                       _make_expression(opt_model, expr.aff))
end

# output
_make_expression (generic function with 8 methods)

```
For simplicity in example, above we assume that only `DistributionSet`s are used,
there are not any `PointVariableRef`s, and all `MeasureRef`s correspond to expectations.
Naturally, a full extension should include checks to enforce that such assumptions
hold.

Now let's extend [`build_optimizer_model!`](@ref) for `DeterministicModel`s.
Such extensions should build an optimizer model in place and in general should
employ the following:
- [`clear_optimizer_model_build!`](@ref InfiniteOpt.clear_optimizer_model_build!(::InfiniteModel))
- [`set_optimizer_model_ready`](@ref).
In place builds without the use of `clear_optimizer_model_build!` are also
possible, but will require some sort of active mapping scheme to update in
accordance with the `InfiniteModel` in the case that the
optimizer model is built more than once. Thus, for simplicity we extend
`build_optimizer_model!` below using an initial clearing scheme:
```jldoctest opt_model; output = false
function InfiniteOpt.build_optimizer_model!(model::InfiniteModel,
                                            key::Val{DetermKey})
    # TODO check that `model` is a stochastic model
    # clear the model for a build/rebuild
    determ_model = InfiniteOpt.clear_optimizer_model_build!(model)

    # add variables
    for vref in all_variables(model)
        dvref = dispatch_variable_ref(vref)
        if dvref isa InfiniteVariableRef # have to handle the infinite variable functional start value
            inf_var = InfiniteOpt._core_variable_object(dvref)
            info = InfiniteOpt.TranscriptionOpt._format_infinite_info(inf_var, zeros(length(raw_parameter_refs(dvref))))
        else
            info = InfiniteOpt._variable_info(dvref)
        end
        new_vref = add_variable(determ_model, ScalarVariable(info),
                                name(dvref)) # TODO update infinite variable names
        deterministic_data(determ_model).infvar_to_detvar[vref] = new_vref
    end

    # add the objective
    set_objective(determ_model, objective_sense(model),
                   _make_expression(determ_model, objective_function(model)))

    # add the constraints
    for cref in all_constraints(model)
        if !InfiniteOpt._is_info_constraint(cref)
            constr = constraint_object(cref)
            new_constr = build_constraint(error, _make_expression(determ_model, constr.func),
                                          constr.set)
            new_cref = add_constraint(determ_model, new_constr, name(cref))
            deterministic_data(determ_model).infconstr_to_detconstr[cref] = new_cref
        end
    end

    # update the status
    set_optimizer_model_ready(model, true)
    return
end

# output

```
Now we can build our optimizer model to obtain a `DeterministicModel` which can
be leveraged to call `optimize!`
```jldoctest opt_model
optimize!(model)
print(optimizer_model(model))

# output
Min x + y[1] + y[2]
Subject to
 2 y[1] - x ≤ 42.0
 y[2]² = 1.5
 y[1] ≥ 0.0
 y[2] ≥ 0.0
```
Note that batter variable naming could be used with the reformulated infinite
variables. Moreover, in general extensions of [`build_optimizer_model!`](@ref)
should account for the possibility that `InfiniteModel` contains `HoldVariable`s
and/or `ScalarConstraint`s that contain [`ParameterBounds`](@ref) as accessed via
[`parameter_bounds`](@ref).

Now that we have optimized out `InfiniteModel` via the use the of a
`DeterministicModel`, we probably will want to access the results. All queries
are enabled when we extend [`optimizer_model_variable`](@ref), 
[`optimizer_model_expression`](@ref), and [`optimizer_model_constraint`](@ref) 
to return the variable(s)/expression(s)/constraint(s) in the
optimizer model corresponding to their `InfiniteModel` counterparts. These will
use the `mutable struct` of mapping data and should error if no mapping can be
found, Let's continue our example using `DeterministicModel`s:
```jldoctest opt_model; output = false
function InfiniteOpt.optimizer_model_variable(vref::GeneralVariableRef,
                                              key::Val{DetermKey};
                                              label = All)
    model = optimizer_model(JuMP.owner_model(vref))
    map_dict = deterministic_data(model).infvar_to_detvar
    haskey(map_dict, vref) || error("Variable $vref not used in the optimizer model.")
    return map_dict[vref]
end

function InfiniteOpt.optimizer_model_expression(expr::JuMP.AbstractJuMPScalar,
                                                key::Val{DetermKey};
                                                label = All)
    model = optimizer_model(JuMP.owner_model(vref))
    return _make_expression(model, expr)
end

function InfiniteOpt.optimizer_model_constraint(cref::InfOptConstraintRef,
                                                key::Val{DetermKey};
                                                label = All)
    model = optimizer_model(JuMP.owner_model(cref))
    map_dict = deterministic_data(model).infconstr_to_detconstr
    haskey(map_dict, cref) || error("Constraint $cref not used in the optimizer model.")
    return map_dict[cref]
end

# output

```
With these extensions we can now access all the result queries. For example,
```jldoctest opt_model
julia> termination_status(model)
LOCALLY_SOLVED::TerminationStatusCode = 4

julia> result_count(model)
1

julia> value.(y)
2-element Array{Float64,1}:
 0.0
 1.224744871391589

julia> optimizer_index(x)
MathOptInterface.VariableIndex(3)
```

!!! note
    If [`optimizer_model_variable`](@ref), [`optimizer_model_expression`](@ref), 
    and/or [`optimizer_model_constraint`](@ref)
    cannot be extended due to the nature of the reformulation then please refer
    to step 10 of the extension steps listed at the beginning of this section.

Furthermore, if appropriate for the given reformulation the following should be
extended:
- [`InfiniteOpt.variable_supports`](@ref) to enable `supports` on variables)
- [`InfiniteOpt.expression_supports`](@ref) to enable `supports` on expressions)
- [`InfiniteOpt.constraint_supports`](@ref) to enable `supports` on constraints)

That's it!

## Wrapper Packages
`InfiniteOpt` provides a convenient modular interface for defining infinite
dimensional optimization problems, implementing many tedious `JuMP` extensions
such as facilitating mixed variable expressions. Thus, `InfiniteOpt` can serve
as a base package for specific types of infinite dimensional problems and/or
solution techniques. These extension packages can implement any of the extensions
shown above and likely will want to introduce wrapper functions and macros to
use package specific terminology (e.g., using random variables instead of
infinite variables).

TODO refer to example packages that do this (e.g., an update of FlexibilityAnalysis.jl)
