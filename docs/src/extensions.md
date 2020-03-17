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
template is provided in `./InfiniteOpt/test/extensions/infinite_set.jl`. The
extension steps employed are:
1. Define the new `struct` infinite set type (only thing required as bare minimum)
2. Extend [`InfiniteOpt.supports_in_set`](@ref) (enables error checking of supports)
3. Extend [`InfiniteOpt.generate_support_values`](@ref) (enables support generation via `num_supports` keyword arguments)
4. Extend [`InfiniteOpt.MeasureEvalMethods.generate_supports_and_coeffs`](@ref) and register the new set type via [`register_eval_method`](@ref) (enables measure evaluation methods)
5. If a lower bound and upper bound can be reported, extend `JuMP` lower bound and upper bound methods

As an example, let's create a disjoint interval set as an infinite set type.
This corresponds to the set ``[lb_1, ub_1] \cup [lb_2, ub_2]`` where
``ub_1 \leq lb_2``. First, we need to create the `DataType` with inheritance from [`AbstractInfiniteSet`](@ref):
```jldoctest set_ext; output = false
using InfiniteOpt, JuMP

struct DisjointSet <: InfiniteOpt.AbstractInfiniteSet
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
- The [`DiscreteMeasureData`](@ref) must be provided explicitly to evaluate measures
- Set bounds cannot be queried.
However, all of these limitations can be eliminated by extending a few functions
as outlined above.

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
ERROR: In `@infinite_parameter(model, set = DisjointSet(0, 1, 3, 4), supports = 2)`: Supports violate the set domain bounds.
```

To enable automatic support generation via the `num_supports` keyword and with
functions such as [`fill_in_supports!`](@ref), we will extend
[`InfiniteOpt.generate_support_values`](@ref):
```jldoctest set_ext; output = false
function InfiniteOpt.generate_support_values(set::DisjointSet;
                                             num_supports::Int = 50,
                                             sig_fig::Int = 5)::Array
    length_ratio = (set.ub1 - set.lb1) / (set.ub1 - set.lb1 + set.ub2 - set.lb2)
    num_supports1 = Int64(ceil(length_ratio * num_supports))
    num_supports2 = num_supports - num_supports1
    supports1 = collect(range(set.lb1, stop = set.ub1, length = num_supports1))
    supports2 = collect(range(set.lb2, stop = set.ub2, length = num_supports2))
    return round.([supports1; supports2], sigdigits = sig_fig)
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
 0.33333
 0.66667
 1.0
 1.3333
 1.6667
 2.0
 3.0
 3.5
 4.0
```
Note that in some cases [`generate_support_values`](@ref) might not be extendable
since returning supports for a single infinite parameter may not be appropriate.
This is typically true when a single set is used for multi-dimensional infinite
parameters. In these cases, [`generate_and_add_supports!`](@ref) will need to
be extended. This is exemplified with multivariate distribution sets that
correspond to an array of parameters where `generate_and_add_supports!` is
extended to generate the supports and then add the appropriate slices to each
dependent infinite parameter.

Next we can enable typical measure definition by extending
[`InfiniteOpt.MeasureEvalMethods.generate_supports_and_coeffs`](@ref) for our
new set type. This function creates support values and corresponding coefficients
for a measure using the `method` argument if appropriate. Also, if the `method`
argument is used then we'll need to register valid method calls for our new set
via [`register_eval_method`](@ref). If we wish to ignore the `method` then we'll
need to set `check_method = false` when calling [`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::Union{ParameterRef, AbstractArray{<:ParameterRef}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing}, ::Union{Number, AbstractArray{<:Number}, Nothing}))
or using [`set_measure_defaults`](@ref). Continuing our example, we obtain:
```jldoctest set_ext; output = false
const JuMPC = JuMP.Containers

function InfiniteOpt.MeasureEvalMethods.generate_supports_and_coeffs(
    set::DisjointSet,
    params::Union{InfiniteOpt.ParameterRef, AbstractArray{<:InfiniteOpt.ParameterRef}},
    num_supports::Int,
    lb::Union{Number, JuMPC.SparseAxisArray, Nothing},
    ub::Union{Number, JuMPC.SparseAxisArray, Nothing},
    method::Function)::Tuple
    length_ratio = (set.ub1 - set.lb1) / (set.ub1 - set.lb1 + set.ub2 - set.lb2)
    num_supports1 = Int64(ceil(length_ratio * num_supports))
    num_supports2 = num_supports - num_supports1
    supports1, coeffs1 = method(set.lb1, set.ub1, num_supports1)
    supports2, coeffs2 = method(set.lb2, set.ub2, num_supports2)
    return ([supports1; supports2], [coeffs1; coeffs2])
end

register_eval_method(model, DisjointSet, [mc_sampling, gauss_legendre])

# output


```
Now we can define measures as normal:
```jldoctest set_ext
julia> mref = measure(t^2, t)
measure(t²)
```

Finally, we can extend the appropriate `JuMP` upper and lower bound functions
if desired which are:
- [`JuMP.has_lower_bound`](@ref JuMP.has_lower_bound(::AbstractInfiniteSet))
- [`JuMP.lower_bound`](@ref JuMP.lower_bound(::AbstractInfiniteSet))
- [`JuMP.set_lower_bound`](@ref JuMP.set_lower_bound(::AbstractInfiniteSet, ::Number))
- [`JuMP.has_upper_bound`](@ref JuMP.has_upper_bound(::AbstractInfiniteSet))
- [`JuMP.upper_bound`](@ref JuMP.upper_bound(::AbstractInfiniteSet))
- [`JuMP.set_upper_bound`](@ref JuMP.set_upper_bound(::AbstractInfiniteSet, ::Number))
However, if we want `has_lower_bound = false` and `has_upper_bound = false` then
no extension is needed. For our current example we won't do this since lower
and upper bounds aren't exactly clear for a disjoint interval. Please refer to
the template in `./InfiniteOpt/test/extensions/infinite_set.jl` to see how this
is done.

## Measure Evaluation Techniques

## Measure Data
Measures are used to evaluate over infinite domains. Users may wish to employ
measure abstractions that cannot be readily represented with coefficients and
discretized supports, and thus may wish to extend `InfiniteOpt`'s
measure framework to accommodate other paradigms. This can be accomplished my
implementing a user-defined measure data structure that inherits from
[`AbstractMeasureData`](@ref). A template for how such an extension is
accomplished is provided in `./InfiniteOpt/test/extensions/measure_data.jl`. The
extension steps employed are:
1. Define the new data struct inheriting from [`AbstractMeasureData`](@ref) (required)
2. Extend [`InfiniteOpt.parameter_refs`](@ref parameter_refs(::AbstractMeasureData)) (required)
3. Extend [`InfiniteOpt.expand_measure`](@ref) (required)
4. Extend [`InfiniteOpt.supports`](@ref supports(::AbstractMeasureData)) (required if parameter supports are employed in any way)
5. Extend [`InfiniteOpt.measure_data_in_hold_bounds`](@ref) (enables hold variable bound checking with measures)
6. Extend [`InfiniteOpt.measure_name`](@ref) (enables meaningful measure naming)
7. Make simple measure constructor wrapper of [`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::AbstractMeasureData)) to ease definition.

To illustrate how this process can be done, let's consider extending `InfiniteOpt`
to include measure for assessing the variance of random expressions. The
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
    parameter_refs::Union{ParameterRef, JuMPC.SparseAxisArray{<:ParameterRef}}
    supports::Vector
    name::String
    # constructor
    function DiscreteVarianceData(
        parameter_refs::Union{ParameterRef, AbstractArray{<:ParameterRef}},
        supports::Vector, name::String = "Var")
        # convert input as necessary to proper array format
        if parameter_refs isa AbstractArray
            parameter_refs = convert(JuMPC.SparseAxisArray, parameter_refs)
            supports = [convert(JuMPC.SparseAxisArray, arr) for arr in supports]
        end
        return new(parameter_refs, supports, name)
    end
end

# output


```

We have defined our data type, so let's extend the measure data query
methods to enable its definition. These include
[`parameter_refs`](@ref parameter_refs(::AbstractMeasureData)),
[`supports`](@ref supports(::AbstractMeasureData)), and
[`measure_name`](@ref) (optional):
```jldoctest measure_data; output = false
function InfiniteOpt.parameter_refs(data::DiscreteVarianceData)
    return data.parameter_refs
end

function InfiniteOpt.supports(data::DiscreteVarianceData)::Vector
    return data.supports
end

function InfiniteOpt.measure_name(data::DiscreteVarianceData)::String
    return data.name
end

# output


```
Note that extending `supports` is not needed for abstractions that don't involve
discretization of the infinite parameter(s), such as the case for outer certain
outer approximation techniques. Our extension is now sufficiently constructor to
allow us to define out new variance measure via
[`measure`](@ref measure(::JuMP.AbstractJuMPScalar, ::AbstractMeasureData)). For
example,
```jldoctest measure_data; setup = :(using Random; Random.seed!(42))
# Setup the infinite model
model = InfiniteModel()
@infinite_parameter(model, xi in Normal(), num_supports = 2) # few for simplicity
@infinite_variable(model, y(xi))
@hold_variable(model, z)

# Define out new variance measure
data = DiscreteVarianceData(xi, supports(xi))
mref = measure(2y + z, data)

# output
Var(2 y(xi) + z)
```
Thus, we can define measure references that employ this our new data type.

We can define variance measures now, but now let's extend
[`expand_measure`](@ref) so that they can be expanded into finite expressions:
```jldoctest measure_data; output = false
function InfiniteOpt.expand_measure(expr::JuMP.AbstractJuMPScalar,
                                    data::DiscreteVarianceData,
                                    write_model::JuMP.AbstractModel,
                                    point_mapper::Function)::JuMP.AbstractJuMPScalar
    # define the expectation data
    expect_data = DiscreteMeasureData(
                      data.parameter_refs,
                      1 / length(data.supports) * ones(length(data.supports)),
                      data.supports, name = data.name)
    # define the mean
    mean = measure(expr, expect_data)
    # return the expansion of the variance using the data mean
    return expand_measure((copy(expr) - mean)^2, expect_data, write_model, point_mapper)
end

# output


```
Notice that we reformulated our abstraction in terms of measures with
[`DiscreteMeasureData`](@ref) so that we could leverage the existing
[`expand_measure`](@ref) library. Now, new measure type can be expanded and
moreover infinite models using this new type can be optimized. Let's try
expanding the measure we already defined:
```jldoctest measure_data
julia> expand(mref)
2 y(-0.55603)² + 2 y(-0.44438)² + 2 z*y(-0.55603) + 2 z*y(-0.44438) - 4 y(-0.55603)² - 4 y(-0.44438)*y(-0.55603) - 4 z*y(-0.55603) + 0 z² - 2 y(-0.55603)*z - 2 y(-0.44438)*z + y(-0.55603)² + y(-0.44438)²
```

Finally, as per recommendation let's make a wrapper method to make defining
variance measures more convenient:
```jldoctest measure_data; output = false
function variance(expr::Union{JuMP.GenericAffExpr, GeneralVariableRef},
                  params::Union{ParameterRef, AbstractArray{<:ParameterRef}};
                  name::String = "Var", num_supports::Int = 10,
                  use_existing::Bool = false)::MeasureRef
    # get the supports
    if use_existing
        supps = supports.(params)
    else
        supps = generate_support_values(infinite_set(first(params)),
                                           num_supports = num_supports)
    end
    # make the data
    data = DiscreteVarianceData(params, supps, name)
    # built the measure
    return measure(expr, data)
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
2 y(-0.55603)² + 2 y(-0.44438)² + 2 z*y(-0.55603) + 2 z*y(-0.44438) - 4 y(-0.55603)² - 4 y(-0.44438)*y(-0.55603) - 4 z*y(-0.55603) + 0 z² - 2 y(-0.55603)*z - 2 y(-0.44438)*z + y(-0.55603)² + y(-0.44438)²
```

We have done it! Now go and extend away!

## [Optimizer Models] (@id extend_optimizer_model)
`InfiniteOpt` provides a convenient interface and abstraction for modeling
infinite dimensional optimization problems. By default, `InfiniteModel`s are
reformulated into a solvable `JuMP.Model` (referred to as an optimizer model)
via `TranscriptionOpt` which discretizes the model in accordance with the
infinite parameter supports. However, users may wish to employ some other
reformulation method to produce the optimizer model. This section will explain
how this can be done in `InfiniteOpt`. A template for implementing this
extension is provided in `./InfiniteOpt/test/extensions/optimizer_model.jl`.
Our default sub-package `InfiniteOpt.TranscriptionOpt` also serves as a good
example.

A new reformulation method and its corresponding optimizer model can be
extended using the following steps:
1. Define a `mutable struct` for variable/constraint mappings and other needed info (required)
2. Define a `JuMP.Model` constructor that uses (1.) in `Model.ext[:my_ext_key]` (recommended)
3. Extend [`build_optimizer_model!`](@ref) to in accordance with the new optimizer model (required)
4. Extend [`optimizer_model_variable`](@ref) if possible (enables result queries)
5. Extend [`optimizer_model_constraint`](@ref) if possible (enables result queries)
6. Extend [`variable_supports`](@ref) if appropriate
7. Extend [`constraint_supports`](@ref) and [`constraint_parameter_refs`](@ref) if appropriate
8. If steps 4 & 5 are skipped then extend the following:
    - [`map_value`](@ref) (enables `JuMP.value`)
    - [`map_optimizer_index`](@ref) (enables `JuMP.optimizer_index`)
    - [`map_dual`](@ref) (enables `JuMP.dual`)
    - [`map_shadow_price`](@ref) (enables `JuMP.shadow_price`)
    - [`map_lp_rhs_perturbation_range`](@ref) (enables `JuMP.lp_rhs_perturbation_range`)
    - [`map_lp_objective_perturbation_range`](@ref) (enables `JuMP.lp_objective_perturbation_range`)

For the sake of example, let's suppose we want to define a reformulation method
for `InfiniteModel`s that are 2-stage stochastic programs (i.e., only
`DistributionSet`s are used, infinite variables are random 2nd stage variables,
and hold variables are 1st stage variables). In particular, let's make a simple
method that replaces the infinite parameters with their mean values, giving us
the deterministic mean-valued problem.

First, let's define the `mutable struct` that will be used to store our variable
and constraint mappings. This this case it is quite simple since our deterministic
model will have a 1-to-1 mapping:
```jldoctest opt_model; output = false
using InfiniteOpt, JuMP, Distributions

mutable struct DeterministicData
    # variable and constraint mapping
    infvar_to_detvar::Dict{GeneralVariableRef, VariableRef}
    infconstr_to_detconstr::Dict{GeneralConstraintRef, ConstraintRef}
    # constructor
    function DeterministicData()
        return new(Dict{GeneralVariableRef, VariableRef}(),
                   Dict{GeneralConstraintRef, ConstraintRef}())
    end
end

# output

```

Now let's define a constructor for optimizer models that will use
`DeterministicData` and let's define a method to access that data:
```jldoctest opt_model; output = false
function DeterministicModel(args...; kwargs...)::Model
    # initialize the JuMP Model
    model = Model(args...; kwargs...)
    model.ext[:DetermData] = DeterministicData()
    return model
end

function deterministic_data(model::Model)::DeterministicData
    haskey(model.ext, :DetermData) || error("Model is not a DeterministicModel.")
    return model.ext[:DetermData]
end

# output
deterministic_data (generic function with 1 method)

```

!!! note
    The use of an extension key such as `:DetermData` is key since it used to
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
Min x + expect(y[1](ξ) + y[2](ξ))
Subject to
 y[1](ξ) ≥ 0.0
 y[2](ξ) ≥ 0.0
 2 y[1](ξ) - x ≤ 42.0
 y[2](ξ)² + ξ = 2.0
 ξ ∈ Uniform(a=0.0, b=1.0)
```

We have defined our `InfiniteModel`, but now we need to specify how to
reformulate it into a `DeterministicModel`. This is accomplished by extending
[`build_optimizer_model!`](@ref). This will enable the use of `optimize!`. First,
let's define an internal function `_make_expression` that will use dispatch to
convert and `InfiniteOpt` expression into a `JuMP` expression using the mappings
stored in `opt_model` in its `DeterministicData`:
```jldoctest opt_model; output = false
## Make dispatch methods for converting InfiniteOpt expressions
# ParameterRef
function _make_expression(opt_model::Model, expr::ParameterRef)
    return mean(infinite_set(expr).distribution) # assuming univariate
end
# InfOptVariableRef
function _make_expression(opt_model::Model, expr::InfOptVariableRef)
    return deterministic_data(opt_model).infvar_to_detvar[expr]
end
# MeasureRef --> assume is expectation
function _make_expression(opt_model::Model, expr::MeasureRef)
    return _make_expression(opt_model, measure_function(expr))
end
# AffExpr
function _make_expression(opt_model::Model, expr::GenericAffExpr)
    new_expr = zero(AffExpr)
    for (vref, coef) in expr.terms
        add_to_expression!(new_expr, coef, _make_expression(opt_model, vref))
    end
    new_expr.constant += expr.constant
    return new_expr
end
# QuadExpr
function _make_expression(opt_model::Model, expr::GenericQuadExpr)
    new_expr = zero(QuadExpr)
    for (pair, coef) in expr.terms
        add_to_expression!(new_expr, coef, _make_expression(opt_model, pair.a),
                           _make_expression(opt_model, pair.b))
    end
    new_expr.aff = _make_expression(opt_model, expr.aff)
    return new_expr
end

# output
_make_expression (generic function with 5 methods)

```
For simplicity in example, above we assume that only `DistributionSet`s are used,
there are `PointVariableRef`s, and all `MeasureRef`s correspond to expectations.
Naturally, a full extension should include checks to enforce that such assumptions
hold.

Now let's extend [`build_optimizer_model!`](@ref) for `DeterministicModel`s.
Such extensions should build an optimizer model in place and in general should
employ the following:
- [`clear_optimizer_model_build!`](@ref )
- [`set_optimizer_model_ready`](@ref).
In place builds without the use of `clear_optimizer_model_build!` are also
possible, but will require some sort of active mapping scheme to update in
accordance with the `InfiniteModel` in the case that the
optimizer model is built more than once. Thus, for simplicity we extend
`build_optimizer_model!` below using an initial clearing scheme:
```jldoctest opt_model; output = false
function InfiniteOpt.build_optimizer_model!(model::InfiniteModel,
                                            key::Val{:DetermData})
    # TODO check that `model` is a stochastic model
    # clear the model for a build/rebuild
    determ_model = clear_optimizer_model_build!(model)

    # add variables
    for vref in all_variables(model)
        new_vref = add_variable(determ_model,
                                ScalarVariable(model.vars[index(vref)].info),
                                name(vref)) # TODO update infinite variable names
        deterministic_data(determ_model).infvar_to_detvar[vref] = new_vref
    end

    # add the objective
    set_objective(determ_model, objective_sense(model),
                   _make_expression(determ_model, objective_function(model)))

    # add the constraints
    for cref in all_constraints(model)
        if !model.constr_in_var_info[index(cref)]
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
Min x + y[1](ξ) + y[2](ξ)
Subject to
 2 y[1](ξ) - x ≤ 42.0
 y[2](ξ)² = 1.5
 y[1](ξ) ≥ 0.0
 y[2](ξ) ≥ 0.0
```
Note that batter variable naming could be used with the reformulated infinite
variables. Moreover, in general extensions of [`build_optimizer_model!`](@ref)
should account for the possibility that `InfiniteModel` contains `HoldVariable`s
and/or `ScalarConstraint`s that contain [`ParameterBounds`](@ref) as accessed via
[`parameter_bounds`](@ref).

Now that we have optimized out `InfiniteModel` via the use the of a
`DeterministicModel`, we probably will want to access the results. All queries
are enabled when we extend [`optimizer_model_variable`](@ref) and
[`optimizer_model_constraint`](@ref) return the variable(s)/constraint(s) in the
optimizer model corresponding to their `InfiniteModel` counterparts. These will
use the `mutable struct` of mapping data and should error if no mapping can be
found, Let's continue our example using `DeterministicModel`s:
```jldoctest opt_model; output = false
function InfiniteOpt.optimizer_model_variable(vref::InfOptVariableRef,
                                              key::Val{:DetermData})
    model = optimizer_model(JuMP.owner_model(vref))
    map_dict = deterministic_data(model).infvar_to_detvar
    haskey(map_dict, vref) || error("Variable $vref not used in the optimizer model.")
    return map_dict[vref]
end

function InfiniteOpt.optimizer_model_constraint(cref::GeneralConstraintRef,
                                                key::Val{:DetermData})
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
    If [`optimizer_model_variable`](@ref) and [`optimizer_model_constraint`](@ref)
    cannot be extended due to the nature of the reformulation then please refer
    to step 8 of the extension steps listed at the beginning of this section.

Furthermore, if appropriate for the given reformulation the following should be
extended:
- [`variable_supports`](@ref) to enable [`supports`](@ref supports(::InfiniteVariableRef))
- [`constraint_supports`](@ref) to enable [`supports`](@ref supports(::GeneralConstraintRef))
- [`constraint_parameter_refs`](@ref) to enable [`parameter_refs`](@ref parameter_refs(::GeneralConstraintRef))

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
