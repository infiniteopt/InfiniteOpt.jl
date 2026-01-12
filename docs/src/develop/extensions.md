```@meta
DocTestFilters = [r"≤|<=", r"≥|>=", r" == | = ", r" ∈ | in ", 
                  r"MathOptInterface|MOI", r" for all | ∀ ", r"d|∂", 
                  r"E|𝔼", r"integral|∫", r"-?[0-9\.]+.*"]
```

# Extensions
Here we provide guidance to various ways `InfiniteOpt` can be extended.

## Overview
Extendibility is one of the core ideas of `InfiniteOpt` so that it can serve as a 
convenient tool for those developing and implementing advanced techniques for 
infinite dimensional optimization problems. Thus, `InfiniteOpt` is developed in 
a modular manner to readily accommodate user-defined functionality and/or to 
serve as useful base in writing a `JuMP` extension. Admittedly, this modularity 
is not perfect and comments/suggestions are welcomed to help us improve this.

## Infinite Domains
Infinite domains are used to characterize the behavior of infinite parameters and 
used to govern the behavior of supports in `InfiniteOpt`. Here we walk through 
how user-defined domains can be added to various degrees of functionality. A 
template is provided in 
[`./test/extensions/infinite_domain.jl`](https://github.com/infiniteopt/InfiniteOpt.jl/blob/master/test/extensions/infinite_domain.jl). 
The extension steps employed are:
1. Define the new `struct` infinite domain type (only thing required as bare minimum)
2. Extend[`InfiniteOpt.round_domain`](@ref) (enables safe use of significant digit rounding)
3. Extend [`InfiniteOpt.supports_in_domain`](@ref) (enables error checking of supports)
4. Extend [`InfiniteOpt.generate_support_values`](@ref) (enables support generation via `num_supports` keyword arguments)
5. If a lower bound and upper bound can be reported, extend `JuMP` lower bound and upper bound methods (enables automatic bound detection in `integral`)
6. Extend [`InfiniteOpt.MeasureToolbox.generate_expect_data`](@ref) (enables the use of `expect`) 

As an example, let's create a univariate disjoint interval domain as an infinite 
domain type. This corresponds to the domain ``[lb_1, ub_1] \cup [lb_2, ub_2]`` 
where ``ub_1 \leq lb_2``. First, we need to create the `DataType` with 
inheritance from [`InfiniteScalarDomain`](@ref):
```jldoctest domain_ext; output = false
using InfiniteOpt

struct DisjointDomain <: InfiniteOpt.InfiniteScalarDomain
    lb1::Float64
    ub1::Float64
    lb2::Float64
    ub2::Float64
    # constructor
    function DisjointDomain(lb1::Number, ub1::Number, lb2::Number, ub2::Number)
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
all we have to do to add a domain in `InfiniteOpt`.

We can now define infinite parameters using this domain via 
[`@infinite_parameter`](@ref) both anonymously and explicitly:
```jldoctest domain_ext
julia> model = InfiniteModel();

julia> t = @infinite_parameter(model, domain = DisjointDomain(0, 1, 3, 4), base_name = "t")
t

julia> @infinite_parameter(model, t in DisjointDomain(0, 1, 3, 4))
t
```
Once defined (without further extension), these parameters can be used as normal 
with the following limitations:
- Supports must be specified manually (`num_supports` is not enabled)
- Supports will not be checked if they are in the domain of the infinite domain
- Domain bounds cannot be queried.
- The [`DiscreteMeasureData`](@ref) or [`FunctionalDiscreteMeasureData`](@ref) 
  must be provided explicitly to evaluate measures
However, all of these limitations except for the last one can be eliminated by 
extending a few functions as outlined above. To address the last one, we need 
to extend [`generate_integral_data`](@ref). See [`Measure Evaluation Techniques`] 
for details. 

To enable support domain checking which is useful to avoid strange bugs, we will 
extend [`InfiniteOpt.round_domain`](@ref) which rounds the domain to use proper 
significant digits and [`InfiniteOpt.supports_in_domain`](@ref) which returns a 
`Bool` whether a vector of supports is in the domain:
```jldoctest domain_ext; output = false
function InfiniteOpt.round_domain(
    domain::DisjointDomain,
    sig_digits::Int
    )
    lb1 = round(domain.lb1, sigdigits = sig_digits)
    ub1 = round(domain.ub1, sigdigits = sig_digits)
    lb2 = round(domain.lb2, sigdigits = sig_digits)
    ub2 = round(domain.ub2, sigdigits = sig_digits)
    return DisjointDomain(lb1, ub1, lb2, ub2)
end

function InfiniteOpt.supports_in_domain(
    supports::Union{Number, Vector{<:Number}},
    domain::DisjointDomain
    )
    return all((domain.lb1 .<= supports .<= domain.ub1) .| (domain.lb2 .<= supports .<= domain.ub2))
end

# output


```
Now the checks are enabled, so the following would yield an error because the 
support is not in the domain:
```jldoctest domain_ext
julia> @infinite_parameter(model, domain = DisjointDomain(0, 1, 3, 4), supports = 2)
ERROR: At none:1: `@infinite_parameter(model, domain = DisjointDomain(0, 1, 3, 4), supports = 2)`: Supports violate the domain bounds.
```

To enable automatic support generation via the `num_supports` keyword and with 
functions such as [`fill_in_supports!`](@ref), we will extend 
[`InfiniteOpt.generate_support_values`](@ref):
```jldoctest domain_ext; output = false
struct DisjointGrid <: InfiniteOpt.PublicLabel end

function InfiniteOpt.generate_support_values(
    domain::DisjointDomain;
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    sig_digits::Int = InfiniteOpt.DefaultSigDigits
    )
    length_ratio = (domain.ub1 - domain.lb1) / (domain.ub1 - domain.lb1 + domain.ub2 - domain.lb2)
    num_supports1 = Int64(ceil(length_ratio * num_supports))
    num_supports2 = num_supports - num_supports1
    supports1 = collect(range(domain.lb1, stop = domain.ub1, length = num_supports1))
    supports2 = collect(range(domain.lb2, stop = domain.ub2, length = num_supports2))
    return round.([supports1; supports2], sigdigits = sig_digits), DisjointGrid
end

# output


```
Now automatic support generation is enabled, for example:
```jldoctest domain_ext
julia> par = @infinite_parameter(model, domain = DisjointDomain(0, 2, 3, 4), num_supports = 10)
noname

julia> supports(par)
10-element Vector{Float64}:
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

We can extend the appropriate `JuMP` upper and lower bound functions 
if desired which are:
- [`JuMP.has_lower_bound`](@ref JuMP.has_lower_bound(::AbstractInfiniteDomain))
- [`JuMP.lower_bound`](@ref JuMP.lower_bound(::AbstractInfiniteDomain))
- [`JuMP.set_lower_bound`](@ref JuMP.set_lower_bound(::AbstractInfiniteDomain, ::Union{Real, Vector{<:Real}}))
- [`JuMP.has_upper_bound`](@ref JuMP.has_upper_bound(::AbstractInfiniteDomain))
- [`JuMP.upper_bound`](@ref JuMP.upper_bound(::AbstractInfiniteDomain))
- [`JuMP.set_upper_bound`](@ref JuMP.set_upper_bound(::AbstractInfiniteDomain, ::Union{Real, Vector{<:Real}}))
However, if we want `has_lower_bound = false` and `has_upper_bound = false` then 
no extension is needed. For our current example we won't do this since lower 
and upper bounds aren't exactly clear for a disjoint interval. Please refer to 
the template in `./InfiniteOpt/test/extensions/infinite_domain.jl` to see how 
this is done.

Finally, we can optionally enable the use of [`expect`](@ref) taken with respect 
to infinite parameters with this new domain type by extending 
[`InfiniteOpt.MeasureToolbox.generate_expect_data`](@ref):
```jldoctest domain_ext; output = false
function InfiniteOpt.MeasureToolbox.generate_expect_data(domain::DisjointDomain, 
    pref::GeneralVariableRef, 
    num_supports::Int; 
    kwargs...
    )
    for (k, _) in kwargs
        error("Keyword argument `$k` not supported for expectations over ",
              "disjoint domains.")
    end
    coeff_func = (supps) -> ones(size(supps)[end]) ./ size(supps)[end] 
    return InfiniteOpt.FunctionalDiscreteMeasureData(pref, coeff_func, 0, All)
end

# output


```
The above implementation simply sums over all the supports associated with `pref` 
and divides by the total number. Now we can use `expect`:
```jldoctest domain_ext
julia> @variable(model, y, Infinite(t))
y(t)

julia> expect(y, t)
𝔼{t}[y(t)]
```

## Derivative Evaluation Methods 
Derivative evaluation methods are used to dictate how we form the auxiliary 
derivative evaluation equations (derivative constraints) when we evaluate 
derivatives in InfiniteOpt. Users may wish to implement their own methods beyond 
the finite difference and orthogonal collocation ones we natively provide. Thus, 
we provide an API to do just this. A complete template is provided in 
[`./test/extensions/derivative_method.jl`](https://github.com/infiniteopt/InfiniteOpt.jl/blob/master/test/extensions/derivative_method.jl) 
to help streamline this process. The extension steps are:
1. Define the new method `struct` that inherits from the correct 
   [`AbstractDerivativeMethod`](@ref) subtype
2. Extend [`InfiniteOpt.allows_high_order_derivatives`](@ref)
3. Extend [`InfiniteOpt.generative_support_info`](@ref InfiniteOpt.generative_support_info(::AbstractDerivativeMethod)) 
   if the method is a [`GenerativeDerivativeMethod`](@ref)
4. Extend [`InfiniteOpt.derivative_expr_data`](@ref)
5. Extend [`InfiniteOpt.make_indexed_derivative_expr`](@ref).

To exemplify this process let's implement 1st order explicit Euler which is already 
implemented via `FiniteDifference(Forward())`, but let's make our own anyway for 
the sake of example. For a first order derivative ``\frac{d y(t)}{dt}`` explicit 
Euler is expressed:
```math
y(t_{n+1}) = y(t_n) + (t_{n+1} - t_{n})\frac{d y(t_n)}{dt}, \ \forall n = 0, 1, \dots, k-1
```

Let's get started with step 1 and define our new method struct:
```jldoctest deriv_ext; output = false
using InfiniteOpt

struct ExplicitEuler <: NonGenerativeDerivativeMethod end

# output


```
Notice that our method `ExplicitEuler` inherits from 
[`NonGenerativeDerivativeMethod`](@ref) since explicit Euler uses the existing 
support scheme without adding any additional supports. If our desired method 
needed to add additional supports (e.g., orthogonal collocation over finite 
elements) then we would need to have used [`GenerativeDerivativeMethod`](@ref).

Now we need to decide if this method will directly support higher order derivatives. 
In this case, let's say it won't and define:
```jldoctest deriv_ext; output = false
InfiniteOpt.allows_high_order_derivatives(::ExplicitEuler) = false

# output


```
Conversely, we could set the output to `true` if we wanted to directly support higher 
order derivatives. In which case, we would need to take the order into account in steps 4 and 5.

Since, this is a `NonGenerativeDerivativeMethod` we skip step 3. This is 
however exemplified in the extension template.

For step 4, we extend [`InfiniteOpt.derivative_expr_data`](@ref). 
This function generates all the needed data to make the
expressions necessary to build the derivative evaluation equations (derivative 
constraints). We assume these relations to be of the form ``h = 0`` where ``h`` 
is a vector of expressions. Thus, mathematically ``h`` should be of the form:
```math
\begin{aligned}
&&& y(t_{2}) - y(t_{1}) - (t_{2} - t_{1})\frac{d y(t_1)}{dt} \\
&&& \vdots \\
&&& y(t_{n+1}) - y(t_n) - (t_{n+1} - t_{n})\frac{d y(t_n)}{dt} \\
\end{aligned}
```
The required data must include the support indices used for each derivative 
variable and then any other constants needed. In this case, we will need the
indices ``\{1, \dots, n\}`` and no additional data (additional data is exemplified 
in the extension template). With this in mind let's now extend 
`InfiniteOpt.derivative_expr_data`:
```jldoctest deriv_ext; output = false
function InfiniteOpt.derivative_expr_data(
    dref::GeneralVariableRef, 
    order::Int,
    supps::Vector{Float64},
    method::ExplicitEuler
    )
    # generate the support indices to be used for each call of `make_indexed_derivative_expr`
    idxs = 1:length(supps)-1
    # return the indexes and the other iterators
    return (idxs, ) # output must be a tuple
end

# output


```

Finally, we just need to extend [`InfiniteOpt.make_indexed_derivative_expr`](@ref).
This will be used to create derivative expressions for each index determined (and additional datum) produced by `derivative_expr_data`.
```jldoctest deriv_ext; output = false
function InfiniteOpt.make_indexed_derivative_expr(
    dref::GeneralVariableRef,
    vref::GeneralVariableRef,
    pref::GeneralVariableRef,
    order::Int,
    idx,
    supps::Vector{Float64}, # ordered
    write_model::Union{InfiniteModel, AbstractTransformationBackend},
    ::ExplicitEuler,
    # put extra data args here (none in this case)
    )
    # generate the derivative expression h corresponding to the equation of 
    # the form h = 0
    d = InfiniteOpt.make_reduced_expr(dref, pref, supps[idx], write_model)
    v1 = InfiniteOpt.make_reduced_expr(vref, pref, supps[idx], write_model)
    v2 = InfiniteOpt.make_reduced_expr(vref, pref, supps[idx + 1], write_model)
    return JuMP.@expression(write_model, -(supps[idx+1] - supps[idx]) * d + v2 - v1)
end

# output


```
We used [`InfiniteOpt.make_reduced_expr`](@ref) as a convenient helper function 
to generate the semi-infinite variables/expressions we need to generate at each 
support point.

!!! note
    If your new derivative method is not compatible can not be broken 
    up into the `derivative_expr_data`-`make_indexed_derivative_expr` 
    workflow, then you can instead extend [`InfiniteOpt.evaluate_derivative`](@ref).
    This is discouraged where possible since it may make your method incompatible 
    with backends that depend on the preferred workflow.

Now that we have completed all the necessary steps, let's try it out! 
```jldoctest deriv_ext
julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10], num_supports = 3, 
                           derivative_method = ExplicitEuler());

julia> @variable(model, y, Infinite(t));

julia> dy = deriv(y, t);

julia> evaluate(dy)

julia> derivative_constraints(dy)
2-element Vector{InfOptConstraintRef}:
 -5 d/dt[y(t)](0) + y(5) - y(0) = 0
 -5 d/dt[y(t)](5) + y(10) - y(5) = 0
```
We implemented explicit Euler and it works! Now go and extend away!

## Measure Evaluation Techniques
Measure evaluation methods are used to dictate how to evaluate measures. Users 
may wish to apply evaluation methods other than Monte Carlo sampling and/or 
Gaussian quadrature methods. To create multiple measures using the same new 
evaluation methods, users may want to embed the new evaluation method under the 
[`integral`](@ref) function that does not require explicit construction of 
[`AbstractMeasureData`](@ref).

### Creating a DiscreteMeasureData Object
The basic way to do that is to write a function that creates 
[`DiscreteMeasureData`](@ref) object and pass the object to [`measure`](@ref). 
This considers a measure approximation of the form:
```math
\sum_{i \in I} \alpha_i f(\tau_i) w(\tau_i)
```
where ``\alpha_i`` are coefficients, ``f(\cdot)`` is the expression being measured, 
``w(\cdot)`` is a weighting function, and ``i \in I`` indexes the support points. 
Let's consider defining a function that enables the definition of a 
uniform grid for a univariate infinite parameter in [`IntervalDomain`](@ref). 
This example approximation uses a uniformly spaced supports ``\tau_i`` with 
``\alpha_i = \frac{ub - lb}{|I|}``:
```jldoctest measure_eval; output = false, setup = :(using InfiniteOpt)
function uniform_grid(param, num_supports)
    lb = lower_bound(param)
    ub = upper_bound(param)
    supps = collect(LinRange(lb, ub, num_supports))
    coeffs = ones(num_supports) / num_supports * (ub - lb)
    return DiscreteMeasureData(param, coeffs, supps, lower_bound = lb, upper_bound = ub)
end

# output
uniform_grid (generic function with 1 method)

```
It is necessary to pass the infinite parameter reference since the 
construction of measure data object needs parameter information. Now let's 
apply the new `uniform_grid` function to infinite parameters in 
intervals. We consider a time parameter `t` and 2D spatial parameter `x`, and 
two variables `f(t)` and `g(x)` parameterized by `t` and `x`, respectively:
```jldoctest measure_eval
julia> m = InfiniteModel();

julia> @infinite_parameter(m, t in [0, 5]);

julia> @variable(m, y, Infinite(t));
```
Now we can use `uniform_grid` to construct a [`DiscreteMeasureData`](@ref) and 
create a measure using the measure data, as shown below:

```jldoctest measure_eval
julia> tdata = uniform_grid(t, 6);

julia> y_meas = measure(y, tdata)
measure{t ∈ [0, 5]}[y(t)]

julia> expand(y_meas)
0.8333333333333333 y(0) + 0.8333333333333333 y(1) + 0.8333333333333333 y(2) + 0.8333333333333333 y(3) + 0.8333333333333333 y(4) + 0.8333333333333333 y(5)
```

### Integral Evaluation Methods
For integrals, we can implement a new approximation method via the extension of 
[`InfiniteOpt.MeasureToolbox.generate_integral_data`](@ref). This will 
allow users to use their custom measure evaluation methods in the 
[`integral`](@ref) function that does not explicitly require a measure data 
object. A template for how such an extension is accomplished is provided in 
[`./test/extensions/measure_eval.jl`](https://github.com/infiniteopt/InfiniteOpt.jl/blob/master/test/extensions/measure_eval.jl).
In general, such an extension can be created as follows: 
1. Define a new empty `struct` (e.g. `my_new_fn`) that dispatches your function
2. Extend [`InfiniteOpt.MeasureToolbox.generate_integral_data`](@ref), 
   where `method` is of the type `my_new_fn`, and `domain` needs to be a subtype 
   of [`AbstractInfiniteDomain`](@ref) that you wish to apply the new evaluation 
   method to.
Note that this procedure can be used to generate new measure evaluation methods 
not only for existing infinite domains, but also for user-defined infinite 
domains. 

For example, an extension of 
[`InfiniteOpt.MeasureToolbox.generate_integral_data`](@ref) that implements 
uniform grid for univariate and multivariate parameters in 
[`IntervalDomain`](@ref) can be created as follows:

```jldoctest measure_eval; output = false
struct UnifGrid <: InfiniteOpt.MeasureToolbox.AbstractUnivariateMethod end

function InfiniteOpt.MeasureToolbox.generate_integral_data(
    pref::InfiniteOpt.GeneralVariableRef,
    lower_bound::Real,
    upper_bound::Real,
    method::UnifGrid;
    num_supports::Int = InfiniteOpt.DefaultNumSupports,
    weight_func::Function = InfiniteOpt.default_weight
    )
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

We create measure for `y` using the `uniform_grid` method:
```jldoctest measure_eval
julia> y_int = integral(y, t, num_supports = 6, eval_method = UnifGrid())
∫{t ∈ [0, 5]}[y(t)]

julia> expand(y_int)
0.8333333333333333 y(0) + 0.8333333333333333 y(1) + 0.8333333333333333 y(2) + 0.8333333333333333 y(3) + 0.8333333333333333 y(4) + 0.8333333333333333 y(5)
```
Here we go! We can freely use `UnifGrid` for infinite parameters residing in 
[`IntervalDomain`](@ref)s now.

## [Measure Data](@id meas_data_ext)
Measures are used to evaluate over infinite domains. Users may wish to employ 
measure abstractions that cannot be readily represented with coefficients and 
discretized supports, and thus may wish to extend `InfiniteOpt`'s 
measure framework to accommodate other paradigms. This can be accomplished by  
implementing a user-defined measure data structure that inherits from 
[`AbstractMeasureData`](@ref). A template for how such an extension is 
accomplished is provided in 
[`./test/extensions/measure_data.jl`](https://github.com/infiniteopt/InfiniteOpt.jl/blob/master/test/extensions/measure_data.jl). 
The extension steps employed are:
1. Define the new data struct inheriting from [`AbstractMeasureData`](@ref) (required)
2. Extend [`InfiniteOpt.parameter_refs`](@ref parameter_refs(::AbstractMeasureData)) (required)
3. Extend [`InfiniteOpt.expand_measure`](@ref) (required)
4. Extend [`InfiniteOpt.supports`](@ref supports(::AbstractMeasureData)) (required if parameter supports are employed in any way)
5. Extend [`InfiniteOpt.add_supports_to_parameters`](@ref) (required if parameter supports are employed in measure evaluation)
6. Extend [`InfiniteOpt.coefficients`](@ref) (useful getter method if applicable)
7. Extend [`InfiniteOpt.weight_function`](@ref) (useful getter method if applicable)
8. Extend [`InfiniteOpt.support_label`](@ref) (needed to enable deletion if supports are added.)
9. Extend [`InfiniteOpt.generative_support_info`](@ref) (Needed if the measure will cause the creation of generative supports)
8. Make simple measure constructor wrapper of [`measure`](@ref) to ease definition.

To illustrate how this process can be done, let's consider extending `InfiniteOpt` 
to include measure support for assessing the variance of random expressions. The 
variance of an expression ``f(x, \xi)`` where ``x \in \mathbb{R}^n`` are finite 
variables and ``\xi \in \mathbb{R}^m`` are random infinite parameters is defined:
```math
\mathbb{V}[f(x, \xi)] = \mathbb{E}\left[(f(x, \xi) - \mathbb{E}[f(x, \xi)])^2 \right].
```
Note, we could just accomplish this by nested use of [`expect`](@ref), but we 
implement this example to illustrate the mechanics of extension.

First, let's define our new `struct` inheriting from `AbstractMeasureData`:
```jldoctest measure_data; output = false
using InfiniteOpt, Distributions

struct DiscreteVarianceData <: AbstractMeasureData
    parameter_refs::Union{GeneralVariableRef, Vector{GeneralVariableRef}}
    supports::Vector
    label::DataType
    # constructor
    function DiscreteVarianceData(
        parameter_refs::Union{GeneralVariableRef, Array{<:GeneralVariableRef}},
        supports::Vector,
        label::DataType = InfiniteOpt.generate_unique_label()
        )
        # convert input as necessary to proper array format
        if parameter_refs isa Array
            parameter_refs = convert(Vector, parameter_refs)
            supports = [convert(Vector, arr) for arr in supports]
        end
        return new(parameter_refs, supports, label)
    end
end

# output


```

We have defined our data type, so let's extend the measure data query 
methods to enable its definition. These include:
- [`parameter_refs`](@ref parameter_refs(::AbstractMeasureData))
- [`supports`](@ref supports(::AbstractMeasureData))
- [`support_label`](@ref support_label(::AbstractMeasureData))
```jldoctest measure_data; output = false
function InfiniteOpt.parameter_refs(data::DiscreteVarianceData)
    return data.parameter_refs
end

function InfiniteOpt.supports(data::DiscreteVarianceData)
    return data.supports
end

function InfiniteOpt.support_label(data::DiscreteVarianceData)
    return data.label
end

# output


```

We also need to extend [`InfiniteOpt.add_supports_to_parameters`](@ref) 
since support points will be used for measure evaluation later:
```jldoctest measure_data; output = false
function InfiniteOpt.add_supports_to_parameters(data::DiscreteVarianceData)
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
outer approximation techniques. Our extension is now sufficiently constructed to 
allow us to define out the new variance measure via [`measure`](@ref). For 
example:
```jldoctest measure_data; setup = :(using Random; Random.seed!(42))
# Setup the infinite model
model = InfiniteModel()
@infinite_parameter(model, xi ~ Normal(), num_supports = 2) # few for simplicity
@variable(model, y, Infinite(xi))
@variable(model, z)

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
function InfiniteOpt.expand_measure(
    expr::JuMP.AbstractJuMPScalar,
    data::DiscreteVarianceData,
    write_model::Union{InfiniteModel, AbstractTransformationBackend}
    )
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
y(-0.556026876146)² + 0 z*y(-0.556026876146) - 2 y(-0.44438335711)*y(-0.556026876146) + 0 z² + 0 z*y(-0.44438335711) + y(-0.44438335711)²
```

Finally, as per recommendation let's make a wrapper method to make defining 
variance measures more convenient:
```jldoctest measure_data; output = false
function variance(
    expr::Union{JuMP.GenericAffExpr, GeneralVariableRef},
    params::Union{GeneralVariableRef, Array{GeneralVariableRef}};
    name::String = "Var", 
    num_supports::Int = 10,
    use_existing::Bool = false
    )
    # get the supports
    if use_existing
        supps = supports.(params)
    else
        supps = generate_support_values(infinite_domain(first(params)),
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
expressions. (This could be overcome by defining a place holder variable 
for `expr`.

Now let's use our constructor to repeat the above measure example:
```jldoctest measure_data
julia> expand(variance(2y + z, xi, use_existing = true))
y(-0.556026876146)² + 0 z*y(-0.556026876146) - 2 y(-0.44438335711)*y(-0.556026876146) + 0 z² + 0 z*y(-0.44438335711) + y(-0.44438335711)²
```

We have done it! Now go and extend away!

## Generative Support Information 
As discussed in the [Generative Supports](@ref gen_supp_docs) section, generative 
supports help enable measure and/or derivative evaluation techniques that require 
the creation of generative supports (e.g., orthogonal collocation). Natively, we 
provide [`UniformGenerativeInfo`](@ref) to help accomplish this which works for 
creating generative supports uniformly over finite elements as is the case for 
orthogonal collocation (note this includes scaling them as need to the size of 
each finite element). However, more complex generative support schemes can be 
enabled by defining a new concrete [`AbstractGenerativeInfo`](@ref) subtype. This 
section will detail how this can be accomplished in `InfiniteOpt`. A template for 
implementing this is provided in 
[`./test/extensions/generative_info.jl`](https://github.com/infiniteopt/InfiniteOpt.jl/blob/master/test/extensions/generative_info.jl).

A new generative support information type can be created via the following:
1. Define a concrete subtype of [`AbstractGenerativeInfo`](@ref) (required)
2. Make a unique support label that inherits [`InternalLabel`](@ref) (recommended)
3. Extend [`InfiniteOpt.support_label`](@ref) (required)
4. Extend [`InfiniteOpt.make_generative_supports`](@ref) (required).

For the sake of example, let's suppose we want to make a method that generates a 
certain amount of random supports for each finite element. First, let's define 
our struct `RandomGenerativeInfo`:
```jldoctest info_model; output = false
using InfiniteOpt, Random

struct RandomGenerativeInfo <: InfiniteOpt.AbstractGenerativeInfo
    amount::Int # amount of random supports per finite element
end

# output

```
With that done, let's define a unique support label `RandomInternal` for these 
types of supports and extend `support_label`:
```jldoctest info_model; output = false
struct RandomInternal <: InternalLabel end

function InfiniteOpt.support_label(info::RandomGenerativeInfo)
    return RandomInternal
end

# output

```
Finally, let's extend `make_generative_supports` to create a vector of the 
generative supports based on a `RandomGenerativeInfo` and the existing model 
supports which are passed in the function as input:
```jldoctest info_model; output = false
function InfiniteOpt.make_generative_supports(info::RandomGenerativeInfo, pref, supps)
    num_existing = length(supps)
    num_existing <= 1 && error("`$pref` doesn't have enough supports.")
    num_internal = info.attr
    gen_supps = Float64[]
    for i = 1:num_existing-1 
        lb = supps[i]
        ub = supps[i+1]
        append!(gen_supps, rand(num_internal) * (ub - lb) .+ lb)
    end
    return gen_supps
end

# output

```
Our extension is done and now `RandomGenerativeInfo` can be incorporated by a 
`GenerativeDerivativeMethod` we create or an `AbstractMeasureData` object of our 
choice like `FunctionalDiscreteMeasureData`. 

## [Transformation Backends](@id extend_backends)
`InfiniteOpt` provides a convenient interface and abstraction for modeling 
infinite-dimensional optimization problems. By default, `InfiniteModel`s are 
reformulated into a solvable `JuMP.Model` via `TranscriptionOpt.TranscriptionBackend` 
which discretizes the model in accordance with the infinite parameter supports. 
However, users may wish to employ some other transformation method to produce 
the transformation backend. This section will explain how this can be done in 
`InfiniteOpt`. A template for implementing this extension is provided in 
[`./test/extensions/backend.jl`](https://github.com/infiniteopt/InfiniteOpt.jl/blob/master/test/extensions/backend.jl). 
Our default sub-module `InfiniteOpt.TranscriptionOpt` also serves as a good 
example.

A new transformation approach and its corresponding transformation backend can be 
extended using the following steps:
1. Define a `mutable struct` for variable/constraint mappings and other needed info (required)
2. Define an [`AbstractTransformationBackend`](@ref) (required)
3. Extend [`Base.empty!`](@ref Base.empty!(::AbstractTransformationBackend)) for the backend (required)
4. Extend [`build_transformation_backend!`](@ref build_transformation_backend!(::InfiniteModel, ::AbstractTransformationBackend)) (required)
5. If appropriate and NOT a [`JuMPBackend`](@ref), extend the following:
    - [`transformation_model`](@ref transformation_model(::AbstractTransformationBackend))
    - [`transformation_data`](@ref transformation_data(::AbstractTransformationBackend))
    - [`JuMP.set_attribute`](@ref JuMP.set_attribute(::AbstractTransformationBackend, ::Any, ::Any)) (including the suggested attributes)
    - [`JuMP.get_attribute`](@ref JuMP.get_attribute(::AbstractTransformationBackend, ::Any)) (including the suggested attributes)
    - [`JuMP.optimize!`](@ref JuMP.optimize!(::AbstractTransformationBackend))
    - [`JuMP.set_optimizer`](@ref JuMP.set_optimizer(::AbstractTransformationBackend, ::Any))
    - [`JuMP.bridge_constraints`](@ref JuMP.bridge_constraints(::AbstractTransformationBackend))
    - [`JuMP.add_bridge`](@ref JuMP.add_bridge(::AbstractTransformationBackend, ::Any))
    - [`JuMP.print_active_bridges`](@ref JuMP.print_active_bridges(::IO,::AbstractTransformationBackend))
    - [`JuMP.print_active_bridges`](@ref JuMP.print_active_bridges(::IO,::AbstractTransformationBackend))
    - [`JuMP.compute_conflict!`](@ref JuMP.compute_conflict!(::AbstractTransformationBackend))
    - [`JuMP.copy_conflict`](@ref JuMP.copy_conflict(::AbstractTransformationBackend))
    - [`JuMP.mode`](@ref JuMP.mode(::AbstractTransformationBackend))
    - [`JuMP.backend`](@ref JuMP.backend(::AbstractTransformationBackend))
    - [`JuMP.unsafe_backend`](@ref JuMP.unsafe_backend(::AbstractTransformationBackend))
    - [`JuMP.variable_ref_type`](@ref JuMP.variable_ref_type(::AbstractTransformationBackend))
6. Extend the following, if possible (also enables result queries for `JuMPBackend`s):
    - [`transformation_variable`](@ref transformation_variable(::GeneralVariableRef, ::AbstractTransformationBackend))
    - [`transformation_expression`](@ref transformation_expression(::Any, ::AbstractTransformationBackend))
    - [`transformation_constraint`](@ref transformation_constraint(::InfOptConstraintRef, ::AbstractTransformationBackend))
7. Extend the following, if appropriate:
    - [`InfiniteOpt.variable_supports`](@ref)
    - [`InfiniteOpt.expression_supports`](@ref)
    - [`InfiniteOpt.constraint_supports`](@ref)
8. As appropriate and if NOT a `JuMPBackend`, extend the following:
    - The remaining result related attributes listed in [`JuMP.get_attribute`](@ref JuMP.get_attribute(::AbstractTransformationBackend, ::Any))
    - [`JuMP.lp_sensitivity_report`](@ref JuMP.lp_sensitivity_report(::AbstractTransformationBackend))
    - [`InfiniteOpt.warmstart_backend_start_values`](@ref) to enable automated warmstarts
9. If Step 6 was skipped and/or the backend is NOT a `JuMPBackend` then extend the following:
    - [`InfiniteOpt.map_value`](@ref) (enables `JuMP.value`)
    - [`InfiniteOpt.map_infinite_parameter_value`](@ref) (enables `JuMP.value` for infinite parameters)
    - [`InfiniteOpt.map_optimizer_index`](@ref) (enables `JuMP.optimizer_index`)
    - [`InfiniteOpt.map_reduced_cost`](@ref) (enables `JuMP.reduced_cost`)
    - [`InfiniteOpt.map_shadow_price`](@ref) (enables `JuMP.shadow_cost`)
    - [`InfiniteOpt.map_dual`](@ref) (enables `JuMP.dual`)
10. Extend [`InfiniteOpt.add_point_variable`](@ref) and 
    [`InfiniteOpt.add_semi_infinite_variable`](@ref) to use 
    [`expand_measure`](@ref) without modifying the infinite model.
11. Extend [`InfiniteOpt.update_parameter_value`](@ref) to enable incremental parameter updates.
12. Extend [`map_value_to_start`](@ref) if the transformation does not return discretized values.

This may seem like a lot a work, but the majority of the above steps can be
skipped for [`JuMPBackend`](@ref)s as exemplified below. A complete extension,
showing all the above is provided in the extension template file.

For the sake of example, let's suppose we want to define a reformulation method 
for `InfiniteModel`s that are 2-stage stochastic programs (i.e., only 
`DistributionDomain`s are used, infinite variables are random 2nd stage variables, 
and finite variables are 1st stage variables). In particular, let's make a simple 
method that replaces the infinite parameters with their mean values, giving us 
the deterministic mean-valued problem.

First, let's define the (potentially mutable) `struct` that will be used to store 
our variable and constraint mappings. This case it is quite simple since our 
deterministic model will have a 1-to-1 mapping:
```jldoctest opt_model; output = false
using InfiniteOpt, Distributions

struct DeterministicData
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

Now let's define the transformation backend based on [`JuMPBackend`](@ref)
that will use a tag `Deterministic`:
```jldoctest opt_model; output = false
struct Deterministic <: AbstractJuMPTag end

const DeterministicBackend = JuMPBackend{Deterministic, Float64, DeterministicData}

# Constructor
function DeterministicBackend(; kwargs...)
    return JuMPBackend{Deterministic}(Model(; kwargs...), DeterministicData())
end
function DeterministicBackend(optimizer_constructor; kwargs...)
    backend = DeterministicBackend(; kwargs...)
    set_optimizer(backend.model, optimizer_constructor)
    return backend
end

# output
JuMPBackend{Deterministic, Float64, DeterministicData}

```
With the constructor we can now specify that a given `InfiniteModel` uses a 
`DeterministicBackend` instead of a `TranscriptionBackend` or via 
[`set_transformation_backend`](@ref):
```jldoctest opt_model; output = false
using Ipopt

# Make model using Ipopt and DeterministicModels
dbackend = DeterministicBackend(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
model = InfiniteModel(dbackend)

# Or equivalently
model = InfiniteModel()
set_transformation_backend(model, DeterministicBackend())
set_optimizer(model, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))

# output

```
Now `model` uses a `DeterministicBackend` as its transformation backend! With that we can
build our `InfiniteModel` as normal, for example:
```jldoctest opt_model
@infinite_parameter(model, ξ ~ Uniform())
@variable(model, y[1:2] >= 0, Infinite(ξ))
@variable(model, z)
@objective(model, Min, z + expect(y[1] + y[2], ξ))
@constraint(model, 2y[1] - z <= 42)
@constraint(model, y[2]^2 + ξ == 2)
@constraint(model, sin(z) >= -1)
print(model)

# output
Min z + 𝔼{ξ}[y[1](ξ) + y[2](ξ)]
Subject to
 y[1](ξ) ≥ 0.0, ∀ ξ ~ Uniform
 y[2](ξ) ≥ 0.0, ∀ ξ ~ Uniform
 2 y[1](ξ) - z ≤ 42.0, ∀ ξ ~ Uniform
 y[2](ξ)² + ξ = 2.0, ∀ ξ ~ Uniform
 sin(z) - -1 ≥ 0.0
```

We have defined our `InfiniteModel`, but now we need to specify how to 
reformulate it into a `DeterministicBackend`. This is accomplished by extending 
[`build_transformation_backend!`](@ref build_transformation_backend!(::InfiniteModel, ::AbstractTransformationBackend))
which will enable the use of `optimize!`. A necessary preliminary step though, is to
define `Base.empty!` for `DeterministicData`:
```jldoctest opt_model; output = false
function Base.empty!(data::DeterministicData)
    empty!(data.infvar_to_detvar)
    empty!(data.infconstr_to_detconstr)
    return data
end

# output

```
This enables the backend to be cleared out before it is rebuilt which is
necessary to allow for modifications to the model.
Now, let's define an internal function `_make_expression` that will use dispatch to 
convert an `InfiniteOpt` expression into a `JuMP` expression using the mappings 
stored in `backend`'s `DeterministicData`:
```jldoctest opt_model; output = false
## Make dispatch methods for converting InfiniteOpt expressions
# GeneralVariableRef
function _make_expression(backend::DeterministicBackend, expr::GeneralVariableRef)
    return _make_expression(backend, expr, index(expr))
end
# IndependentParameterRef
function _make_expression(
    backend::DeterministicBackend, 
    expr::GeneralVariableRef, 
    ::IndependentParameterIndex
    )
    return mean(infinite_domain(expr).distribution) # assuming univariate
end
# FiniteParameterRef
function _make_expression(
    backend::DeterministicBackend, 
    expr::GeneralVariableRef, 
    ::FiniteParameterIndex
    )
    return parameter_value(expr)
end
# DependentParameterRef
function _make_expression(
    backend::DeterministicBackend, 
    expr::GeneralVariableRef, 
    ::DependentParameterIndex
    )
    return mean(infinite_domain(expr).distribution) # assuming valid distribution
end
# DecisionVariableRef
function _make_expression(
    backend::DeterministicBackend, 
    expr::GeneralVariableRef, 
    ::Union{InfiniteVariableIndex, FiniteVariableIndex}
    )
    return backend.data.infvar_to_detvar[expr]
end
# MeasureRef --> assume is expectation
function _make_expression(
    backend::DeterministicBackend, 
    expr::GeneralVariableRef,
    ::MeasureIndex
    )
    return _make_expression(backend, measure_function(expr))
end
# AffExpr/QuadExpr/NonlinearExpr
function _make_expression(backend::DeterministicBackend, expr::Union{GenericAffExpr, GenericQuadExpr, GenericNonlinearExpr})
    return map_expression(v -> _make_expression(backend, v), expr)
end

# output
_make_expression (generic function with 7 methods)

```
For simplicity in example, above we assume that only `DistributionDomain`s are 
used, there are not any `PointVariableRef`s, and all `MeasureRef`s correspond to 
expectations. Naturally, a full extension should include checks to enforce that 
such assumptions hold. Notice that [`map_expression`](@ref) is useful for 
converting expressions.

Now let's extend
[`build_transformation_backend!`](@ref build_transformation_backend!(::InfiniteModel, ::AbstractTransformationBackend))
for `DeterministicBackend`s. This should build out the backend in-place and thus we should also be
sure to have it clear out any previous builds with `Base.empty!`:
```jldoctest opt_model; output = false
function InfiniteOpt.build_transformation_backend!(
    model::InfiniteModel,
    backend::DeterministicBackend
    )
    # TODO check that `model` is a stochastic model
    # empty the model for a build/rebuild
    empty!(backend)
    backend.model.operator_counter = 0 # clears out any previous user defined operators

    # add user-defined nonlinear operators if there are any
    add_operators_to_jump(backend.model, model)

    # add variables
    for vref in all_variables(model)
        if index(vref) isa InfiniteVariableIndex
            start = NaN # simple hack for sake of example
        else
            start = start_value(vref)
            start = isnothing(start) ? NaN : start
        end
        lb = has_lower_bound(vref) ? lower_bound(vref) : NaN
        ub = has_upper_bound(vref) ? upper_bound(vref) : NaN
        if is_fixed(vref)
            lb = fix_value(vref)
        end
        info = VariableInfo(!isnan(lb), lb, !isnan(ub), ub, is_fixed(vref), lb, 
                            !isnan(start), start, is_binary(vref), is_integer(vref))
        new_vref = add_variable(backend.model, ScalarVariable(info), name(vref))
        backend.data.infvar_to_detvar[vref] = new_vref
    end

    # add the objective
    obj_func = _make_expression(backend, objective_function(model))
    set_objective(backend.model, objective_sense(model), obj_func)

    # add the constraints
    for cref in all_constraints(model, Union{GenericAffExpr, GenericQuadExpr, GenericNonlinearExpr})
        constr = constraint_object(cref)
        new_func = _make_expression(backend, constr.func)
        new_constr = build_constraint(error, new_func, constr.set)
        new_cref = add_constraint(backend.model, new_constr, name(cref))
        backend.data.infconstr_to_detconstr[cref] = new_cref
    end
    return
end

# output

```
Note that Step 5 can be skipped since we are using the `JuMPBackend` API which inherits 
all the needed methods. Now we can build our backend automatically and enable the use of 
`optimize!`:
```jldoctest opt_model
optimize!(model)
print(transformation_model(model))

# output
Min z + y[1] + y[2]
Subject to
 sin(z) - -1.0 ≥ 0
 2 y[1] - z ≤ 42
 y[2]² = 1.5
 y[1] ≥ 0
 y[2] ≥ 0
```
Note that better variable naming could be used with the reformulated infinite 
variables. Moreover, in general extensions of [`build_transformation_backend!`](@ref) 
should account for the possibility that `InfiniteModel` contains constraints with 
[`DomainRestriction`](@ref) as accessed via [`domain_restriction`](@ref).

Now that we have optimized our `InfiniteModel` via the use the of a 
`DeterministicBackend`, we probably want to access the results. All queries 
are enabled via Step 6 where we extend:
- [`transformation_variable`](@ref transformation_variable(::GeneralVariableRef, ::AbstractTransformationBackend))
- [`transformation_expression`](@ref transformation_expression(::Any, ::AbstractTransformationBackend))
- [`transformation_constraint`](@ref transformation_constraint(::InfOptConstraintRef, ::AbstractTransformationBackend))
to return the variable(s)/expression(s)/constraint(s) in the backend.
These will use the `DeterministicData` and should error if no mapping can be 
found.
```jldoctest opt_model; output = false
function InfiniteOpt.transformation_variable(
    vref::GeneralVariableRef,
    backend::DeterministicBackend
    )
    map_dict = backend.data.infvar_to_detvar
    haskey(map_dict, vref) || error("Variable $vref not used in the transformation backend.")
    return map_dict[vref]
end

function InfiniteOpt.transformation_expression(
    expr::JuMP.AbstractJuMPScalar,
    backend::DeterministicBackend
    )
    return _make_expression(backend, expr)
end

function InfiniteOpt.transformation_constraint(
    cref::InfOptConstraintRef,
    backend::DeterministicBackend
    )
    map_dict = backend.data.infconstr_to_detconstr
    haskey(map_dict, cref) || error("Constraint $cref not used in the transformation backend.")
    return map_dict[cref]
end

# output

```
With these extensions we can now access all the result queries (skipping Steps 8-9). For example:
```jldoctest opt_model
julia> termination_status(model)
LOCALLY_SOLVED::TerminationStatusCode = 4

julia> result_count(model)
1

julia> value.(y)
2-element Vector{Float64}:
 -9.164638781941642e-9
  1.224744871391589

julia> optimizer_index(z)
MathOptInterface.VariableIndex(3)
```
We also skip steps 7 and 10 since these are not applicable to this particular 
example.

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

Please reach out to us via the 
[discussion forum](https://github.com/infiniteopt/InfiniteOpt.jl/discussions) to 
discuss your plans before starting this on your own.
