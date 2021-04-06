# [Measure Operators] (@id measure_page)
A guide and manual for the definition and use of measures in `InfiniteOpt`.
The Datatypes and Methods sections at the end comprise the manual, and the
above sections comprise the guide.  

## Overview
Measure operators are objects that capture the evaluation of an expression with respect
to parameters, which is a distinct feature of optimization problems with
infinite decision spaces. In dynamic optimization measures can represent integral
terms such as the total cost over time, and in stochastic optimization measures
can represent integrals over the uncertain parameters, such as expectations. In
`InfiniteOpt`, measures are general operators that can be uni-variate or 
multi-variate. Natively we employ measure abstractions that employ discretization 
schemes, which evaluate the expression at a set of points over the parameter space and
approximates the measures based on the expression values at these points. However, 
we support the use of alternative measure operator paradigms.

## [Basic Usage] (@id measure_basic_usage)
First, we consider a dynamic optimization problem with the time parameter `t`
from 0 to 10. We also consider a state variable `y(t)` and a control variable
`u(t)` that are parameterized by `t`:
```jldoctest meas_basic; setup = :(using InfiniteOpt, JuMP, Random; Random.seed!(0); model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10])
t

julia> @infinite_variable(model, y(t))
y(t)

julia> @infinite_variable(model, u(t))
u(t)
```

Now suppose we want to evaluate the integral ``\int_{2}^{8}y(t)^2 + u(t)^2 dt``.
We can construct a measure to represent this integral using the
[`integral`](@ref) function
```jldoctest meas_basic
julia> mref1 = integral(y^2 + u^2, t, 2, 8)
∫{t ∈ [2, 8]}[y(t)² + u(t)²]
```
The four positional arguments of [`integral`](@ref) are the integrand expression, 
the parameter of integration, the lower bound, and the upper bound, respectively. 
Specifying the integrand expression and the parameter of integration is required.
If the lower and upper bounds are not specified, then the integration will
be over the entire domain, which is ``[0, 10]`` in this case.

The `integral` function uses trapezoid rule as the default discretization scheme
for univariate parameters in finite `IntervalSet`s. In addition, the user can also 
use quadrature methods for univariate parameters in all `IntervalSet`s by setting
the keyword argument `eval_method` as `Quadrature()`:
```jldoctest meas_basic
julia> mref2 = integral(y^2 + u^2, t, eval_method = Quadrature())
∫{t ∈ [0, 10]}[y(t)² + u(t)²]
```

The `integral` function also allows for specifying the number of points for the
discretization scheme using the keyword argument `num_supports`. The default
value of `num_supports` is 10.
```jldoctest meas_basic
julia> mref3 = ∫(y^2 + u^2, t, num_supports = 20)
∫{t ∈ [0, 10]}[y(t)² + u(t)²]
```
Notice here how we used [`∫`](@ref) in place of `integral` as a convenient wrapper.

Two other explicit measure type methods include [`expect`](@ref) for expectations 
and [`support_sum`](@ref) for summing an expression over the support points of 
selected infinite parameters. The syntax for these is analogous to that of `integral` 
except that there are no lower/upper bounds. For example, we can define the following 
expectation of a random expression:
```jldoctest; setup = :(using InfiniteOpt, JuMP, Distributions)
julia> m = InfiniteModel();

julia> @infinite_parameter(m, ξ in Normal(), num_supports = 100);

julia> @infinite_variable(m, x(ξ));

julia> expect_x = expect(x^2, ξ)
𝔼{ξ}[x(ξ)²]
```

!!! note
    For integrals, expectations, and support sums involving moderate to large 
    expressions, the macro versions [`@integral`](@ref), [`@expect`](@ref), and 
    [`@support_sum`](@ref) should be used instead of their functional equivalents 
    for better performance.

!!! note 
    For convenience in compact representation we can use [`∫`](@ref), [`@∫`](@ref), 
    [`𝔼`](@ref), and [`@𝔼`](@ref) as wrappers for [`integral`](@ref), 
    [`@integral`](@ref), [`expect`](@ref), and [`@expect`](@ref), respectively.

Other measure paradigms can be implemented via [`measure`](@ref) as described in 
the sections further below.

Depending on the type of measures created, support points may be generated
at the time of creating the measures. In these cases, the new support points
will be added to the support list of the integrated parameter.

Once a measure is created, the evaluation of that measure is stored in a
measure data object. Users can query the measure data object using the
[`measure_data`](@ref) function as follows
```jldoctest meas_basic
julia> measure_data(mref2)
DiscreteMeasureData{GeneralVariableRef,1,Float64}(t, [0.333357, 0.747257, 1.09543, 1.34633, 1.47762, 1.47762, 1.34633, 1.09543, 0.747257, 0.333357], [0.130467, 0.674683, 1.60295, 2.83302, 4.25563, 5.74437, 7.16698, 8.39705, 9.32532, 9.86953], UniqueMeasure{Symbol("##809")}, InfiniteOpt.default_weight, 0.0, 10.0, false)

julia> measure_data(mref3)
FunctionalDiscreteMeasureData{GeneralVariableRef,Float64,NoGenerativeSupports}(t, InfiniteOpt.MeasureToolbox._trapezoid_coeff, 0, All, NoGenerativeSupports(), InfiniteOpt.default_weight, 0.0, 10.0, false)
```
Natively in `InfiniteOpt`, two types of measure data objects are used to store the measure
data information depending on the nature of the measures created: `DiscreteMeasureData` and
`FunctionalDiscreteMeasureData`. For more details on the measure data object, 
refer to [Measure Data Generation](@ref).

Similarly, one can also query the expression the measure operates on using 
[`measure_function`](@ref):
```jldoctest meas_basic
julia> measure_function(mref3)
y(t)² + u(t)²
```

In addition to `eval_method` and `num_supports` as shown above, `integral` function 
also accepts `weight_func` as keyword argument, which dictates the weight function
of the measure. The default value of these keyword arguments can be queried using
[`uni_integral_defaults`](@ref) and [`multi_integral_defaults`](@ref) as follows:
```jldoctest meas_basic
julia> uni_integral_defaults()
Dict{Symbol,Any} with 3 entries:
  :num_supports => 10
  :eval_method  => Automatic()
  :weight_func  => default_weight

julia> multi_integral_defaults()
Dict{Symbol,Any} with 3 entries:
  :num_supports => 10
  :eval_method  => Automatic()
  :weight_func  => default_weight
```
`Automatic` dictates that the integral is created using the default method depending
on the type of integral, and `default_weight` is assigning weights of 1 for all points.

Now suppose we want to create multiple measures that share the same keyword argument 
values that are different from the defaults. We don't have to input the keyword argument
values every time we construct a new measure. Instead, we can modify the
default values of measure keyword arguments, and construct measures using the new 
default values. To do that, we use the functions
[`set_uni_integral_defaults`](@ref) and [`set_multi_integral_defaults`](@ref).
Adding new keyword arguments will be useful if users want to extend the measure 
functions with their custom representation/evaluation schemes that need to take 
additional arguments somehow. See [Extensions](@ref) for more details.

Now we can add integrals to the constraints and objective functions in our
model using these measures. For more detailed information, please review the
information below.

## Theoretical Abstraction
In `InfiniteOpt`, measures denote operators ``M_\ell`` that operate on some infinite 
expression ``y`` over the infinite domain ``\mathcal{D}_\ell`` associated with 
the infinite parameter ``\ell``:
```math
M_{\ell}y : \mathcal{D}_{-\ell} \mapsto \mathbb{R}^{n_y}
```
Such a paradigm can capture a wide variety of mathematical operators commonly 
encountered in infinite-dimensional programming such as integrals, expectations, 
risk measures, and chance constraints.

Currently, `InfiniteOpt` natively contains programmatic objects for measures that 
can be represented as integrals of the form:
```math
\int_{\tau \in \mathcal{T}} f(\tau)w(\tau) d\tau
```
where ``\tau`` is a (possibly multivariate) infinite parameter, ``f(\tau)`` is an 
expression parameterized by ``\tau``, ``w(\tau)`` is a weight function, and 
``\mathcal{T}`` is a subset of the domain of ``\tau``. The measures approximate 
the integrals by taking a discretization scheme
```math
\int_{\tau \in \mathcal{T}} f(\tau)w(\tau) d\tau \approx \sum_{i=1}^N \alpha_i f(\tau_i) w(\tau_i)
```
where ``\tau_i`` are the grid points where the expression ``f(\tau)`` is
evaluated, and ``N`` is the total number of points taken.

This is the abstraction behind both [`DiscreteMeasureData`](@ref) and 
[`FunctionalDiscreteMeasureData`](@ref) which are the native measure data types 
in InfiniteOpt. The [Measure Data Generation](@ref) section below details how 
these can be implemented to enable schemes that fit this mathematical paradigm, but 
lie out of the realm of the supported features behind `integral`, `expect`, and 
`support_sum`.

More complex measure paradigms can also be implemented by creating concrete 
subtype of [`AbstractMeasureData`](@ref) as detailed in [Measure Data](@ref) Section 
on our extensions page.

## Measure Data Generation
The general [`measure`](@ref) function takes two arguments: the argument expression and
a measure data object that contains the details of the measure representation.
Measure data objects can be constructed using [`DiscreteMeasureData`](@ref),
where the parameter of integration, the coefficients ``\alpha_i``, and the
support points need to be defined explicitly. For example, if we want to
evaluate a function at each integer time point between 0 and 10, we
can construct the following measure data object to record this discretization
scheme:
```jldoctest meas_basic
julia> md_t = DiscreteMeasureData(t, ones(10), [i for i in 1:10])
DiscreteMeasureData{GeneralVariableRef,1,Float64}(t, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], UniqueMeasure{Symbol("##813")}, InfiniteOpt.default_weight, NaN, NaN, false)
```
The arguments of [`DiscreteMeasureData`](@ref) are parameter, coefficients, and
supports. The default weight function is ``w(\tau) = 1`` for
any ``\tau``, which can be overwritten by the keyword argument `weight_function`.
The `weight_function` should take a function that returns a number for any
value that is well defined for the integrated infinite parameter. The data type
is [`DiscreteMeasureData`](@ref), which is a subtype of the abstract data type
[`AbstractMeasureData`](@ref).

With [`DiscreteMeasureData`](@ref), a measure can be generated in a custom and
quick manner. For example, using the measure data above, we can define a measure
for ``y^2`` as follows:
```jldoctest meas_basic
julia> mref = measure(y^2, md_t)
measure{t}[y(t)²]
```
In the same way, we can define measure data for multi-variate infinite parameters.
For example, we can define a discretization scheme for a 2D position parameter
``x \in [0, 1] \times [0, 1]`` as follows:
```jldoctest meas_basic
julia> @infinite_parameter(model, x[1:2] in [0, 1])
2-element Array{GeneralVariableRef,1}:
 x[1]
 x[2]

julia> md_x = DiscreteMeasureData(x, 0.25 * ones(4), [[0.25, 0.25], [0.25, 0.75], [0.75, 0.25], [0.75, 0.75]])
DiscreteMeasureData{Array{GeneralVariableRef,1},2,Array{Float64,1}}(GeneralVariableRef[x[1], x[2]], [0.25, 0.25, 0.25, 0.25], [0.25 0.25 0.75 0.75; 0.25 0.75 0.25 0.75], UniqueMeasure{Symbol("##818")}, InfiniteOpt.default_weight, [NaN, NaN], [NaN, NaN], false)
```
where `md_x` cuts the domain into four 0.5-by-0.5 squares, and evaluates the
integrand on the center of these squares. Note that for multivariate parameters, 
each support point should be an `AbstractArray` that stores the value at each dimension.

In addition to the intuitive [`DiscreteMeasureData`], another type of measure data 
object is [`FunctionalDiscreteMeasureData`](@ref). This type captures measure data
where the support points are not known at the time of measure data creation. Instead of 
storing the specific support and coefficient values, `FunctionalDiscreteMeasureData`
stores the minimum number of supports required for the measure, and a coefficient function
that maps supports to coefficients. When the measure is built on a `FunctionalDiscreteMeasureData` 
is evaluated (expanded), supports will be generated based on the functions stored in 
the data object. The method of support generation is recorded as a `label` in the
measure object. 

For example, suppose we want to uniformly generate at least 20 Monte Carlo samples 
over the interval that `t` is in. A build-in label `UniformGrid` can be used to 
signify the use of this method. A `FunctionalDiscreteMeasureData` can be created as follows:
```jldoctest meas_basic
julia> coeff_f(supports) = [(10 - 0) / length(supports) for i in supports]
coeff_f (generic function with 1 method)

julia> fmd_t = FunctionalDiscreteMeasureData(t, coeff_f, 20, UniformGrid)
FunctionalDiscreteMeasureData{GeneralVariableRef,Float64,NoGenerativeSupports}(t, coeff_f, 20, UniformGrid, NoGenerativeSupports(), InfiniteOpt.default_weight, NaN, NaN, false)
```
For more details see [`FunctionalDiscreteMeasureData`](@ref). 

Our higher-level measure methods, such as [`integral`](@ref), do not require 
explicit construction of the measure data object and instead serve as wrappers 
that construct the appropriate data object and then call [`measure`](@ref).

## Evaluation Methods
The [`integral`](@ref) function calls [`generate_integral_data`](@ref) under the hood
to construct the measure data object. [`generate_integral_data`](@ref) takes as
positional arguments the integrated parameter, lower bound, upper bound, and method, and returns
a measure data object of type [`AbstractMeasureData`](@ref).

[`generate_integral_data`](@ref) applies multiple dispatch to encode different
support generation methods depending on the input `eval_method`. Each dispatch is distinguished by 
the `method`, which takes a concrete subtype of `AbstractIntegralMethod`. 
Each dispatch of `generate_integral_data` implements the specified method and returns
the resulting measure data, which will be used by [`@integral`](@ref) to create the measure.
A table of available `method` options in our package is listed below.
Each method is limited on the dimension of parameter and/or the type of set
that it can apply for. For the details of what each method type means, refer to the corresponding
docstrings.

The eval method `FEGaussLobatto` creates finite elements by decomposing over the supports that have been added to the integral parameter and then approximates the integral over each finite element via Lobatto quadrature using num_nodes. 
All other Gauss quadrature methods do not incorporate any existing supports nor do they decompose the integral into finite elements, 
but instead generate their quadrature node points over the entire integral domain. 
See the A Note on Support Management Section for more information."
This method will take in the user supports, and create generative supports along each interval 
and match them with corresponding coefficients. An example of what this looks like can be seen below. 

![Image](C:\Users\david\.julia\dev\InfiniteOpt\docs\src\assets\FEGaussLobatto2.png)

``\int_{x_1}^{x_3} f(x) dx = \int_{x_1}^{x_2} f(x) dx + \int_{x_2}^{x_3} f(x) dx``

``\approx \sum_{i=1}^{n} \alpha_{ia} f(\tau_{ia}) + \sum_{i=1}^{n} \alpha_{ib} f(\tau_{ib})``

`` = \sum_{i=1}^{n} (\alpha_{ia} f(\tau_{ia}) + \alpha_{ib} f(\tau_{ib}))``

where ``\tau_{ia}`` and ``\tau_{ib}`` are the discrete nodes for the two intervals 

and ``\alpha_{ia}`` and ``\alpha_{ib}`` are the coefficients.

```jldoctest meas_basic
julia> mref_lob = integral(y^2 + u^2, t, num_nodes = 3, eval_method = FEGaussLobatto())
∫{t ∈ [0, 10]}[y(t)² + u(t)²]

julia> expand(mref_lob)

```

We pass set `num_nodes` to be equal to three, to give three generative supports inbetween each of the user defined supports. 




| Evaluation Method                                | Uni/Multi-Variate? | Weight Function              | Set Type                            |
|:------------------------------------------------:|:------------------:|:----------------------------:|:-----------------------------------:|
| [`Automatic()`](@ref)                            | Both               |   ``1``                      | Any                                 |
| [`UniTrapezoid()`](@ref)                         | Both               |   ``1``                      | [`IntervalSet`](@ref)               |
| [`UniMCSampling()`](@ref)                        | Univariate         |   ``1``                      | Finite [`IntervalSet`](@ref)        |
| [`UniIndepMCSampling()`](@ref)                   | Univariate         |   ``1``                      | Finite [`IntervalSet`](@ref)        |
| [`Quadrature()`](@ref)                           | Univariate         |   ``1``                      | [`IntervalSet`](@ref)               |
| [`GaussLegendre()`](@ref)                        | Univariate         |   ``1``                      | Finite [`IntervalSet`](@ref)        |
| [`GaussRadau()`](@ref)                           | Univariate         |   ``1``                      | Finite [`IntervalSet`](@ref)        |
| [`GaussJacobi()`](@ref)                          | Univariate         | `` (1-x)^\alpha (1+x)^\beta``| Finite [`IntervalSet`](@ref)        |
| [`GaussLobatto()`](@ref)                         | Univariate         |   ``1``                      | Finite [`IntervalSet`](@ref)        |
| [`FEGaussLobatto()`](@ref)                       | Univariate         |   ``1``                      | Finite [`IntervalSet`](@ref)        |
| [`GaussChebyshev(1)`](@ref GaussChebyshev)       | Univariate         | ``\frac{1}{\sqrt{1-x^2}}  `` | Finite [`IntervalSet`](@ref)        |
| [`GaussChebyshev(2)`](@ref GaussChebyshev)       | Univariate         | ``\sqrt{1-x^2}``             | Finite [`IntervalSet`](@ref)        |
| [`GaussChebyshev(3)`](@ref GaussChebyshev)       | Univariate         | ``\sqrt{(1+x)/(1-x)}``       | Finite [`IntervalSet`](@ref)        |
| [`GaussChebyshev(4)`](@ref GaussChebyshev)       | Univariate         | ``\sqrt{(1-x)/(1+x)}``       | Finite [`IntervalSet`](@ref)        |
| [`GaussLaguerre()`](@ref)                        | Univariate         | ``e^{-x}``                   | Semi-infinite [`IntervalSet`](@ref) |
| [`GaussHermite()`](@ref)                         | Univariate         | ``e^{-x^2}``                 | Infinite [`IntervalSet`](@ref)      |
| [`MultiMCSampling()`](@ref)                      | Multivariate       |   ``1``                      | Finite [`IntervalSet`](@ref)        |
| [`MultiIndepMCSampling()`](@ref)                 | Multivariate       |   ``1``                      | Finite [`IntervalSet`](@ref)        |

In summary, we natively support trapezoid rule, Gaussian quadrature methods for univariate parameters,
and Monte Carlo sampling for both univariate and multivariate parameters.
For extension purposes, users may define their own [`generate_integral_data`](@ref)
to encode custom evaluation methods. See [Extensions](@ref) for more details.

### A Note on Support Management
There is a difference in how supports are considered using `UniTrapezoid()` and 
the other schemes. Namely, the other schemes will NOT incorporate other supports 
specified elsewhere in the model. Consider the following example with 3 equidistant 
supports and an integral objective function that uses `UniTrapezoid()` (the default):
```jldoctest support_manage; setup = :(using InfiniteOpt, JuMP), output = false
# Create a model, with one variable and an infinite parameter with a given number of supports
m = InfiniteModel()
@infinite_parameter(m, t in [0, 2], num_supports = 3)
@infinite_variable(m, u(t))

# Create an objective function with the default trapezoid integration
@objective(m, Min, integral(u^2, t))

# Get the transcribed model to check how the supports are taken into account
build_optimizer_model!(m)
trans_m = optimizer_model(m);

# output
A JuMP Model
Minimization problem with:
Variables: 3
Objective function type: GenericQuadExpr{Float64,VariableRef}
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

```

If we look at how many supports there are, how the variable `u` is transcribed, 
and how the objective function of the transcribed model looks like, we notice that 
the same supports are used in both the objective function and the transcribed 
variable:
```jldoctest support_manage
julia> supports(t) 
3-element Array{Float64,1}:
 0.0
 1.0
 2.0

julia> transcription_variable(u)  
3-element Array{VariableRef,1}:
 u(support: 1)
 u(support: 2)
 u(support: 3)

julia> objective_function(trans_m) 
0.5 u(support: 1)² + u(support: 2)² + 0.5 u(support: 3)²
```
Thus, the integral incorporates the 3 supports generated outside of the `integral` 
declaration.

Then we readjust the model to use Gauss-Legendre quadrature via `GaussLegendre()` 
that uses 2 quadrature nodes:
```jldoctest support_manage; output = false
# Set the new objective and update the TranscriptionModel
set_objective_function(m, integral(u^2, t, eval_method = GaussLegendre(), num_supports = 2))
build_optimizer_model!(m)
trans_m = optimizer_model(m);

# output
A JuMP Model
Minimization problem with:
Variables: 5
Objective function type: GenericQuadExpr{Float64,VariableRef}
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

```
Now let's look again at the number of supports, the transcription of `u`, and the 
new objective function:
```jldoctest support_manage
julia> supports(t) 
5-element Array{Float64,1}:
 0.0
 0.42264973081
 1.0
 1.57735026919
 2.0

julia> transcription_variable(u)  
5-element Array{VariableRef,1}:
 u(support: 1)
 u(support: 2)
 u(support: 3)
 u(support: 4)
 u(support: 5)

julia> objective_function(trans_m) 
u(support: 2)² + u(support: 4)²
```
The supports used in the objective function are different from the supports used 
in the transcription of `u`. The integral objective function has been transcribed 
using the 2 quadrature supports, but does not include the other supports since 
they cannot be incorporated into the Gaussian quadrature approximation. Whereas,
`u` is defined over all the supports and thus certain realizations of `u` will 
excluded from the objective function which will affect the behavior of the 
optimization and lead to unexpected results.

However, this behavior is avoided if we let the integral add the supports and 
not add supports elsewhere (for convenience we'll use `set_uni_integral_defaults`):
```jldoctest support_manage; output = false
# Define a new model, parameter, and variable
m = InfiniteModel()
@infinite_parameter(m, t in [0, 2])
@infinite_variable(m, u(t))

# Update the integral default keyword arguments for convenience 
set_uni_integral_defaults(eval_method = GaussLegendre(), num_supports = 2)

# Set the objective with our desired integral
@objective(m, Min, integral(u^2, t))

# Build the transcribed model 
build_optimizer_model!(m)
trans_m = optimizer_model(m);

# output
A JuMP Model
Minimization problem with:
Variables: 2
Objective function type: GenericQuadExpr{Float64,VariableRef}
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

```
Then we get the supports are consistent for `u` and the integral:
```jldoctest support_manage
julia> supports(t) 
2-element Array{Float64,1}:
 0.42264973081
 1.57735026919

julia> transcription_variable(u)  
2-element Array{VariableRef,1}:
 u(support: 1)
 u(support: 2)

julia> objective_function(trans_m) 
u(support: 1)² + u(support: 2)²
```
Therefore, using quadratures other than trapezoid requires careful analysis if 
there are user-defined supports in the problem. 

## Expansion
In a model, each measure records the integrand expression and an evaluation
scheme that details the discretization scheme to approximate the integral.
The model will not expand the measures until the transcription stage, at which
a `JuMP.AbstractJuMPScalar` is created for each measure to represent how
the measure is modeled in a transcription model based on the stored
discretization scheme (see [Model Transcription](@ref transcription_docs) for
details on transcription). Additional point variables will be created in the
expansion process if the measure is evaluated at infinite parameter points that
do not have corresponding point variables yet.

Sometimes for extension purposes, one might want to expand a specific measure
before reaching the transcription stage. Alternatively, one might want to use
custom reformulation instead of the transcription encoded in this package, in
which expanding measures will also be useful. This can be done using the [`expand`](@ref)
function, which takes a [`MeasureRef`](@ref) object and returns a `JuMP.AbstractJuMPScalar`
based on the [`AbstractMeasureData`](@ref). For example, suppose we want to
integrate ``y^2`` in ``t``, with two supports ``t = 2.5`` and ``t = 7.5``.
We can set up and expand this measure as follows:
```jldoctest meas_basic
julia> tdata = DiscreteMeasureData(t, [5, 5], [2.5, 7.5])
DiscreteMeasureData{GeneralVariableRef,1,Float64}(t, [5.0, 5.0], [2.5, 7.5], UniqueMeasure{Symbol("##831")}, InfiniteOpt.default_weight, NaN, NaN, false)

julia> mref4 = measure(y^2, tdata)
measure{t}[y(t)²]

julia> expanded_measure = expand(mref4)
5 y(2.5)² + 5 y(7.5)²

julia> typeof(expanded_measure)
GenericQuadExpr{Float64,GeneralVariableRef}
```
In the expand call, two point variables, `y(2.5)` and `y(7.5)`, are created
because they are not defined in the model before the expand call. One can use
the [`expand_all_measures!`](@ref) function to expand all measures in a model,
which simply applies the [`expand`](@ref) to all measures stored in the model.

## Semi-Infinite Variables
!!! note 
    Previous versions of `InfiniteOpt` referred to semi-infinite variables as 
    "reduced infinite variables". 
     
Expanding measures that cover a subset of infinite parameter dependencies present 
in an expression will introduce semi-infinite variables to the
model. To see what this means, suppose we have an infinite variable that is
parameterized by multiple infinite parameters defined as follows:
```jldoctest meas_basic
julia> @infinite_variable(model, T(x, t))
T(x, t)
```
Now say we want to integrate `T` over `t`. We can define a measure for the
integral similar to how we have defined other measures:
```jldoctest meas_basic
julia> mref5 = measure(T, tdata)
measure{t}[T(x, t)]
```
Now if we expand this measure, the measure data object `tdata` records the
supports for `t`, but no supports for `x` because `T` is not evaluated over
`x` in this measure. Therefore, point variables cannot be defined in the
measure expansion.

Instead of point variables, each new variable in the measure expansion will be
represented using semi-infinite variables. Semi-infinite variables are
reduced from their original infinite variables in that they are parameterized
by less infinite parameters. In the example above, in the expansion each
semi-infinite variable for `T` should only be parameterized by `x` since
the value of `t` is fixed. The expanded measure now looks like this:
```jldoctest meas_basic
julia> expanded_measure = expand(mref5)
5 T([x[1], x[2]], 2.5) + 5 T([x[1], x[2]], 7.5)
```
where the expanded measure is a `JuMP.GenericAffExpr` that takes in its terms [`GeneralVariableRef`](@ref)s
pointing to [`SemiInfiniteVariable`](@ref)s created on the fly. [`SemiInfiniteVariable`](@ref)
refers to the information of the semi-infinite variable stored in its model. The semi-infinite variable
records a reference for its original infinite variable, and the value of the
fixed infinite parameter. One can query this information using
[`infinite_variable_ref`](@ref) and [`eval_supports`](@ref) function as follows:
```jldoctest meas_basic
julia> T1 = first(keys(expanded_measure.terms))
T([x[1], x[2]], 2.5)

julia> infinite_variable_ref(T1)
T(x, t)

julia> eval_supports(T1)
Dict{Int64,Float64} with 1 entry:
  3 => 2.5
```
All the `JuMP` functions extended for infinite variables are also extended for
semi-infinite variables, e.g. [`JuMP.lower_bound`](@ref JuMP.lower_bound(::SemiInfiniteVariableRef)).

## Datatypes
```@index
Pages   = ["measure.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
AbstractMeasureData
DiscreteMeasureData
FunctionalDiscreteMeasureData
Measure
MeasureIndex
MeasureData
MeasureRef
SemiInfiniteVariable
SemiInfiniteVariableIndex
SemiInfiniteVariableRef
```

## Methods
```@index
Pages   = ["measure.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
default_weight
DiscreteMeasureData(::GeneralVariableRef, ::Vector{<:Real}, ::Vector{<:Real})
DiscreteMeasureData(::AbstractArray{GeneralVariableRef}, ::Vector{<:Real}, ::Vector{<:AbstractArray{<:Real}})
FunctionalDiscreteMeasureData(::GeneralVariableRef,::Function,::Int,::Type{<:AbstractSupportLabel})
FunctionalDiscreteMeasureData(::AbstractArray{GeneralVariableRef},::Function,::Int,::Type{<:AbstractSupportLabel})
parameter_refs(::AbstractMeasureData)
support_label(::AbstractMeasureData)
generative_support_info(::AbstractMeasureData)
JuMP.lower_bound(::AbstractMeasureData)
JuMP.upper_bound(::AbstractMeasureData)
supports(::AbstractMeasureData)
num_supports(::AbstractMeasureData)
min_num_supports(::AbstractMeasureData)
coefficient_function(::AbstractMeasureData)
coefficients(::AbstractMeasureData)
weight_function(::AbstractMeasureData)
build_measure
InfiniteOpt.measure_data_in_finite_var_bounds(::AbstractMeasureData,::ParameterBounds)
add_measure
InfiniteOpt.add_supports_to_parameters(::AbstractMeasureData)
measure_function
measure_data
is_analytic
parameter_refs(::MeasureRef)
measure
@measure
used_by_constraint(::MeasureRef)
used_by_measure(::MeasureRef)
used_by_objective(::MeasureRef)
is_used(::MeasureRef)
JuMP.name(::MeasureRef)
JuMP.set_name(::MeasureRef, ::String)
num_measures
all_measures
JuMP.delete(::InfiniteModel, ::MeasureRef)
expand
expand_all_measures!
InfiniteOpt.expand_measure
InfiniteOpt.analytic_expansion
InfiniteOpt.expand_measures
make_point_variable_ref
make_semi_infinite_variable_ref
add_measure_variable(::JuMP.Model, ::Any, ::Any)
delete_internal_semi_infinite_variable
delete_semi_infinite_variable(::JuMP.Model, ::Any, ::Any)
internal_semi_infinite_variable
JuMP.build_variable(::Function, ::GeneralVariableRef,::Dict{Int, Float64})
JuMP.add_variable(::InfiniteModel, ::SemiInfiniteVariable,::String)
infinite_variable_ref(::SemiInfiniteVariableRef)
eval_supports(::SemiInfiniteVariableRef)
parameter_refs(::SemiInfiniteVariableRef)
parameter_list(::SemiInfiniteVariableRef)
raw_parameter_refs(::SemiInfiniteVariableRef)
JuMP.set_name(::SemiInfiniteVariableRef,::String)
JuMP.has_lower_bound(::SemiInfiniteVariableRef)
JuMP.lower_bound(::SemiInfiniteVariableRef)
JuMP.LowerBoundRef(::SemiInfiniteVariableRef)
JuMP.has_upper_bound(::SemiInfiniteVariableRef)
JuMP.upper_bound(::SemiInfiniteVariableRef)
JuMP.UpperBoundRef(::SemiInfiniteVariableRef)
JuMP.is_fixed(::SemiInfiniteVariableRef)
JuMP.fix_value(::SemiInfiniteVariableRef)
JuMP.FixRef(::SemiInfiniteVariableRef)
start_value_function(::SemiInfiniteVariableRef)
JuMP.is_binary(::SemiInfiniteVariableRef)
JuMP.BinaryRef(::SemiInfiniteVariableRef)
JuMP.is_integer(::SemiInfiniteVariableRef)
JuMP.IntegerRef(::SemiInfiniteVariableRef)
```

## MeasureToolbox Datatypes
```@index
Pages   = ["measure.md"]
Modules = [InfiniteOpt.MeasureToolbox]
Order   = [:type]
```
```@docs
InfiniteOpt.MeasureToolbox.AbstractIntegralMethod
InfiniteOpt.MeasureToolbox.Automatic
InfiniteOpt.MeasureToolbox.AbstractUnivariateMethod
InfiniteOpt.MeasureToolbox.UniTrapezoid
InfiniteOpt.MeasureToolbox.UniMCSampling
InfiniteOpt.MeasureToolbox.UniIndepMCSampling
InfiniteOpt.MeasureToolbox.Quadrature
InfiniteOpt.MeasureToolbox.GaussHermite
InfiniteOpt.MeasureToolbox.GaussLegendre
InfiniteOpt.MeasureToolbox.GaussRadau
InfiniteOpt.MeasureToolbox.GaussLobatto
InfiniteOpt.MeasureToolbox.GaussJacobi
InfiniteOpt.MeasureToolbox.FEGaussLobatto
InfiniteOpt.MeasureToolbox.GaussChebyshev
InfiniteOpt.MeasureToolbox.GaussLaguerre
InfiniteOpt.MeasureToolbox.AbstractMultivariateMethod
InfiniteOpt.MeasureToolbox.MultiMCSampling
InfiniteOpt.MeasureToolbox.MultiIndepMCSampling
InfiniteOpt.MeasureToolbox.InternalGaussLobatto
```

## MeasureToolbox Methods
```@index
Pages   = ["measure.md"]
Modules = [InfiniteOpt.MeasureToolbox]
Order   = [:macro, :function]
```
```@docs
InfiniteOpt.MeasureToolbox.@integral
InfiniteOpt.MeasureToolbox.@∫
InfiniteOpt.MeasureToolbox.integral(::JuMP.AbstractJuMPScalar,::InfiniteOpt.GeneralVariableRef,::Real,::Real)
InfiniteOpt.MeasureToolbox.integral(::JuMP.AbstractJuMPScalar,::AbstractArray{InfiniteOpt.GeneralVariableRef},::Union{Real, AbstractArray{<:Real}},::Union{Real, AbstractArray{<:Real}})
InfiniteOpt.MeasureToolbox.∫(::JuMP.AbstractJuMPScalar,::InfiniteOpt.GeneralVariableRef,::Real,::Real)
InfiniteOpt.MeasureToolbox.∫(::JuMP.AbstractJuMPScalar,::AbstractArray{InfiniteOpt.GeneralVariableRef},::Union{Real, AbstractArray{<:Real}},::Union{Real, AbstractArray{<:Real}})
InfiniteOpt.MeasureToolbox.@expect
InfiniteOpt.MeasureToolbox.@𝔼
InfiniteOpt.MeasureToolbox.expect
InfiniteOpt.MeasureToolbox.𝔼
InfiniteOpt.MeasureToolbox.@support_sum
InfiniteOpt.MeasureToolbox.support_sum
InfiniteOpt.MeasureToolbox.uni_integral_defaults
InfiniteOpt.MeasureToolbox.set_uni_integral_defaults
InfiniteOpt.MeasureToolbox.multi_integral_defaults
InfiniteOpt.MeasureToolbox.set_multi_integral_defaults
InfiniteOpt.MeasureToolbox.generate_integral_data
```
