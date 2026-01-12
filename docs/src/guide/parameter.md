# [Infinite Parameters](@id inf_par_docs)
A guide for infinite parameters in `InfiniteOpt`. See the respective 
[technical manual](@ref inf_par_manual) for more details.

## Overview
Infinite parameters are what live on the infinite domains of infinite dimensional  
mathematical optimization problems. In dynamic optimization this 
corresponds to time and in stochastic optimization this to uncertain parameters 
that follow a certain underlying statistical distribution. `InfiniteOpt` 
considers natively two kinds of infinite parameters, ones defined over 
continuous intervals and ones characterized by a distribution (others can be 
added by defining a user-defined type). These can be used to parameterize 
infinite variables, semi-infinite variables, point variables, derivatives,  
measures, and can be used directly inside constraints. 

## Basic Usage
First, we need to initialize and add infinite parameters to our `InfiniteModel`. 
This can be accomplished using [`@infinite_parameter`](@ref). For example, let's 
define a parameter for time in a time interval from 0 to 10:
```jldoctest basic
julia> using InfiniteOpt

julia> model = InfiniteModel();

julia> @infinite_parameter(model, t in [0, 10])
t
```
Now `t` is a Julia variable that stores a [`GeneralVariableRef`](@ref) which 
points to where the time parameter is stored in `model`. It can now be used with 
infinite variables, derivatives, measures, and constraints as described in their  
respective user guide sections.

When the model is optimized, `t` will be transcribed (discretized) over its domain 
following its support points (assuming the transformation backend relies on
discretization). Users can specify support points via the  
`num_supports` or `supports` keyword arguments. For example, if we desire to 
have only 10 equidistant supports then we could have instead defined `t`: 
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel()) 
julia> @infinite_parameter(model, t in [0, 10], num_supports = 10)
t
```
More complex support schemes can be specified via `supports` such as:
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10], supports = [0, 2, 7, 10])
t
```
Where we specified `t` to use 4 supports: 0, 2, 7, and 10.

We can also add supports after `t` has been initialized. This can be accomplished 
with [`add_supports`](@ref). For example, consider the initial case where `t` has 
no supports and we now wish to add 4 supports:
```jldoctest basic
julia> add_supports(t, [0., 2.5, 7.5, 10.])

julia> supports(t)
4-element Vector{Float64}:
  0.0
  2.5
  7.5
 10.0
```
Here only 4 supports are specified for the sake of example. Alternatively, we 
could have initialized the parameter and added supports in just one step using 
the `supports` keyword argument:
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10], supports = [0., 2.5, 7.5, 10.])
t
```

We could also define a random parameter described by a distribution. This can be 
accomplished using [`@infinite_parameter`](@ref) in combination with a 
distribution from 
[`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl). For 
example let's define a vector of independent random parameters described by a 
Normal distribution:
```jldoctest basic
julia> using Distributions

julia> @infinite_parameter(model, ξ[i = 1:3] ~ Normal(), independent = true)
3-element Vector{GeneralVariableRef}:
 ξ[1]
 ξ[2]
 ξ[3]
```
Note that we use `~` instead of `in` when specifying distributions. We also could have 
used `i` as an index to assign a different distribution to each parameter. 
Supports can also be specified for each parameter as shown above. Similarly, the 
`num_supports` keyword is used to generate random supports.

More interestingly, we can also define multi-variate random parameters, for example:
```jldoctest basic
julia> @infinite_parameter(model, θ[1:2] ~ MvNormal([0, 0], [1, 1]))
2-element Vector{GeneralVariableRef}:
 θ[1]
 θ[2]
```

Now we have infinite parameters `t` and `ξ` that are ready to be used in 
defining infinite variables and constraints. We also mention here that the 
[`@infinite_parameter`](@ref) macro is designed to closely emulate 
`JuMP.@variable` and thus handles arrays and keyword arguments in the same way. 
This is described in more detail below.

## Parameter Definition
Defining/initializing an infinite parameter principally involves the following 
steps (these are typically automated by [`@infinite_parameter`](@ref)):
1. Define an [`AbstractInfiniteDomain`](@ref)
2. Define support points within the domain to later discretize the parameter
3. Construct an [`InfOptParameter`](@ref) to store this information
4. Add the `InfOptParameter` object to an `InfiniteModel` and assign a name
5. Create a [`GeneralVariableRef`](@ref)(s) that points to the parameter object

Infinite domain definition is described above in the 
[Infinite Domains](@ref infinite_domains_docs) section. The supports should be 
a vector of finite numbers that are drawn from the domain of the infinite domain. 
These supports will be used to transcribe the `InfiniteModel` in preparation for 
it to be optimized. If desired, the supports can be specified after the parameter 
is defined and the support container of the defined parameter will be temporarily 
empty.

[`InfOptParameter`](@ref) is an abstract data type that encompasses all concrete 
infinite parameter types. The concrete type for individual infinite parameters 
is [`IndependentParameter`](@ref), since these parameters are independent of 
other parameters. On the other hand, [`DependentParameters`](@ref) handle 
multivariate infinite parameters, within which each individual parameter is not 
independent. These are useful for characterizing, for example, parameters 
subject to multivariate distribution.

Regardless of the specific concrete type, the [`build_parameter`](@ref) function 
is used to construct an `InfOptParameter`. For example, let's create a time 
parameter ``t \in [0, 10]`` with supports `[0, 2, 5, 7, 10]`:
```jldoctest time_define; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> domain = IntervalDomain(0, 10)
[0, 10]

julia> t_param = build_parameter(error, domain, supports = [0, 2, 5, 7, 10]);
```  
Now that we have a `InfOptParameter` that contains an `IntervalDomain` and supports, 
let's now add `t_param` to our `InfiniteModel` using [`add_parameter`](@ref) 
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
distributions from `Distributions.jl` to define a [`UniDistributionDomain`](@ref). 
For example, let's consider a random variable ``x \in \mathcal{N}(0,1)`` with 
supports `[-0.5, 0.5]`:
```jldoctest rand_define; setup = :(using InfiniteOpt, Distributions; model = InfiniteModel())
julia> dist = Normal(0., 1.)
Normal{Float64}(μ=0.0, σ=1.0)

julia> domain = UniDistributionDomain(dist)
Normal{Float64}(μ=0.0, σ=1.0)

julia> x_param = build_parameter(error, domain, supports = [-0.5, 0.5]);
```
Again, we use [`add_parameter`](@ref) to add `x_param` to the `InfiniteModel` and 
assign it the name `x`:
```jldoctest rand_define
julia> x_ref = add_parameter(model, x_param, "x")
x
```
Note that `add_parameter` does not register the name of the parameters into the 
model that it adds to. As shown in [Macro Definition](@ref param_macro), the 
macro definition does not allow for multiple parameters sharing the same name and 
will throw an error if it happens.

For dependent parameters, we do not provide a publicly available `build_parameter`  
method due to inherent complexities. Thus, it is recommended to construct these  
using [`@infinite_parameter`](@ref). However, these can be constructed manually 
via the basic constructor for [`DependentParameters`](@ref) and then invoking 
[`add_parameters`](@ref). Note that this should be done with caution since most 
error checking will be omitted in this case.

## [Macro Definition](@id param_macro)
### One-Dimensional Parameters
One user-friendly way of defining infinite parameters is by macro 
[`@infinite_parameter`](@ref). The macro executes the same process as the 
manual definition (steps listed in [Parameter Definition](@ref)), but allows 
the users to manipulate several features of the defined infinite parameters. 
Again, let's consider a time parameter ``t \in [0, 10]`` with supports 
`[0, 2, 5, 7, 10]`. We use `in` (or `∈`) to define the domain that an infinite 
parameter is subject to (any [`InfiniteScalarDomain`](@ref) for single 
parameters). For example, we can define ``t \in [0, 10]``:
```jldoctest macro_define; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10], supports = [0, 2, 5, 7, 10])
t
```
Similarly, we can define a random infinite parameter subject to some 
distribution using `~` as the operator. For example, a Gaussian infinite 
parameter with mean 0 and standard deviation 1 can be defined:
```jldoctest macro_define
julia> using Distributions

julia> dist = Normal(0., 1.)
Normal{Float64}(μ=0.0, σ=1.0)

julia> @infinite_parameter(model, ξ ~ dist, num_supports = 10)
ξ
```

For anonymous definition, we use either the `domain` or `distribution` keywords:
```jldoctest macro_define
julia> t = @infinite_parameter(model, domain = [0, 10], supports = [0, 2, 5, 7, 10], 
                               base_name = "t")
t

julia> ξ = @infinite_parameter(model, distribution = dist, num_supports = 10, 
                               base_name = "ξ")
ξ
```

All the definitions above return a [`GeneralVariableRef`](@ref) that refers to 
the parameter object.

### Multi-Dimensional Parameters
We can also define multi-dimensional infinite parameters in a concise way. For 
example, consider a position parameter ``x \in [0, 1]^3``:
```jldoctest macro_define
julia> @infinite_parameter(model, x[1:3] in [0, 1], independent = true, num_supports = 3)
3-element Vector{GeneralVariableRef}:
 x[1]
 x[2]
 x[3]
```
Here we used `independent = true` to signify that each `x[i]` can be treated 
independently. Hence, the overall infinite domain is the cartesian product of 
their individual domains. In this example, we defined 3 supports for each `x[i]` 
such that there will be ``3^3 = 27`` supports for the overall domain:
```jldoctest macro_define
julia> supports(x[1])
3-element Vector{Float64}:
 0.0
 0.5
 1.0

julia> supports(x)
3×27 Matrix{Float64}:
 0.0  0.5  1.0  0.0  0.5  1.0  0.0  0.5  …  1.0  0.0  0.5  1.0  0.0  0.5  1.0
 0.0  0.0  0.0  0.5  0.5  0.5  1.0  1.0     0.0  0.5  0.5  0.5  1.0  1.0  1.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     1.0  1.0  1.0  1.0  1.0  1.0  1.0
```

For multivariate random infinite parameters, we commonly wish their individual 
domains to not be treated independently. Hence, we'll want `independent = false` 
(the default value). For example, a 2-dimensional parameter 
``\theta \sim \mathcal{N}(\mu, \Sigma)`` is defined:
```jldoctest macro_define; setup = :(using Random; Random.seed!(0))
julia> dist = MvNormal([0, 0], [1 0; 0 2])
FullNormal(
dim: 2
μ: [0.0, 0.0]
Σ: [1.0 0.0; 0.0 2.0]
)

julia> @infinite_parameter(model, θ[1:2] ~ dist, num_supports = 3)
2-element Vector{GeneralVariableRef}:
 θ[1]
 θ[2]
```
Here 3 supports are generated for all the parameters simultaneously according to 
the distribution. 
```jldoctest macro_define
julia> supports(θ)
2×3 Matrix{Float64}:
 0.942971  1.52507   -1.20577
 0.189395  0.175223   0.440976
```
We refer to groups of parameters defined this way as dependent infinite 
parameters. In principle, nonrandom infinite parameter types can be made 
dependent as well when users wish to avoid the Cartesian product of their 
supports.

Anonymous groups of parameters can be defined as follows:
```jldoctest macro_define
julia> x = @infinite_parameter(model, [1:3] in [0, 1], independent = true, 
                               num_supports = 3, base_name = "x")
3-element Vector{GeneralVariableRef}:
 x[1]
 x[2]
 x[3]

julia> θ = @infinite_parameter(model, [1:2] ~ dist, num_supports = 3, base_name = "θ")
2-element Vector{GeneralVariableRef}:
 θ[1]
 θ[2]
```

### Containers for Multi-Dimensional Parameters
Unlike other modelling objects, only `Array` containers are supported for infinite parameters. 

## Supports
For an infinite parameter, its supports are a finite set of points that the 
parameter will take (or possibly take, if the parameter is random). Assuming
the underlying transformation backend relies on a discretization strategy (as is
the case with [`TranscriptionBackend`](@ref)), these supports will be used
when building the transcription backend to approximate all variables/expressions
parameterized by the infinite parameter over a grid of points.

Once an infinite parameter is defined, users can access the supports using 
[`supports`](@ref) function:
```jldoctest supports; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10], supports = [0, 2, 5, 7, 10])
t

julia> supports(t)
5-element Vector{Float64}:
  0.0
  2.0
  5.0
  7.0
 10.0
```

!!! note
    Most support query functions have a keyword argument `label` that is used 
    to specify the type of supports that will be involved in the query. By default, 
    this will be `PublicLabel` which will correspond to any supports that are 
    reported to the user by default, but will exclude any supports that have 
    `InternalLabel`s (e.g., internal collocation nodes). The full set can always 
    be obtained via `label = All`. We can also query more specific subsets of 
    support information with more specific labels such as `label = UniformGrid`.

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
7-element Vector{Float64}:
  0.0
  2.0
  3.0
  5.0
  7.0
  8.0
 10.0
```
At times, we might want to change the supports completely. In those cases, the 
function [`set_supports`](@ref) resets the supports for a certain parameter with 
new supports provided:
```jldoctest supports
julia> set_supports(t, [0,3,5,8,10], force = true)

julia> supports(t)
5-element Vector{Float64}:
  0.0
  3.0
  5.0
  8.0
 10.0
```
Note that the keyword argument [`force`] must be set as [`true`] if the 
parameter has been assigned with supports. Users can also delete all the 
supports of a parameter with [`delete_supports`](@ref).

### Automatic Support Generation During Parameter Definition
For the examples in the [Parameter Definition](@ref), we have seen how to 
manually add supports to an infinite parameter. For a quick automatic 
generation of support points, though, users do not have to input the support 
points. Instead, the number of support points generated is supplied.

For an infinite parameter subject to an [`IntervalDomain`](@ref), uniformly spaced 
supports including both ends are generated across the interval. For example, 
defining a time parameter ``t \in [0, 10]`` with 4 supports using 
[`build_parameter`](@ref) gives
```jldoctest; setup = :(using InfiniteOpt)
julia> domain = IntervalDomain(0, 10)
[0, 10]

julia> t_param = build_parameter(error, domain, num_supports = 4, sig_digits = 3);
```  
Using macro definition we have
```jldoctest; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10], num_supports = 4, sig_digits = 3)
t

julia> supports(t)
4-element Vector{Float64}:
  0.0   
  3.33
  6.67
 10.0   

```
Note that the user can use the keyword argument `sig_digits` to dictate the 
significant figures for the supports. The default value of `sig_digits` is 12.

For an infinite parameter that follows a univariate distribution, 
supports are sampled from the underlying distribution. For example, we can 
define an infinite parameter subject to a normal distribution with mean 0 and 
variance 1:
```jldoctest; setup = :(using InfiniteOpt, Distributions, Random; Random.seed!(0); model = InfiniteModel())
julia> @infinite_parameter(model, x ~ Normal(), num_supports = 4)
x

julia> supports(x)
4-element Vector{Float64}:
 0.123901231206
 0.133922757653
 0.942970533446
 1.52506890851  
```
For multivariate distributions, though, we require support points are provided 
in the definition. However, we can use [`fill_in_supports!`](@ref) to generate 
supports for parameters following multivariate distributions. See 
[Automatic Support Generation For Defined Parameters](@ref) for details.

### Automatic Support Generation For Defined Parameters
So far, we have seen that in both definition methods it is allowed to initialize 
a parameter with no supports. This is done by not specifying `supports` and 
`num_supports`. However, infinite parameters would not be allowed at the 
transcription step since it needs information about how to discretize the 
infinite parameters. In previous examples, we have shown that users can add 
supports to a defined parameter using methods [`add_supports`](@ref) and 
[`set_supports`](@ref).

In this section we introduce automatic support generation for defined 
parameters with no associated supports. This can be done using the 
[`fill_in_supports!`](@ref) functions. [`fill_in_supports!`](@ref) can take as 
argument a [`GeneralVariableRef`](@ref) or an 
`Array{<:GeneralVariableRef}`, in which case it will generate 
supports for the associated infinite parameter. Alternatively, 
[`fill_in_supports!`](@ref) can also take an [`InfiniteModel`](@ref) as an 
argument, in which case it will generate supports for all infinite parameters 
of the [`InfiniteModel`](@ref) with no supports.

The [`fill_in_supports!`](@ref) method allows users to specify integer keyword 
arguments `num_supports` and `sig_digits`. `num_supports` dictates the 
number of supports to be generated, and `sig_digits` dictates the significant 
figures of generated supports desired. The default values are 10 and 12, 
respectively.

The ways by which supports are automatically generated are as follows. If the 
parameter is in an [`IntervalDomain`](@ref), then we generate an array of supports 
that are uniformly distributed along the interval, including the two ends. For 
example, consider a 3D position parameter `x` distributed in the unit cube 
`[0, 1]`. We can generate supports for that point in the following way:
```jldoctest supp_gen_defined; setup = :(using InfiniteOpt, Distributions, Random; Random.seed!(0); model = InfiniteModel())
julia> @infinite_parameter(model, x[1:3] in [0, 1], independent = true);

julia> fill_in_supports!.(x, num_supports = 3);

julia> supports.(x)
3-element Vector{Vector{Float64}}:
 [0.0, 0.5, 1.0]
 [0.0, 0.5, 1.0]
 [0.0, 0.5, 1.0]
```
Note that the dot syntax because [`fill_in_supports!`](@ref) takes single 
[`GeneralVariableRef`](@ref) as argument. In each dimension, three equally spaced 
supports (`[0.0, 0.5, 1.0]`) are generated. Since the `independent` keyword is 
set as `true`, the transcription stage will create a three-dimensional grid for 
all variables parameterized by `x`, with each point separated by 0.5 units in 
each dimension. We can view this grid by simply invoking `supports` without the  
vectorized syntax:
```jldoctest supp_gen_defined
julia> supports(x)
3×27 Matrix{Float64}:
 0.0  0.5  1.0  0.0  0.5  1.0  0.0  0.5  …  1.0  0.0  0.5  1.0  0.0  0.5  1.0
 0.0  0.0  0.0  0.5  0.5  0.5  1.0  1.0     0.0  0.5  0.5  0.5  1.0  1.0  1.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     1.0  1.0  1.0  1.0  1.0  1.0  1.0
```

If the parameter is in a [`UniDistributionDomain`](@ref) or 
[`MultiDistributionDomain`](@ref), [`fill_in_supports!`](@ref) 
samples `num_supports` supports from the distribution. Recall that support 
generation is not allowed for parameters under multivariate distribution during 
parameter definition. However, if the parameter is defined first without 
supports, [`fill_in_supports!`](@ref) allows for supports generation. For 
example, for a 2D random variable `ξ` under a multivariate Gaussian 
distribution, we can generate supports for it in the following way:
```jldoctest supp_gen_defined
julia> dist = MvNormal([0., 0.], [1. 0.; 0. 2.])
FullNormal(
dim: 2
μ: [0.0, 0.0]
Σ: [1.0 0.0; 0.0 2.0]
)


julia> @infinite_parameter(model, ξ[1:2] ~ dist);

julia> fill_in_supports!(ξ, num_supports = 3)

julia> supports(ξ)
2×3 Matrix{Float64}:
 0.942971  1.52507   -1.20577
 0.189395  0.175223   0.440976
```
Note that [`fill_in_supports!`](@ref) only fill in supports for parameters with no 
associated supports. To modify the supports of parameters already associated 
with some supports, refer to [Supports](@ref) for how to do that.

## Parameter Queries
In addition to the modeling framework, this package provides many functions for 
users to access information about the model. This section will go over basic 
functions for accessing parameter information.

Once a (possibly large-scale) `InfiniteModel` is built, the users might want to 
check if an infinite parameter is actually used in any way. This could be 
checked by [`is_used`](@ref) function as follows:
```jldoctest param_queries; setup = :(using InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, x in [0, 1])
x

julia> is_used(x)
false

```
This function checks if the parameter is used by any constraint, measure, or 
variable. Similarly, functions [`used_by_constraint`](@ref), 
[`used_by_measure`](@ref) and [`used_by_infinite_variable`](@ref) can be applied to 
find out any dependency of specific types on the infinite parameter.

In addition, sometimes we need to check if a certain [`GeneralVariableRef`](@ref) 
for an infinite parameter is valid with an `InfiniteModel` model, meaning that the 
parameter reference actually refers to some parameter associated with the model. We 
extend the [`JuMP.is_valid`](@ref) function from `JuMP` for that purpose. To see how 
to use this, for example,
```jldoctest param_queries
julia> pref1 = GeneralVariableRef(model, 1, IndependentParameterIndex);

julia> pref2 = GeneralVariableRef(model, 2, IndependentParameterIndex);

julia> is_valid(model, pref1)
true

julia> is_valid(model, pref2)
false
```
The second call of [`is_valid`](@ref) returns `false` because the model does 
not have parameter with index 2 yet.

We can also access different information about the domain that the infinite 
parameter is in. This is given by [`infinite_domain`](@ref), which takes a 
[`GeneralVariableRef`] as argument. For example, we have
```jldoctest param_queries
julia> infinite_domain(x)
[0, 1]
```
[`infinite_domain`](@ref) might be more useful if the infinite parameter is in a 
[`UniDistributionDomain`](@ref) or [`MultiDistributionDomain`](@ref), by which users 
can access information about the underlying distribution. On the other hand, 
if we already know that the parameter is in an interval domain, we can use 
[`JuMP.has_lower_bound`](@ref), [`JuMP.lower_bound`](@ref), 
[`JuMP.has_upper_bound`](@ref), [`JuMP.upper_bound`](@ref) to retrieve information 
about the interval domain in a more specific way:
```jldoctest param_queries
julia> has_lower_bound(x)
true

julia> lower_bound(x)
0.0

julia> has_upper_bound(x)
true

julia> upper_bound(x)
1.0
```

A quick way for users to obtain a [`GeneralVariableRef`](@ref) for a parameter with 
a known name would be through [`parameter_by_name`](@ref) function. This 
function takes an [`InfiniteModel`](@ref) and the parameter name in string, 
and returns a [`GeneralVariableRef`](@ref) for that parameter. For example, 
```jldoctest param_queries
julia> pref = parameter_by_name(model, "x")
x
```
If there is no parameter associated with that name, the function would return 
nothing. Otherwise, if multiple parameters share the same name, the function 
would throw an error.

Now we introduce two additional functions that we can use to access parameter 
information for an [`InfiniteModel`](@ref). The function 
[`num_parameters`](@ref) returns the number of infinite parameters associated 
with a model, while [`all_parameters`](@ref) returns the list of all infinite 
parameter references in the model. For a quick example: 
```jldoctest param_queries
julia> @infinite_parameter(model, y[1:2] in [0, 5])
2-element Vector{GeneralVariableRef}:
 y[1]
 y[2]

julia> num_parameters(model)
3

julia> all_parameters(model)
3-element Vector{GeneralVariableRef}:
 x   
 y[1]
 y[2]
```

## Parameter Modification
In this section we introduce a few shortcuts for users to modify defined 
infinite parameters.

First, once an infinite parameter is defined, we can change its name by calling 
the [`JuMP.set_name`] function, which takes the [`GeneralVariableRef`] of 
the parameter that needs a name change and the name string as arguments. For 
example, to change the parameter `x` to `t` we can do:
```jldoctest param_queries
julia> JuMP.set_name(x, "t")

julia> all_parameters(model)
3-element Vector{GeneralVariableRef}:
 t   
 y[1]
 y[2]
```
Similarly, we can also change the infinite domain that the parameter is in 
using the [`set_infinite_domain`](@ref) function as follows:
```jldoctest param_queries
julia> t = parameter_by_name(model, "t")
t

julia> set_infinite_domain(t, IntervalDomain(0, 5))

julia> infinite_domain(t)
[0, 5]
```

For parameters in an [`IntervalDomain`](@ref), we extend 
[`JuMP.set_lower_bound`](@ref) and [`JuMP.set_upper_bound`](@ref) functions 
for users to modify the lower bounds and upper bounds. For example,
```jldoctest param_queries
julia> JuMP.set_lower_bound(t, 1)

julia> JuMP.set_upper_bound(t, 4)

julia> infinite_domain(t)
[1, 4]
```
We do not support setting lower bounds and upper bounds for random parameters 
in a [`UniDistributionDomain`](@ref) and will throw an error if users attempt to do 
so. If users want to set lower bound and upper bound for a random infinite 
parameter, consider using `Distributions.Truncated`, which creates a truncated 
distribution from a univariate distribution.

## [Generative Supports](@id gen_supp_docs)
Generative supports denote supports that are generated based on existing supports  
(treated as finite elements). These are important for enabling certain measure  
and derivative evaluation schemes. Examples of such supports include internal  
collocation nodes and quadrature supports generated for quadrature methods that 
decompose the infinite domain such that existing supports are incorporated. Users 
shouldn't modify these directly, but extension writers will need to utilize the 
generative support API when developing measures and/or derivative evaluation 
methods that need to generate supports based on existing ones (e.g., adding 
a new orthogonal collocation method). More information about extension writing 
for either case is given on the [Extensions](@ref) page. For enhanced context, we 
outline the general API below.

Information about producing generative supports are stored via concrete subtypes 
of [`AbstractGenerativeInfo`](@ref). Each `IndependentParameter` stores one of 
these objects (the default being [`NoGenerativeSupports`](@ref)). Hence, a 
particular independent parameter can only be associated with 1 generative support 
scheme. We currently provide 1 concrete generative subtype of 
`AbstractGenerativeInfo` which is [`UniformGenerativeInfo`](@ref). 
`UniformGenerativeInfo` stores the necessary information to make generative 
supports that are uniformly applied to each finite element formed by the existing 
supports. For example, let's say we want to use a generative support scheme that 
adds 1 generative support exactly in the middle of each finite element with a 
unique support label to we'll call `MyGenLabel`:
```jldoctest; setup = :(using InfiniteOpt)
julia> struct MyGenLabel <: InfiniteOpt.InternalLabel end;

julia> UniformGenerativeInfo([0.5], MyGenLabel)
UniformGenerativeInfo([0.5], MyGenLabel)
```
Users can make other generative support schemes as described on the [Extensions](@ref) 
page. 

These `AbstractGenerativeInfo` objects are added to parameters as needed via the 
addition of measures and/or derivative methods that require generative supports. 
We can always check what generative information is currently associated with a 
particular parameter via [`generative_support_info`](@ref generative_support_info(::IndependentParameterRef)). 
The generation of these supports is handled automatically at the appropriate 
times via [`add_generative_supports`](@ref). We can always check if generative 
supports have been created for a particular parameter with 
[`has_generative_supports`](@ref).
