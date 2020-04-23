# [Infinite Parameters] (@id inf_par_page)
A guide and manual to the definition and use of infinite parameters in
`InfiniteOpt`. The Datatypes and Methods sections at the end comprise the manual,
and the above sections comprise the guide.  

## Overview
Infinite parameters are what parameterize the infinite decision spaces for
infinite dimensional mathematical programs. In dynamic optimization this
corresponds to time and in stochastic optimization this to uncertain parameters
that follow a certain underlying statistical distribution. `InfiniteOpt`
considers principally two kinds of infinite parameters, ones defined over
continuous intervals and ones characterized by a distribution (others can be
added by defining a user-defined type). These can be used to parameterize
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
where the time parameter is stored in `model`. It can now be used with infinite
variables, measures, and constraints as described in the appropriate
sections.

When the model is optimized, `t` will be transcribed (discretized) over its domain
following its support points. If none are specified by the user than 50
support points are generated that are equidistant over the interval. However,
users may wish to employ their own support scheme. This can be done by using the
`num_supports` or `supports` keyword arguments. For example, if we desire to
have only 10 equidistant supports then we could have instead defined `t`:
```jldoctest; setup = :(using JuMP, InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10], num_supports = 10)
t
```
More complex support schemes can be specified via `supports` such as:
```jldoctest; setup = :(using JuMP, InfiniteOpt; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10], supports = [0; 2; 7; 10])
t
```
Where we specified `t` to use 4 supports: 0, 2, 7, and 10.

We can also add supports after `t` has been initialized. This can be accomplished
with [`add_supports`](@ref). For example, consider the initial case where `t` has
no supports and we now wish to add 4 supports:
```jldoctest basic
julia> add_supports(t, [0., 2.5, 7.5, 10.])

julia> supports(t)
4-element Array{Number,1}:
  0.0
  2.5
  7.5
 10.0
```
Here only 4 supports are specified for the sake of example. Alternatively, we
could have initialized the parameter and added supports in just one step using
the `supports` keyword argument:
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, t in [0, 10], supports = [0., 2.5, 7.5, 10.])
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
to each parameter. Supports can also be specified for each parameter as shown
above. The key difference is that by default 50 supports are generated randomly
following the distribution. Similarly, the `num_supports` keyword is used to
generate random supports. Now we have infinite parameters `t` and `xi` that are
ready to be used in defining infinite variables and constraints. We also mention
here that the [`@infinite_parameter`](@ref) macro is designed to closely emulate
[`JuMP.@variable`](@ref) and thus handles arrays and keyword arguments in the
same way. This is described in more detail below.

## Parameter Definition
Defining/initializing an infinite parameter principally involves the following
steps:
1. Define an [`AbstractInfiniteSet`](@ref)
2. Define support points within the set to later discretize the parameter
3. Specify if parameter is independent (only for multi-dimensional parameter groups)
4. Construct an [`InfOptParameter`](@ref) to store this information
5. Add the `InfOptParameter` object to an `InfiniteModel` and assign a name
6. Create a [`ParameterRef`](@ref) that points to the parameter object

### Manual Definition
Infinite set definition is described above in the
[Infinite Sets](@ref infinite_sets_normal) section. The
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
[0, 10]

julia> t_param = build_parameter(error, set, supports = [0, 2, 5, 7, 10])
InfOptParameter{IntervalSet}([0, 10], [0, 2, 5, 7, 10], false)
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
Normal{Float64}(μ=0.0, σ=1.0)

julia> x_param = build_parameter(error, set, supports = [-0.5, 0.5])
InfOptParameter{DistributionSet{Normal{Float64}}}(Normal{Float64}(μ=0.0, σ=1.0), [-0.5, 0.5], false)
```
Again, we use [`add_parameter`](@ref) to add `x_param` to the `InfiniteModel` and
assign it the name `x`:
```jldoctest rand_define
julia> x_ref = add_parameter(model, x_param, "x")
x
```
Note that `add_parameter` does not register the name of the parameters into the
model that it adds to. As shown in [Macro Definition](@ref param_macro), the macro definition
does not allow for multiple parameters sharing the same name and will throw an
error if it happens.

### [Macro Definition] (@id param_macro)
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
when we introduce how to define anonymous parameters. See the part for
anonymous parameter definition for more details.

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

For multi-dimensional parameters, the macro calls for [`build_parameter`](@ref)
and [`add_parameter`](@ref) and creates looped codes that construct separate
parameters and references for each dimension. If an array of supports is
provided, the macro will assign that array of supports to all dimensions.
Otherwise, the indexed syntax can be used to feed in different array of
supports to each dimension, similar to [`JuMP.@variable`](@ref). For example:
```jldoctest; setup = :(using InfiniteOpt; model= InfiniteModel())
julia> points = [0.2 0.8; 0.3 0.7]
2×2 Array{Float64,2}:
 0.2  0.8
 0.3  0.7

julia> @infinite_parameter(model, a[i = 1:2] in [0, 1], supports = points[i, :])
2-element Array{ParameterRef,1}:
 a[1]
 a[2]

julia> supports(a[1])
2-element Array{Float64,1}:
 0.2
 0.8

julia> supports(a[2])
2-element Array{Float64,1}:
 0.3
 0.7
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
`SparseAxisArray` from the `JuMP` package for the space
parameter `x`:
```jldoctest; setup = :(using InfiniteOpt, JuMP; model= InfiniteModel())
julia> @infinite_parameter(model, x[1:3] in [0, 1], container = SparseAxisArray)
JuMP.Containers.SparseAxisArray{ParameterRef,1,Tuple{Int64}} with 3 entries:
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
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> x = @infinite_parameter(model, [1:3], lower_bound = 0, upper_bound = 1)
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
julia> @infinite_parameter(model, [1:3], lower_bound = 0, upper_bound = 1, base_name = "x")
3-element Array{ParameterRef,1}:
 x[1]
 x[2]
 x[3]

julia> @infinite_parameter(model, [1:3], lower_bound = -1, upper_bound = 0, base_name = "x")
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
This section is for people who wish to know more about how the macro
[`@infinite_parameter`](@ref) works in the backend. Users who only want to learn
about the setting up the model can skip over this part.

In general, the macro [`@infinite_parameter`](@ref) follows the same steps as
the manual definition. First, it parses the arguments and identifies any
recognizable keyword arguments. Specifically, the first argument must be the
model, and the second argument, if exists, must be an expression that declares
the parameter or simply specify the dimension of the parameter if users choose
to define it anonymously. If the information in the keyword arguments is not
sufficient to define the set the parameter is in, the users also need to specify
the sets in the second argument using expressions like `a <= x <= b` or
`x in set`.

The keyword arguments give users flexibility in how to define their parameters.
As mentioned above, the users can choose to specify the set either in the
second argument (nonanonymous parameter definition only), or in the keyword
arguments. However, the users cannot do both at the same time. The macro will
check this behavior and throw an error if this happens. For example,
```jldoctest; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model,  y in [0, 1], lower_bound = 0, upper_bound = 1)
ERROR: LoadError

julia> @infinite_parameter(model,  y in [0, 1], set = IntervalSet(0, 1))
ERROR: LoadError
```

Once the check on arguments and keyword arguments is done, the macro will create
the [`AbstractInfiniteSet`](@ref) based on given information, and create the
infinite parameter accordingly. If the users create a multi-dimensional
parameter, the macro will create looped code to define individual infinite
parameter for each dimension. The looped code will also incorporate different
supports for different dimensions.

In the end, if the created parameter is not anonymous, the macro will register
the name to the model. In this way, we prevent parameters created by
[`@infinite_parameter`](@ref) non-anonymously to share the same name.


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
parameter has been assigned with supports. Users can also delete all the
supports of a parameter with [`delete_supports`](@ref).

### Automatic Support Generation During Parameter Definition
For the examples in the [Parameter Definition](@ref), we have seen how to
manually add supports to an infinite parameter. For a quick automatic
generation of support points, though, users do not have to input the support
points. Instead, the number of support points generated is supplied.

For an infinite parameter subject to an [`IntervalSet`](@ref), uniformly spaced
supports including both ends are generated across the interval. For example,
defining a time parameter ``t \in [0, 10]`` with 4 supports using
[`build_parameter`](@ref) gives
```jldoctest; setup = :(using InfiniteOpt)
julia> set = IntervalSet(0, 10)
[0, 10]

julia> t_param = build_parameter(error, set, num_supports = 4, sig_fig = 3)
InfOptParameter{IntervalSet}([0, 10], [0.0, 3.33, 6.67, 10.0], false)
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
```jldoctest; setup = :(using InfiniteOpt, Distributions, Random; Random.seed!(0); model = InfiniteModel(); dist = Normal(0., 1.))
julia> @infinite_parameter(model, x in dist, num_supports = 4)
x

julia> supports(x)
4-element Array{Float64,1}:
  0.67911
  0.82841
 -0.35301
 -0.13485

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
argument a [`ParameterRef`](@ref), in which case it will generate
supports for the associated infinite parameter. Alternatively,
[`fill_in_supports!`](@ref) can also take an [`InfiniteModel`](@ref) as an
argument, in which case it will generate supports for all infinite parameters
of the [`InfiniteModel`](@ref) with no supports.

The [`fill_in_supports!`](@ref) method allows users to specify integer keyword
arguments `num_supports` and `sig_fig`. `num_supports` dictates the
number of supports to be generated, and `sig_fig` dictates the significant
figures of generated supports desired. The default values are 50 and 5,
respectively.

The ways by which supports are automatically generated are as follows. If the
parameter is in an [`IntervalSet`](@ref), then we generate an array of supports
that are uniformly distributed along the interval, including the two ends. For
example, consider a 3D position parameter `x` distributed in the unit cube
`[0, 1]`. We can generate supports for that point in the following way:
```jldoctest supp_gen_defined; setup = :(using InfiniteOpt, JuMP, Distributions, Random; Random.seed!(0); model = InfiniteModel())
julia> @infinite_parameter(model, x[1:3] in [0, 1], independent = true);

julia> fill_in_supports!.(x, num_supports = 3);

julia> supports.(x)
3-element Array{Array{Number,1},1}:
 [0.0, 0.5, 1.0]
 [0.0, 0.5, 1.0]
 [0.0, 0.5, 1.0]
```
Note that the dot syntax because [`fill_in_supports!`](@ref) takes single
[`ParameterRef`](@ref) as argument. In each dimension, three equally spaced
supports (`[0.0, 0.5, 1.0]`) are generated. Since the `independent` keyword is
set as `true`, the transcription stage will create a three-dimensional grid for
all variables parameterized by `x`, with each point separated by 0.5 units in
each dimension.

If the parameter is in a [`DistributionSet`](@ref), [`fill_in_supports!`](@ref)
samples `num_supports` supports from the distribution. Recall that support
generation is not allowed for parameters under multivariate distribution during
parameter definition. However, if the parameter is defined first without
supports, [`fill_in_supports!`](@ref) allows for supports generation. For
example, for a 2D random variable `xi` under a multivariate Gaussian
distribution, we can generate supports for it in the following way:
```jldoctest supp_gen_defined
julia> dist = MvNormal([0., 0.], [1. 0.; 0. 2.])
FullNormal(
dim: 2
μ: [0.0, 0.0]
Σ: [1.0 0.0; 0.0 2.0]
)


julia> @infinite_parameter(model, xi[1:2] in dist);

julia> fill_in_supports!(xi[1], num_supports = 3)

julia> supports.(xi)
2-element Array{Array{Number,1},1}:
 [0.67911, -0.35301, 0.58662]
 [1.1716, -0.19071, 0.4205]
```
Note that for parameters under multivariate distribution, the users do not have
to call [`fill_in_supports!`](@ref) for all [`ParameterRef`](@ref) associated
with that parameter. [`fill_in_supports!`](@ref), once called upon one
[`ParameterRef`](@ref), would search for the other [`ParameterRef`](@ref)
and fill in supports for them if they, together with the original
[`ParameterRef`](@ref), make up for the same multi-dimensional parameter under
a multivariate distribution.

Also, [`fill_in_supports!`](@ref) only fill in supports for parameters with no
associated supports. To modify the supports of parameters already associated
with some supports, refer to [Supports](@ref) for how to do that.

## Parameter Queries
In addition to the modeling framework, this package provides many functions for
users to access information about the model. This section will go over basic
functions for accessing parameter information.

Once a (possibly large-scale) `InfiniteModel` is built, the users might want to
check if an infinite parameter is actually used in any way. This could be
checked by [`is_used`](@ref) function as follows:
```jldoctest param_queries; setup = :(using InfiniteOpt, JuMP; model = InfiniteModel())
julia> @infinite_parameter(model, x in [0, 1])
x

julia> is_used(x)
false

```
This function checks if the parameter is used by any constraint, measure, or
variable. In a similar way, functions [`used_by_constraint`](@ref),
[`used_by_measure`](@ref) and [`used_by_variable`](@ref) can be applied to
find out any dependency of specific types on the infinite parameter.

In addition, sometimes we need to check if a certain [`ParameterRef`](@ref) is
valid with an `InfiniteModel` model, meaning that the parameter reference
actually refers to some parameter associated with the model. We extend the
[`JuMP.is_valid`](@ref) function from JuMP for that purpose. To see how to use
this, for example,
```jldoctest param_queries
julia> pref1 = ParameterRef(model, 1);

julia> pref2 = ParameterRef(model, 2);

julia> is_valid(model, pref1)
true

julia> is_valid(model, pref2)
false
```
The second call of [`is_valid`](@ref) returns `false` because the model does
not have parameter with index 2 yet.

We can also access different information about the set that the infinite
parameter is in. This is given by [`infinite_set`](@ref), which takes a
[`ParameterRef`] as argument. For example, we have
```jldoctest param_queries
julia> infinite_set(x)
[0, 1]
```
[`infinite_set`](@ref) might be more useful if the infinite parameter is in a
[`DistributionSet`](@ref), by which users can access information about the
underlying distribution. On the other hand, if we already know that the
parameter is in an interval set, we can use [`JuMP.has_lower_bound`](@ref),
[`JuMP.lower_bound`](@ref), [`JuMP.has_upper_bound`](@ref),
[`JuMP.upper_bound`](@ref) to retrieve information about the interval set in a
more specific way:
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

For multi-dimensional parameters, the definition step would create the
parameter for each dimension an add it into the model separately. However,
users can find out if different parameters belong to the same group using
[`group_id`](@ref) function. The group ID is simply an index assigned at the
parameter definition. See the following example for how to use it:
```jldoctest param_queries
julia> @infinite_parameter(model, y[1:2] in [0, 5])
2-element Array{ParameterRef,1}:
 y[1]
 y[2]

julia> group_id(x)
1

julia> group_id(y)
2

julia> group_id([x, y[1]])
ERROR: Array contains parameters from multiple groups.
```
In the example above, the function returns the group ID if it receives one
[`ParameterRef`](@ref). If the argument is a vector of [`ParameterRef`](@ref),
however, it will throw an error if the parameters do not belong to the same
group. Otherwise, it will return the group ID for all parameters in the vector.
Again, this function could be helpful if the users want to verify if a set of
parameters are correlated.

Another important piece of information about multi-dimensional parameters is
whether different dimensions are independent. This can be checked by
[`is_independent`](@ref):
```jldoctest param_queries
julia> is_independent(y[1])
false
```

A quick way for users to obtain a [`ParameterRef`](@ref) for a parameter with
a known name would be through [`parameter_by_name`](@ref) function. This
function takes an [`InfiniteModel`](@ref) and the parameter name in string,
and returns a [`ParameterRef`](@ref) for that parameter. For example,
```jldoctest param_queries
julia> pref = parameter_by_name(model, "x")
x
```
If there is no parameter associated with that name, the function would return
nothing. Otherwise, if multiple parameters share the same name, the function
would throw an error.

Now we introduce two additional functions that we can use to access parameter
information for an  [`InfiniteModel`](@ref). The function
[`num_parameters`](@ref) returns the number of infinite parameters associated
with a model, while [`all_parameters`](@ref) returns the list of all infinite
parameter references in the model. For a quick example:
```jldoctest param_queries
julia> num_parameters(model)
3

julia> all_parameters(model)
3-element Array{ParameterRef,1}:
 x   
 y[1]
 y[2]
```

## Parameter Modification
In this section we introduce a few shortcuts for users to modify defined
infinite parameters.

First, once an infinite parameter is defined, we can change its name by calling
the [`JuMP.set_name`] function, which takes the [`ParameterRef`] that needs a
name change and the name string as arguments. For example, to change the
parameter `x` to `t` we can do:
```jldoctest param_queries
julia> JuMP.set_name(x, "t")

julia> all_parameters(model)
3-element Array{ParameterRef,1}:
 t   
 y[1]
 y[2]
```
In a similar way, we can also change the infinite set that the parameter is in
using the [`set_infinite_set`](@ref) function as follows:
```jldoctest param_queries
julia> t = parameter_by_name(model, "t")
t

julia> set_infinite_set(t, IntervalSet(0, 5))

julia> infinite_set(t)
[0, 5]
```

For parameters in an [`IntervalSet`](@ref), we extend
[`JuMP.set_lower_bound`](@ref) and [`JuMP.set_upper_bound`](@ref) functions
for users to modify the lower bounds and upper bounds. For example,
```jldoctest param_queries
julia> JuMP.set_lower_bound(t, 1)

julia> JuMP.set_upper_bound(t, 4)

julia> infinite_set(t)
[1, 4]
```
We do not support setting lower bounds and upper bounds for random parameters
in a [`DistributionSet`](@ref) and will throw an error if users attempt to do
so. If users want to set lower bound and upper bound for a random infinite
parameter, consider using `Distributions.Truncated`, which creates a truncated
distribution from a univariate distribution.

We can also modify the `independent` status of a parameter using
[`set_independent`](@ref) and [`unset_independent`](@ref) functions. The
[`set_independent`](@ref) function sets `true` for the `independent` status,
while [`unset_independent`](@ref) function sets `false` for the `independent`
status.
```jldoctest param_queries
julia> set_independent(y[1])

julia> is_independent(y[1])
true

julia> is_independent(y[2])
false

julia> unset_independent(y[1])

julia> is_independent(y[1])
false
```
Note that these two functions takes one single [`ParameterRef`] as argument,
and it will not automatically find the other associated parameters and reset
their status altogether. This gives more flexibility on how we want to
discretize certain parameters. For example, if for a 3D position parameter `x`
we want to have uniform discretization along `x[1]`, but custom supports over
the 2D space `x[2]` and `x[3]`, we can set `x[1]` only as independent.


## Datatypes
```@index
Pages   = ["parameter.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
InfOptParameter
ScalarParameter
IndependentParameter
DependentParameters
ScalarParameterData
MultiParameterData
IndependentParameterIndex
DependentParametersIndex
DependentParameterIndex
IndependentParameterRef
DependentParameterRef
```

## Methods/Macros
```@index
Pages   = ["parameter.md"]
Modules = [InfiniteOpt, JuMP]
Order   = [:macro, :function]
```
```@docs
@infinite_parameter
@independent_parameter
@dependent_parameters
add_parameters
JuMP.name(::DependentParameterRef)
JuMP.set_name(::DependentParameterRef, ::String)
used_by_infinite_variable(::DependentParameterRef)
used_by_measure(::DependentParameterRef)
used_by_constraint(::DependentParameterRef)
is_used(::DependentParameterRef)
infinite_set(::DependentParameterRef)
infinite_set(::AbstractArray{<:DependentParameterRef})
set_infinite_set(::DependentParameterRef,::InfiniteScalarSet)
set_infinite_set(::AbstractArray{<:DependentParameterRef},::InfiniteArraySet)
JuMP.has_lower_bound(::DependentParameterRef)
JuMP.lower_bound(::DependentParameterRef)
JuMP.set_lower_bound(::DependentParameterRef,::Real)
JuMP.has_upper_bound(::DependentParameterRef)
JuMP.upper_bound(::DependentParameterRef)
JuMP.set_upper_bound(::DependentParameterRef,::Real)
num_supports(::DependentParameterRef)
num_supports(::AbstractArray{<:DependentParameterRef})
has_supports(::DependentParameterRef)
has_supports(::AbstractArray{<:DependentParameterRef})
supports(::DependentParameterRef)
supports(::AbstractArray{<:DependentParameterRef})
set_supports(::AbstractArray{<:DependentParameterRef},::AbstractArray{<:Vector{<:Real}})
add_supports(::AbstractArray{<:IndependentParameterRef},::AbstractArray{<:Vector{<:Real}})
delete_supports(::AbstractArray{<:DependentParameterRef})
generate_and_add_supports!(::AbstractArray{<:DependentParameterRef},::InfiniteArraySet)
fill_in_supports!(::AbstractArray{<:DependentParameterRef})
fill_in_supports!(::InfiniteModel)
```
