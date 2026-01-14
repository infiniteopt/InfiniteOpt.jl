# [Extension Packages](@id ext_docs)
Here, we document extension packages that extend the base capabilities of InfiniteOpt.

## [InfiniteInterpolations](@id interpolate)
This extension uses [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)
to automatically convert discretized solution values into continuous functions. This
extension is enabled by importing both InfiniteOpt and Interpolations:
```jldoctest interpolate
julia> using InfiniteOpt, Interpolations
```
Now let's solve a model:
```jldoctest interpolate
julia> using HiGHS; model = InfiniteModel(HiGHS.Optimizer);

julia> @infinite_parameter(model, t ∈ [0, 1], num_supports = 5);

julia> @variable(model, y >= 0, Infinite(t));

julia> @objective(model, Min, ∫(y, t));

julia> @constraint(model, y >= 2t);

julia> set_silent(model); optimize!(model)

julia> discrete_y = value(y)
5-element Vector{Float64}:
 0.0
 0.5
 1.0
 1.5
 2.0
```
where we get vector of five values for `y` since 5 supports are used to discretize ``y(t)``. To obtain a continuous representation, we can use interpolation by specifying one of three supported kinds:
- `constant_interpolation` or `Constant()`
- `linear_interpolation` or `Linear()`
- `cubic_spline_interpolation` or `Cubic()`

For our example, let's choose a linear interpolation:
```jldoctest interpolate
julia> continuous_y = value(y, linear_interpolation);

julia> continuous_y(0.1)
0.2

julia> continuous_y(0.25)
0.5
```
Note we could have equivalently used `value(y, Linear())`.

!!! warning
    There is a type piracy conflict between JuMP and OffsetArrays 
    (a dependency of Interpolations.jl). As a result, type piracy issues may arise 
    when Interpolations is loaded in (typically when `JuMP.DenseAxisArray`s are used). 
    Hence, we recommend using `Array` containers when using this extension.

## [InfiniteMathOptAI](@id mathoptai_guide)
This extension allows us to import machine learning models into `InfiniteModels`
via [MathOptAI](https://lanl-ansi.github.io/MathOptAI.jl/stable/). This is enabled by
importing InfiniteOpt and MathOptAI:
```jldoctest mathoptai
julia> using InfiniteOpt, MathOptAI
```
Now we can incorporate any predictor (i.e., machine learning model) supported by MathOptAI
which includes neural networks, Gaussian processes, decision trees, and generalized linear 
models. For instance, consider the neural ODE with the right-hand side:
```jldoctest mathoptai
julia> using Flux

julia> NN = Flux.Chain(Flux.Dense(2 => 3, Flux.relu), Flux.Dense(3 => 1));
```
Let's create an `InfiniteModel` those poses a neural operator constraint with `NN`. This
is accomplished by using [`MathOptAI.add_predictor`](https://lanl-ansi.github.io/MathOptAI.jl/stable/api/#add_predictor):
```jldoctest mathoptai
julia> model = InfiniteModel();

julia> @infinite_parameter(model, t ∈ [0, 1]);

julia> @variable(model, y >= 0, Infinite(t));

julia> f, formulation = add_predictor(model, NN, [y, t]);

julia> @constraint(model, ∂(y, t) == only(f))
d/dt[y(t)] - moai_Affine[1](t) = 0, ∀ t ∈ [0, 1]
```
Here, `f` is the vector of output variables (only 1 in this case) and `formulation`
is an object that stores all the variables and constraints created to embed `NN` in
`model`. To learn more about MathOptAI's syntax, please visit 
[MathOptAI's documentation](https://lanl-ansi.github.io/MathOptAI.jl/stable/).