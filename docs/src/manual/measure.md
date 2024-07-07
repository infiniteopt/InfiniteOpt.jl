# [Measure Operators](@id measure_manual)
A technical manual for measures in `InfiniteOpt`. See the respective 
[guide](@ref measure_docs) for more information.

## Measure Data
```@docs
AbstractMeasureData
DiscreteMeasureData(::GeneralVariableRef, ::Vector{<:Real}, ::Vector{<:Real})
DiscreteMeasureData(::AbstractArray{GeneralVariableRef}, ::Vector{<:Real}, ::Vector{<:AbstractArray{<:Real}})
DiscreteMeasureData
FunctionalDiscreteMeasureData(::GeneralVariableRef,::Function,::Int,::Type{<:AbstractSupportLabel})
FunctionalDiscreteMeasureData(::AbstractArray{GeneralVariableRef},::Function,::Int,::Type{<:AbstractSupportLabel})
FunctionalDiscreteMeasureData
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
default_weight
```

## Definition
### General
```@docs
measure
@measure
build_measure
Measure
add_measure
InfiniteOpt.add_supports_to_parameters(::AbstractMeasureData)
MeasureIndex
MeasureData
MeasureRef
```

### Integrals
```@docs
InfiniteOpt.MeasureToolbox.integral(::JuMP.AbstractJuMPScalar,::InfiniteOpt.GeneralVariableRef,::Real,::Real)
InfiniteOpt.MeasureToolbox.∫(::JuMP.AbstractJuMPScalar,::InfiniteOpt.GeneralVariableRef,::Real,::Real)
InfiniteOpt.MeasureToolbox.uni_integral_defaults
InfiniteOpt.MeasureToolbox.set_uni_integral_defaults
InfiniteOpt.MeasureToolbox.clear_uni_integral_defaults
InfiniteOpt.MeasureToolbox.integral(::JuMP.AbstractJuMPScalar,::AbstractArray{InfiniteOpt.GeneralVariableRef},::Union{Real, AbstractArray{<:Real}},::Union{Real, AbstractArray{<:Real}})
InfiniteOpt.MeasureToolbox.∫(::JuMP.AbstractJuMPScalar,::AbstractArray{InfiniteOpt.GeneralVariableRef},::Union{Real, AbstractArray{<:Real}},::Union{Real, AbstractArray{<:Real}})
InfiniteOpt.MeasureToolbox.multi_integral_defaults
InfiniteOpt.MeasureToolbox.set_multi_integral_defaults
InfiniteOpt.MeasureToolbox.clear_multi_integral_defaults
InfiniteOpt.MeasureToolbox.@integral
InfiniteOpt.MeasureToolbox.@∫
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
InfiniteOpt.MeasureToolbox.generate_integral_data
```

### Expectations
```@docs
InfiniteOpt.MeasureToolbox.expect
InfiniteOpt.MeasureToolbox.@expect
InfiniteOpt.MeasureToolbox.𝔼
InfiniteOpt.MeasureToolbox.@𝔼
InfiniteOpt.MeasureToolbox.generate_expect_data
```

### Support Sum
```@docs
InfiniteOpt.MeasureToolbox.@support_sum
InfiniteOpt.MeasureToolbox.support_sum
```

## Queries
```@docs
JuMP.name(::MeasureRef)
num_measures
all_measures
measure_function
measure_data
is_analytic
parameter_refs(::MeasureRef)
is_used(::MeasureRef)
used_by_derivative(::MeasureRef)
used_by_constraint(::MeasureRef)
used_by_measure(::MeasureRef)
used_by_objective(::MeasureRef)
parameter_group_int_indices(::MeasureRef)
```

## Modification
```@docs
JuMP.set_name(::MeasureRef, ::String)
JuMP.delete(::InfiniteModel, ::MeasureRef)
```

## Expansion
```@docs
expand
expand_all_measures!
InfiniteOpt.expand_measure
InfiniteOpt.analytic_expansion
InfiniteOpt.expand_measures
make_point_variable_ref
make_semi_infinite_variable_ref
add_point_variable(::AbstractTransformationBackend, ::Any, ::Any)
add_semi_infinite_variable(::AbstractTransformationBackend, ::Any)
internal_semi_infinite_variable
```
