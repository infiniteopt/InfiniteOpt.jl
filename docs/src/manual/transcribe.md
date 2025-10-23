# [TranscriptionOpt](@id transcription_manual)
A technical manual for `InfiniteOpt.TranscriptionOpt` (the default transformation 
backend). See the respective [guide](@ref transcription_docs) for more information.

## Definition
```@docs
InfiniteOpt.TranscriptionOpt.TranscriptionBackend
InfiniteOpt.TranscriptionOpt.TranscriptionData
InfiniteOpt.TranscriptionOpt.Transcription
InfiniteOpt.TranscriptionOpt.set_parameter_supports
InfiniteOpt.TranscriptionOpt.transcribe_finite_parameters!
InfiniteOpt.TranscriptionOpt.transcribe_finite_variables!
InfiniteOpt.TranscriptionOpt.transcribe_infinite_variables!
InfiniteOpt.TranscriptionOpt.transcribe_parameter_functions!
InfiniteOpt.TranscriptionOpt.transcribe_derivative_variables!
InfiniteOpt.TranscriptionOpt.transcribe_semi_infinite_variables!
InfiniteOpt.TranscriptionOpt.transcribe_point_variables!
InfiniteOpt.TranscriptionOpt.transcription_expression
InfiniteOpt.TranscriptionOpt.transcribe_measures!
InfiniteOpt.TranscriptionOpt.transcribe_objective!
InfiniteOpt.TranscriptionOpt.transcribe_constraints!
InfiniteOpt.TranscriptionOpt.transcribe_derivative_evaluations!
InfiniteOpt.TranscriptionOpt.transcribe_variable_collocation_restictions!
InfiniteOpt.TranscriptionOpt.build_transcription_backend!
InfiniteOpt.add_point_variable(::InfiniteOpt.TranscriptionOpt.TranscriptionBackend,::InfiniteOpt.GeneralVariableRef,::Vector{Float64})
InfiniteOpt.add_semi_infinite_variable(::InfiniteOpt.TranscriptionOpt.TranscriptionBackend,::InfiniteOpt.SemiInfiniteVariable)
InfiniteOpt.build_transformation_backend!(::InfiniteOpt.InfiniteModel,::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
```

## Queries
```@docs
InfiniteOpt.TranscriptionOpt.transcription_data
InfiniteOpt.TranscriptionOpt.has_internal_supports
InfiniteOpt.TranscriptionOpt.transcription_variable(::InfiniteOpt.GeneralVariableRef, ::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
InfiniteOpt.transformation_variable(::InfiniteOpt.GeneralVariableRef,::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
InfiniteOpt.variable_supports(::Union{InfiniteOpt.InfiniteVariableRef, InfiniteOpt.SemiInfiniteVariableRef},::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
InfiniteOpt.TranscriptionOpt.lookup_by_support(::InfiniteOpt.GeneralVariableRef,::InfiniteOpt.TranscriptionOpt.TranscriptionBackend,::Vector)
InfiniteOpt.internal_semi_infinite_variable(::InfiniteOpt.SemiInfiniteVariableRef,::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
InfiniteOpt.TranscriptionOpt.transcription_expression(::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr}, ::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
InfiniteOpt.transformation_expression(::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr},::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
InfiniteOpt.expression_supports(::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr}, ::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
InfiniteOpt.TranscriptionOpt.transcription_constraint(::InfiniteOpt.InfOptConstraintRef, ::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
InfiniteOpt.transformation_constraint(::InfiniteOpt.InfOptConstraintRef,::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
InfiniteOpt.constraint_supports(::InfiniteOpt.InfOptConstraintRef,::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
InfiniteOpt.TranscriptionOpt.parameter_supports(::InfiniteOpt.TranscriptionOpt.TranscriptionBackend)
```

## Updates
```@docs
InfiniteOpt.update_parameter_value(::InfiniteOpt.TranscriptionOpt.TranscriptionBackend, ::InfiniteOpt.FiniteParameterRef, ::Real)
InfiniteOpt.update_parameter_value(::InfiniteOpt.TranscriptionOpt.TranscriptionBackend, ::InfiniteOpt.ParameterFunctionRef, ::Function)
```

## Utilities
```@docs
InfiniteOpt.TranscriptionOpt.support_index_iterator
InfiniteOpt.TranscriptionOpt.index_to_support
```