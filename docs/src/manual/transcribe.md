# [Model Transcription](@id transcription_manual)
A technical manual for `InfiniteOpt.TranscriptionOpt`. See the respective 
[guide](@ref transcription_docs) for more information.

## Definition
```@docs
InfiniteOpt.TranscriptionOpt.TranscriptionModel
InfiniteOpt.TranscriptionOpt.TranscriptionData
InfiniteOpt.TranscriptionOpt.set_parameter_supports
InfiniteOpt.TranscriptionOpt.transcribe_finite_variables!
InfiniteOpt.TranscriptionOpt.transcribe_infinite_variables!
InfiniteOpt.TranscriptionOpt.transcribe_derivative_variables!
InfiniteOpt.TranscriptionOpt.transcribe_semi_infinite_variables!
InfiniteOpt.TranscriptionOpt.transcribe_point_variables!
InfiniteOpt.TranscriptionOpt.transcription_expression
InfiniteOpt.TranscriptionOpt.transcribe_measures!
InfiniteOpt.TranscriptionOpt.transcribe_objective!
InfiniteOpt.TranscriptionOpt.transcribe_constraints!
InfiniteOpt.TranscriptionOpt.transcribe_derivative_evaluations!
InfiniteOpt.TranscriptionOpt.transcribe_variable_collocation_restictions!
InfiniteOpt.TranscriptionOpt.build_transcription_model!
InfiniteOpt.add_point_variable(::JuMP.Model,::InfiniteOpt.GeneralVariableRef,::Vector{Float64},::Val{:TransData})
InfiniteOpt.add_semi_infinite_variable(::JuMP.Model,::InfiniteOpt.SemiInfiniteVariable,::Val{:TransData})
InfiniteOpt.build_optimizer_model!(::InfiniteOpt.InfiniteModel,::Val{:TransData})
```

## Queries
```@docs
InfiniteOpt.TranscriptionOpt.is_transcription_model
InfiniteOpt.TranscriptionOpt.transcription_data
InfiniteOpt.TranscriptionOpt.has_internal_supports
InfiniteOpt.TranscriptionOpt.transcription_model
InfiniteOpt.TranscriptionOpt.transcription_variable(::JuMP.Model,::InfiniteOpt.GeneralVariableRef)
InfiniteOpt.optimizer_model_variable(::InfiniteOpt.GeneralVariableRef,::Val{:TransData})
InfiniteOpt.variable_supports(::JuMP.Model,::Union{InfiniteOpt.InfiniteVariableRef, InfiniteOpt.SemiInfiniteVariableRef},::Val{:TransData})
InfiniteOpt.TranscriptionOpt.lookup_by_support(::JuMP.Model,::InfiniteOpt.GeneralVariableRef,::Vector)
InfiniteOpt.internal_semi_infinite_variable(::InfiniteOpt.SemiInfiniteVariableRef,::Val{:TransData})
InfiniteOpt.TranscriptionOpt.transcription_expression(::JuMP.Model,::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr})
InfiniteOpt.optimizer_model_expression(::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr},::Val{:TransData})
InfiniteOpt.expression_supports(::JuMP.Model,::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr}, ::Val{:TransData})
InfiniteOpt.TranscriptionOpt.transcription_constraint(::JuMP.Model,::InfiniteOpt.InfOptConstraintRef)
InfiniteOpt.optimizer_model_constraint(::InfiniteOpt.InfOptConstraintRef,::Val{:TransData})
InfiniteOpt.constraint_supports(::JuMP.Model,::InfiniteOpt.InfOptConstraintRef,::Val{:TransData})
InfiniteOpt.TranscriptionOpt.parameter_supports(::JuMP.Model)
```

## Utilities
```@docs
InfiniteOpt.TranscriptionOpt.support_index_iterator
InfiniteOpt.TranscriptionOpt.index_to_support
InfiniteOpt.TranscriptionOpt.index_to_labels
InfiniteOpt.TranscriptionOpt.make_ndarray
```