# Model Transcription
A guide and manual for transcribing infinite models using `InfiniteOpt`.

## Overview


## Basic Usage


## TranscriptionOpt


## Datatypes
```@index
Pages   = ["transcribe.md"]
Modules = [InfiniteOpt]
Order   = [:type]
```
```@docs
TranscriptionData
```

## Methods
```@index
Pages   = ["transcribe.md"]
Modules = [InfiniteOpt, InfiniteOpt.TranscriptionOpt]
Order   = [:function]
```
```@docs
transcription_model
build_optimizer_model!(::InfiniteModel,::Val{:TransData})
TranscriptionModel()
TranscriptionModel(::InfiniteModel)
is_transcription_model
transcription_data
transcription_variable
InfiniteOpt.supports(::JuMP.Model, ::InfiniteOpt.InfiniteVariableRef)
InfiniteOpt.supports(::InfiniteOpt.InfiniteVariableRef)
transcription_constraint
InfiniteOpt.supports(::JuMP.Model, ::InfiniteOpt.InfiniteConstraintRef)
InfiniteOpt.supports(::InfiniteOpt.InfiniteConstraintRef)
InfiniteOpt.parameter_refs(::JuMP.Model, ::InfiniteOpt.InfiniteConstraintRef)
InfiniteOpt.parameter_refs(::InfiniteOpt.InfiniteConstraintRef)
```
