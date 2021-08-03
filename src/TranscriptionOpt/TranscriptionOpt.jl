module TranscriptionOpt

# Import the necessary packages.
import JuMP
const MOI = JuMP.MOI
const MOIU = MOI.Utilities
using ..InfiniteOpt

include("model.jl")
include("measures.jl")
include("transcribe.jl")
include("optimize.jl")

# Export transcription datatypes
export TranscriptionData, TranscriptionModel

# Export transcription methods
export is_transcription_model, transcription_data, transcription_variable,
transcription_constraint, transcription_model, transcription_expression

end # end module
