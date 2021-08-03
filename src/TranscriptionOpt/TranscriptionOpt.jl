module TranscriptionOpt

using ..InfiniteOpt

# Load in the module files
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
