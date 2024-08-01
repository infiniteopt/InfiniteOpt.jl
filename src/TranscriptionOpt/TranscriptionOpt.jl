module TranscriptionOpt

using ..InfiniteOpt

# Load in the module files
include("model.jl")
include("measures.jl")
include("transcribe.jl")

# Export transcription datatypes
export TranscriptionBackend

end # end module
