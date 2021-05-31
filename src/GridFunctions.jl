module GridFunctions

using ADIOS2
using Random
using StaticArrays
using Test

################################################################################

include("defs.jl")

include("domain.jl")
include("bbox.jl")
include("bboxset.jl")
include("bboxarray.jl")

include("abstract.jl")
include("rectilinear.jl")

include("gridarray.jl")

end
