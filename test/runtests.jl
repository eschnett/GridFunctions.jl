using GridFunctions
using StaticArrays
using Test

@testset "GridArray D=$D" for D in 1:3
    dom = Domain{Float64,D}(fill(0, SVector{D,Float64}), fill(1, SVector{D,Float64}))
    bb = BBox{D}(fill(0, SVector{D,Int}), fill(10, SVector{D,Int}))
    cc = fill(false, SVector{D,Bool})

    ga = GridArray{Float64,Nothing,D}(dom, bb, cc, fill(nothing, Tuple(size(bb))))

    test_AbstractGridComponent(ga)
end
