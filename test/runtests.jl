using GridFunctions
using StaticArrays
using Test

svector(D, x) = SVector(ntuple(_ -> x, D))

@testset "GridArray D=$D" for D in 1:3
    dom = Domain{Float64,D}(svector(D, 0), svector(D, 1))
    bb = BBox{D}(svector(D, 0), svector(D, 10))
    cc = svector(D, false)

    ga = GridArray{Float64,Nothing,D}(dom, bb, cc,
                                      fill(nothing, Tuple(size(bb))))

    test_RectilinearGridFunction(ga)
end
