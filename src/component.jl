# AbstractGridComponent

export AbstractGridComponent
"""
    abstract type AbstractGridComponent{S,T,D} <: AbstractGridMesh{S,T,D} end
"""
abstract type AbstractGridComponent{S,T,D} <: AbstractGridMesh{S,T,D} end

export box
box(::AbstractGridComponent) = error("undefined")

export test_AbstractGridComponent
function test_AbstractGridComponent(gf::AbstractGridComponent)
    test_AbstractGridMesh(gf)

    @testset "$(typeof(gf)) <: AbstractGridComponent" begin
        D = ndims(gf)
        gfbox = box(gf)
        @test gfbox isa BBox{D}

        T = eltype(gf)
        for n in 1:10
            ind = map((a, b) -> rand(a:b), min(gfbox), max(gfbox))
            val = gf[CartesianIndex(Tuple(ind))]
            @test val isa T
        end
    end
end
