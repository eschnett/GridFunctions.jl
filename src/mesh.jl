# AbstractGridMesh

export AbstractGridMesh
"""
    abstract type AbstractGridMesh{S,T,D} <: AbstractGridFunction{S,T,D} end
"""
abstract type AbstractGridMesh{S,T,D} <: AbstractGridFunction{S,T,D} end

export bbox
bbox(::AbstractGridMesh) = error("undefined")

export test_AbstractGridMesh
function test_AbstractGridMesh(gf::AbstractGridMesh)
    test_AbstractGridFunction(gf)

    @testset "$(typeof(gf)) <: AbstractGridMesh" begin
        D = ndims(gf)
        gfbbox = bbox(gf)
        @test gfbbox isa BBox{D}

        gfdomain = domain(gf)

        for n in 1:10
            ind = map((a, b) -> rand(a:b), min(gfbbox), max(gfbbox))
            coord = coordinate(gf, ind)
            @test all(min(gfdomain) .≤ coord .≤ max(gfdomain))
            ind′ = index(gf, coord)
            @test ind′ ≈ ind
        end
    end
end
