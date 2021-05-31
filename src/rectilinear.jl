# RectilinearGridFunction

export RectilinearGridFunction
abstract type RectilinearGridFunction{S,T,D} <: AbstractGridFunction{S,T,D} end

export bbox
bbox(::RectilinearGridFunction) = error("undefined")

export test_RectilinearGridFunction
function test_RectilinearGridFunction(gf::RectilinearGridFunction)
    test_AbstractGridFunction(gf)

    @testset "$(typeof(gf)) <: RectilinearGridFunction" begin
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
