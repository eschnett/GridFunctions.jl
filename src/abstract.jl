# AbstractGridFunction

export AbstractGridFunction
abstract type AbstractGridFunction{S,T,D} end

Base.ndims(::AbstractGridFunction{S,T,D} where {S,T}) where {D} = D

export indextype
indextype(::AbstractGridFunction{S}) where {S} = S

Base.eltype(::AbstractGridFunction{S,T}) where {S,T} = T

export domain
domain(::AbstractGridFunction) = error("undefined method")

export test_AbstractGridFunction
function test_AbstractGridFunction(gf::AbstractGridFunction)
    @testset "$(typeof(gf)) <: AbstractGridFunction" begin
        D = ndims(gf)
        @test D isa Int
        @test D ≥ 0

        S = indextype(gf)
        @test S isa Type

        T = eltype(gf)
        @test T isa Type

        dom = domain(gf)
        @test dom isa Domain
        @test ndims(dom) == D
        @test eltype(dom) == S

        len = length(gf)
        @test len isa Int
        @test len ≥ 0

        rng = MersenneTwister(0)

        for val in [1, 1.0, SVector(1.0f0, 1.0f0)]
            T = typeof(val)

            function test_props(gf)
                @test gf isa AbstractGridFunction
                @test indextype(gf) == S
                @test eltype(gf) == T
                @test ndims(gf) == D
            end

            e = map(_ -> val, gf)
            test_props(e)

            @test reduce(+, e; init=zero(val)) == len * val
            @test mapreduce(sum, +, e; init=sum(zero(val))) == length(val) * len

            c = copy(e)
            s = similar(e)
            test_props(c)
            test_props(s)
            fill!(s, val)
            @test s == c
            @test s ≈ e

            n = zero(e)
            if T <: AbstractFloat
                x = map(v -> T(v[1]), coordinates(e))
            else
                x = rand(rng, e)
            end
            y = rand(rng, e)
            z = rand(rng, e)
            a = rand(rng, eltype(T))
            b = rand(rng, eltype(T))
            if T <: Integer
                # Prevent integer overflow
                map!(v -> v % 1000, x, x)
                map!(v -> v % 1000, y, y)
                map!(v -> v % 1000, z, z)
                a %= 1000
                b %= 1000
            end

            test_props(n)
            test_props(x)
            test_props(y)
            test_props(z)

            @test map(identity, y) == y
            @test map(a -> a .+ 1, map(a -> 2a, y)) ≈
                  map((a -> a .+ 1) ∘ (a -> 2a), y)

            @test +x == x
            @test n + x == x
            @test x + n == x
            @test (x + y) + z ≈ x + (y + z)

            @test x + (-x) == n
            @test x - y == x + (-y)

            @test zero(a) * x == n
            @test one(a) * x == x

            @test a * n == n

            @test a * x == x * a

            @test (a * b) * x ≈ a * (b * x)
            @test (x * a) * b ≈ x * (a * b)

            @test a * (x + y) ≈ a * x + a * y
            @test (a + b) * x ≈ a * x + b * x

            if a ≠ 0
                @test inv(a) * (a * x) ≈ x
                @test a \ x ≈ inv(a) * x
                @test x / a ≈ x * inv(a)
            end
        end
    end
end
