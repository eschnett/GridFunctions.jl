# GridArray

export GridArray
struct GridArray{S,T,D} <: RectilinearGridFunction{S,T,D}
    domain::Domain{S,D}
    bbox::BBox{D}
    cell_centered::SVector{D,Bool}
    data::AbstractArray{T,D}
end

function invariant(ga::GridArray)
    size(bbox) == size(data) || return false
    return true
end

domain(ga::GridArray) = ga.domain
bbox(ga::GridArray) = ga.bbox

Base.size(ga::GridArray) = size(bbox(ga))
Base.length(ga::GridArray) = length(bbox(ga))

function Base.similar(ga::GridArray{S,T,D},
                      element_type::Type{R}=T) where {S,T,D,R}
    return GridArray{S,R,D}(ga.domain, ga.bbox, ga.cell_centered,
                            similar(ga.data, R))
end

function Base.copy(ga::GridArray{S,T,D}) where {S,T,D}
    return GridArray{S,T,D}(ga.domain, ga.bbox, ga.cell_centered, copy(ga.data))
end

function Base.fill!(dest_ga::GridArray{S,R,D}, a) where {S,R,D}
    fill!(dest_ga.data, a)
    return dest_ga
end

function Base.map(f, ga::GridArray{S,T,D} where {T},
                  gas::(GridArray{S,T,D} where {T})...) where {S,D}
    resdata = map(f, ga.data, (ga.data for ga in gas)...)
    R = eltype(resdata)
    return GridArray{S,R,D}(ga.domain, ga.bbox, ga.cell_centered, resdata)
end

function Base.map!(f, dest_ga::GridArray{S,R,D}, ga::GridArray{S,T,D} where {T},
                   gas::(GridArray{S,T,D} where {T})...) where {S,R,D}
    map!(f, dest_ga.data, ga.data, (ga.data for ga in gas)...)
    return dest_ga
end

function Base.mapreduce(f, op, ga::GridArray{S,T,D} where {T},
                        gas::(GridArray{S,T,D} where {T})...; init) where {S,D}
    res = mapreduce(f, op, ga.data, (ga.data for ga in gas)...; init=init)
    return res
end

function Base.reduce(op, ga::GridArray; init)
    res = reduce(op, ga.data; init=init)
    return res
end

function Random.rand!(rng::AbstractRNG, dest_ga::GridArray)
    rand!(rng, dest_ga.data)
    return dest_ga
end
Random.rand(rng::AbstractRNG, ga::GridArray) = rand!(rng, similar(ga))

Random.rand!(dest_ga::GridArray) = rand!(Random.GLOBAL_RNG, dest_ga)
Random.rand(ga::GridArray) = rand(Random.GLOBAL_RNG, ga)

function Base.:(==)(ga1::GridArray{S,T1,D},
                    ga2::GridArray{S,T2,D}) where {S,T1,T2,D}
    return ga1.domain == ga2.domain &&
           ga1.bbox == ga2.bbox &&
           ga1.cell_centered == ga2.cell_centered &&
           ga1.data == ga2.data
end

function find_rtol(::Type{T1}, ::Type{T2}) where {T1,T2}
    T = typeof(zero(T1) + zero(T2))
    T = eltype(T)               # needed if T is an array type
    T <: Union{Integer,Rational} && return 0
    return sqrt(eps(T))
end

function Base.isapprox(ga1::GridArray{S,T1,D}, ga2::GridArray{S,T2,D};
                       atol::Real=0,
                       rtol::Real=atol > 0 ? 0 : find_rtol(T1, T2),
                       nans::Bool=false) where {S,T1,T2,D}
    return ga1.domain == ga2.domain &&
           ga1.bbox == ga2.bbox &&
           ga1.cell_centered == ga2.cell_centered &&
           isapprox(ga1.data, ga2.data; atol=atol, rtol=rtol, nans=nans)
end

export iota!
function iota!(dest_ga::GridArray)
    for i in CartesianIndices(Tuple(size(dest_ga.bbox)))
        dest_ga.data[i] = SVector(Tuple(i)) + min(dest_ga.bbox) .- 1
    end
    return dest_ga
end

export iota
iota(ga::GridArray{S,T,D}) where {S,T,D} = iota!(similar(ga, SVector{D,Int}))

export coordinate
function coordinate(ga::GridArray{S,X,D}, pos::SVector{D,T}) where {S,X,T,D}
    imin = min(ga.bbox) + SVector{D,S}(ga.cell_centered) / 2
    imax = max(ga.bbox) - SVector{D,S}(ga.cell_centered;) / 2
    return linterp(imin, min(ga.domain), imax, max(ga.domain),
                   pos)::SVector{D,S}
end

export coordinates!
function coordinates!(dest_ga::GridArray)
    return map!(i -> coordinate(dest_ga, i), dest_ga, iota(dest_ga))
end

export coordinates
function coordinates(ga::GridArray{S,T,D}) where {S,T,D}
    return map(i -> coordinate(ga, i), iota(ga))
end

export index
function index(ga::GridArray{S,T,D}, coord::SVector{D,S}) where {S,T,D}
    imin = min(ga.bbox) + SVector{D,S}(ga.cell_centered) / 2
    imax = max(ga.bbox) - SVector{D,S}(ga.cell_centered;) / 2
    return linterp(min(ga.domain), imin, max(ga.domain), imax,
                   coord)::SVector{D,S}
end

export zero!
zero!(dest_ga::GridArray) = map!(x -> zero(x), dest_ga)
Base.zero(ga::GridArray) = map(x -> zero(eltype(ga)), ga)

function Base.:+(ga1::GridArray{S,T1,D},
                 ga2s::GridArray{S,T2,D}...) where {S,T1,T2,D}
    return map(+, ga1, ga2s...)
end

Base.:-(ga::GridArray) = map(-, ga)
function Base.:-(ga1::GridArray{S,T1,D},
                 ga2::GridArray{S,T2,D}) where {S,T1,T2,D}
    return map(-, ga1, ga2)
end

Base.:*(α, ga::GridArray) = map(x -> α * x, ga)
Base.:*(ga::GridArray, α) = map(x -> x * α, ga)
Base.:\(α, ga::GridArray) = map(x -> α \ x, ga)
Base.:/(ga::GridArray, α) = map(x -> x / α, ga)
