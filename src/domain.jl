# Domain

export Domain
"""
Bounding box for a domain. The domain does not need to be recangular.
"""
struct Domain{S,D}
    min::SVector{D,S}
    max::SVector{D,S}
end

Base.ndims(::Domain{S,D}) where {S,D} = D

Base.eltype(::Domain{S}) where {S} = S

Base.min(dom::Domain) = dom.min
Base.max(dom::Domain) = dom.max
