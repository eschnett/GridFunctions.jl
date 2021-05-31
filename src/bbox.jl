# BBox

export BBox
"""
Bounding box for a region of integer points. Lower and upper bounds
are inclusive.
"""
struct BBox{D}
    min::SVector{D,Int}
    max::SVector{D,Int}
    BBox{D}() where {D} = new{D}(zero(SVector{D,Int}), -one(SVector{D,Int}))
    function BBox{D}(min::SVector{D,Int}, max::SVector{D,Int}) where {D}
        return new{D}(min, max)
    end
end

BBox(min::SVector{D,Int}, max::SVector{D,Int}) where {D} = BBox{D}(min, max)

Base.min(bbox::BBox) = bbox.min
Base.max(bbox::BBox) = bbox.max

Base.ndims(bbox::BBox{D}) where {D} = D

Base.eltype(bbox::BBox) = Int

Base.isempty(bbox::BBox) = any(max(bbox) > min(bbox))

function Base.size(bbox::BBox{D}) where {D}
    isempty(bbox.max) && return zero(SVector{D,Int})
    return max(bbox) - min(bbox) .+ 1
end

Base.length(bbox::BBox) = prod(size(bbox))

function Base.strides(bbox::BBox{D}) where {D}
    res = zero(SVector{D,Int})
    str = 1
    for d in 1:D
        res = setindex!(res, str, d)
        str *= size(bbox, d)
    end
    return res
end

export linear_index
function linear_index(bbox::BBox{D}, pos::SVector{D,Int}) where {D}
    return sum((pos - min(bbox)) .* strides(bbox))
end

export shift
function shift(bbox::BBox{D}, offset::SVector{D,Int}) where {D}
    isempty(bbox) && return BBox{D}()
    return BBox{D}(min(bbox) + offset, max(bbox) + offset)
end

export extend
function extend(bbox::BBox{D}, lo::SVector{D,Int}, hi::SVector{D,Int}) where {D}
    isempty(bbox) && return BBox{D}()
    return BBox{D}(min(bbox) - lo, max(bbox) + hi)
end

function Base.in(pos::SVector{D,Int}, bbox::BBox{D}) where {D}
    isempty(bbox) && return false
    return all(min(bbox) .≤ pos .≤ max(bbox))
end

function Base.intersect(bbox1::BBox{D}, bbox2::BBox{D}) where {D}
    (isempty(bbox1) || isempty(bbox2)) && return BBox{D}()
    return BBox{D}(max(min(bbox1), min(bbox2)), min(max(bbox1), max(bbox2)))
end

export intersects
function intersects(bbox1::BBox{D}, bbox2::BBox{D}) where {D}
    return !isempty(intersect(bbox1, bbox2))
end
