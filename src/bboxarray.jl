# BBoxArray

export BBoxArray
struct BBoxArray{D} <: BBoxSet{D}
    bboxes::Vector{BBox{D}}
    BBoxArray{D}() where {D} = new(BBox{D}[])
    BBoxArray{D}(bboxes::Vector{BBox{D}}) where {D} = new(bboxes)
end

BBoxArray(bboxes::Vector{BBox{D}}) where {D} = BBoxArray{D}(bboxes)

function invariant(bbarray::BBoxArray)
    any(isempty, bbarray.bboxes) && return false
    for i in 1:length(bbarray.bboxes), j in (i + 1):length(bbarray.bboxes)
        intersects(bbarray.bboxes[i], bbarray.bboxes[j]) && return false
    end
    return true
end

Base.ndims(bbarray::BBoxArray{D}) where {D} = D

Base.eltype(bbarray::BBoxArray) = Int

Base.isempty(bbarray::BBoxArray) = isempty(bbarray.bboxes)

Base.length(bbarray::BBoxArray) = sum(length, bbarray.bboxes)

function shift!(bbarray::BBoxArray{D}, offset::SVector{D,Int}) where {D}
    for i in 1:length(bbarray.bboxes)
        bbarray.bboxes[i] = shift(bbarray.bboxes[i], offset)
    end
    return bbarray
end

function shift(bbarray::BBoxArray{D}, offset::SVector{D,Int}) where {D}
    return shift!(BBoxArray(copy(bbarray.bboxes)))
end

function Base.in(pos::SVector{D,Int}, bbarray::BBoxArray{D}) where {D}
    return any(bbox -> pos in bbox, bbarray.bboxes)
end

# in, intersect, union, setdiff, all with BBox and BBoxArray
