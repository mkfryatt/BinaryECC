#L is the number of blocks
#T is the type of each block
struct StaticUInt{L,T<:Unsigned}
    value::SVector{L,T}
    StaticUInt{L,T}(value::SVector{L,T}) where {L,T} = new(value)
end

function StaticUInt{L,T}(value::Integer) where {L,T}
    if value<0 throw(ArgumentError("Expected a nonnegative value.")) end
    blocksize = sizeof(T)*8
    bitmask = typemax(T)
    generator = (T((value>>(blocksize*i)) & bitmask) for i=0:(L-1))
    return StaticUInt{L,T}(SVector{L,T}(generator))
end

#sec1v2 2.3.6
#convert a hex string to a field element
function StaticUInt{L,T}(s::String) where {L,T}
    blocksize = sizeof(T)*2
    generator = (begin
        range = max(1, length(s)-blocksize*i+1):(length(s)-blocksize*(i-1))
        parse(T, s[range], base=16)
    end for i=1:L)
    return StaticUInt{L,T}(SVector{L,T}(generator))
end

#TODO is this needed?
function copy(x::StaticUInt{L,T}) where {L,T}
    return StaticUInt{L,T}(x.value)
end

function isone(x::StaticUInt{L,T}) where {L,T}
    if x.value[1]!=1
        return false
    end
    for i in 2:L
        if x.value[i]!=0
            return false
        end
    end
    return true
end

function one(::Type{StaticUInt{L,T}}) where {L,T}
    generator = (if i==1 T(1) else T(0) end for i=1:L)
    return StaticUInt{L,T}(SVector{L,T}(generator))
end

function iszero(x::StaticUInt)
    return iszero(x.value)
end

function zero(::Type{StaticUInt{L,T}}) where {L,T}
    return StaticUInt{L,T}(zeros(SVector{L,T}))
end

function getbit(x::StaticUInt{L,T}, i::Integer) where {L,T}
    blocksize = sizeof(T)*8
    bit = i%blocksize
    block = (i÷(8*sizeof(T))) +1
    if block>L return 0 end
    return x.value[block]>>bit & 1
end

function flipbit(x::StaticUInt{L,T}, i::Integer) where {L,T}
    blocksize = sizeof(T)*8
    bit = i%blocksize
    block = (i÷(8*sizeof(T))) +1
    if block>L throw(ArgumentError("Cannot flip bit beyond the size of the number.")) end
    generator = (if i==block x.value[block] ⊻ (T(1)<<bit) else x.value[i] end for i=1:L)
    return StaticUInt{L,T}(SVector{L,T}(generator))
end

#length of a number in bits
function bits(x::StaticUInt{L,T}) where {L,T}
    block = L
    while x.value[block]==0
        block -= 1
        if block==0 return 0 end
    end
    return sizeof(T)*8*block - leading_zeros(x.value[block])
end

function ==(x::StaticUInt{L,T}, y::StaticUInt{L,T}) where {L,T}
    return x.value==y.value
end

function ==(x::StaticUInt{L1,T}, y::StaticUInt{L2,T}) where {L1,L2,T}
    if L1<L2
        x, y = y, x
    end
    for i in 1:L2
        if x.value[i]!=y.value[i]
            return false
        end
    end
    for i in (L2+1):L1
        if !iszero(x.value[i])
            return false
        end
    end
    return true
end

function ⊻(x::StaticUInt{L1,T}, y::StaticUInt{L2,T}) where {L1,L2,T}
    if L1<L2
        x, y = y, x
    end
    generator = (x.value[i] ⊻ y.value[i] for i=1:L1)
    return StaticUInt{L1,T}(SVector{L1,T}(generator))
end

function <<(x::StaticUInt{L,T}, shift::Integer) where {L,T}
    if shift<0 x>>(-shift)
    elseif shift==0 return x
    end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    lowerbits = shift % blocksize
    upperbits = blocksize - lowerbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    generator = (
    if i<=blockshift
        T(0)
    else
        if i-blockshift>1
            ((x.value[i-blockshift]&uppermask)<<lowerbits) ⊻
            ((x.value[i-blockshift-1]>>upperbits)&lowermask)
        else
            (x.value[i-blockshift]&uppermask)<<lowerbits
        end
    end
    for i=1:L
    )

    return StaticUInt{L,T}(SVector{L,T}(generator))
end

function length(x::StaticUInt{L,T}) where {L,T}
    return L
end

function changelength(x::StaticUInt{L,T}, newL::Integer) where {L,T}
    if newL<1 throw(ArgumentError("New length must be positive")) end
    generator = (
    if i<=L
        x.value[i]
    else
        T(0)
    end for i=1:newL)
    return StaticUInt{newL,T}(SVector{newL,T}(generator))
end

#sec1 v2, 2.3.9
function convert(::Type{BigInt}, a::StaticUInt{L,T}) where {L,T}
    b = BigInt(0)
    for i in 1:L
        b += BigInt(a.value[i])<<(8*sizeof(T)*(i-1))
    end
    return b
end

#sec1 v2, 2.3.5
function convert(::Type{String}, a::StaticUInt{L,T}) where {L,T}
    M = ""
    for block in a.value
        Mi = string(block, base=16)
        Mi = "0"^(2*sizeof(T)-length(Mi)) *Mi
        M = Mi*M
    end
    return M
end

function random(::Type{StaticUInt{L,T}}) where {L,T}
    generator = (rand(T) for i=1:L)
    return StaticUInt{L,T}(SVector{L,T}(generator))
end
