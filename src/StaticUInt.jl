#L is the number of blocks
#T is the type of each block
struct StaticUInt{L,T<:Unsigned}
    value::MVector{L,T}
    StaticUInt{L,T}(value::MVector{L,T}) where {L,T} = new(value)
end

function StaticUInt{L,T}(value::Integer) where {L,T}
    if value<0
        throw(ArgumentError("Expected a nonnegative value."))
    end
    blocksize = sizeof(T)*8
    valuevector = zeros(MVector{L,T})
    bitmask = typemax(T)
    for i in 1:L
        valuevector[i] = T(value & bitmask)
        value >>= blocksize
    end
    return StaticUInt{L,T}(valuevector)
end

#sec1v2 2.3.6
#convert a hex string to a field element
function StaticUInt{L,T}(s::String) where {L,T}
    value = zeros(MVector{L,T})
    blocksize = sizeof(T)*2
    for i in 1:L
        range = max(1, length(s)-blocksize*i+1):(length(s)-blocksize*(i-1))
        value[i] = parse(T, s[range], base=16)
    end
    return StaticUInt{L,T}(value)
end

function copy(x::StaticUInt{L,T}) where {L,T}
    return StaticUInt{L,T}(copy(x.value))
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
    value = zeros(MVector{L,T})
    value[1] = T(1)
    return StaticUInt{L,T}(value)
end

function iszero(x::StaticUInt)
    return iszero(x.value)
end

function zero(::Type{StaticUInt{L,T}}) where {L,T}
    return StaticUInt{L,T}(zeros(MVector{L,T}))
end

function getbit(x::StaticUInt{L,T}, i::Integer) where {L,T}
    blocksize = sizeof(T)*8
    bit = i%blocksize
    block = (i÷(8*sizeof(T))) +1
    if block>L
        return 0
    end
    return x.value[block]>>bit & 1
end

function flipbit!(x::StaticUInt{L,T}, i::Integer) where {L,T}
    blocksize = sizeof(T)*8
    bit = i%blocksize
    block = (i÷(8*sizeof(T))) +1
    if block>L
        throw(ArgumentError("Cannot flip bit beyond the size of the number."))
    end
    x.value[block] ⊻= T(1)<<bit
end

#length of a number in bits
function bits(x::StaticUInt{L,T}) where {L,T}
    block = L
    while x.value[block]==0
        block -= 1
        if block==0
            return 0
        end
    end

    blocksize = sizeof(T)*8
    for i in 0:(blocksize-1)
        if x.value[block] == (T(1)<<i)
            return i+1 + blocksize*(block-1)
        elseif x.value[block] < (T(1)<<i)
            return i + blocksize*(block-1)
        end
    end
    return blocksize*block
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
    value = copy(x.value)
    for i in 1:L2
        value[i] ⊻= y.value[i]
    end
    return StaticUInt{L1,T}(value)
end

function <<(x::StaticUInt{L,T}, shift::Integer) where {L,T}
    if shift<0
        throw(ArgumentError("Cannot shift by a negative amount"))
    elseif shift==0
        return x
    end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    upperbits = shift % blocksize
    lowerbits = blocksize - upperbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    value = zeros(MVector{L,T})

    for block in 1:(L-blockshift)
        value[block+blockshift] += (x.value[block]&lowermask)<<upperbits
        if block+blockshift<L
            value[block+blockshift+1] += (x.value[block]>>lowerbits)&uppermask
        end
    end

    return StaticUInt{L,T}(value)
end

#returns x ⊻ (y<<shift)
#TODO
function shiftedxor(x::StaticUInt{L1,T}, y::StaticUInt{L2,T}, shift::Integer) where {L1,L2,T}
    if shift<0
        throw(ArgumentError("Cannot shift by a negative amount"))
    elseif shift==0
        return x ⊻ y
    end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    upperbits = shift % blocksize
    lowerbits = blocksize - upperbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    value = copy(a.value)

    for block in 1:(L1-blockshift)
        value[block+blockshift] ⊻= (x.value[block]&lowermask)<<upperbits
        if block+blockshift<L
            value[block+blockshift+1] ⊻= (x.value[block]>>lowerbits)&uppermask
        end
    end

    return StaticUInt{L1,T}(value)
end

function >>(x::StaticUInt{L,T}, shift::Integer) where {L,T}
    if shift<0
        throw(ArgumentError("Cannot shift by a negative amount"))
    elseif shift==0
        return x
    end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    lowerbits = shift % blocksize
    upperbits = blocksize - lowerbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    value = zeros(MVector{L,T})

    for block in (1+blockshift):L
        value[block-blockshift] = (x.value[block]>>lowerbits)&uppermask
        if block-blockshift>1
            value[block-blockshift-1] += (x.value[block]&lowermask)<<upperbits
        end
    end
    return StaticUInt{L,T}(value)
end

function length(x::StaticUInt{L,T}) where {L,T}
    return L
end

function changelength(x::StaticUInt{L,T}, newL::Integer) where {L,T}
    if newL<1
        throw(ArgumentError("New length must be positive"))
    end
    if L==newL
        return x
    end
    value = zeros(MVector{newL,T})
    for i in 1:min(L,newL)
        value[i] = x.value[i]
    end
    return StaticUInt{newL,T}(value)
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
    value = zeros(MVector{L,T})
    for i in 1:L
        value[i] = rand(T)
    end
    return StaticUInt{L,T}(value)
end
