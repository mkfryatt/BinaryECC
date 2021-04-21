#sec1v2 2.3.6
#create a vector from a hex string
function MVector{L,T}(s::String)::MVector{L,T} where {L,T}
    value = zeros(MVector{L,T})
    blocksize = sizeof(T)*2
    for i in 1:L
        range = max(1, length(s)-blocksize*i+1):(length(s)-blocksize*(i-1))
        value[i] = parse(T, s[range], base=16)
    end
    return value
end
function SVector{L,T}(s::String)::SVector{L,T} where {L,T}
    value = zeros(MVector{L,T})
    blocksize = sizeof(T)*2
    for i in 1:L
        range = max(1, length(s)-blocksize*i+1):(length(s)-blocksize*(i-1))
        value[i] = parse(T, s[range], base=16)
    end
    return SVector(value)
end

function isone(x::Union{SVector{L,T},MVector{L,T}})::Bool where {L,T}
    return isone(x[1]) && iszero(x[2:L])
end
function isone(x::Union{SVector{1,T},MVector{1,T}})::Bool where T
    return isone(x[1])
end

function one(::Type{MVector{L,T}})::MVector{L,T} where {L,T}
    value = zeros(MVector{L,T})
    value[1] = T(1)
    return value
end
function one(::Type{SVector{L,T}})::SVector{L,T} where {L,T}
    return SVector(one(MVector{L,T}))
end

#returns the bit (as a 0 or 1 of type T) in position i of the number stored in x
function getbit(x::Union{SVector{L,T},MVector{L,T}}, i::Int)::Int where {L,T}
    blocksize = sizeof(T)*8
    bit = i%blocksize
    block = (i÷(8*sizeof(T))) +1
    if block>L return 0 end
    return x[block]>>bit & 1
end

#TODO doesn't work over block boundaries
function getbits(x::Union{SVector{L,T},MVector{L,T}}, start::Int, len::Int)::Int where {L,T}
    blocksize = sizeof(T)*8
    bit = start%blocksize
    block = (start÷blocksize) +1
    if block>L return 0 end
    bitmask = (1<<len) -1
    return (x[block]>>bit) & bitmask
end

#updates x to have the bit at positon i flipped
function flipbit!(x::MVector{L,T}, i::Integer) where {L,T}
    blocksize = sizeof(T)*8
    bit = i%blocksize
    block = (i÷blocksize) +1
    if block>L
        throw(ArgumentError("Cannot flip bit beyond the size of the number."))
    end
    x[block] ⊻= T(1)<<bit
end

#length in bits of the number stored in x
function bits(x::Union{SVector{L,T},MVector{L,T}})::Int where {L,T}
    block = L
    while x[block]==0
        block -= 1
        if block==0 return 0 end
    end
    return bitsize(T)*block - leading_zeros(x[block])
end

#returns true if the numbers stored in x and y are the same
function ==(x::Union{SVector{L1,T},MVector{L1,T}}, y::Union{SVector{L2,T},MVector{L2,T}})::Bool where {L1,L2,T}
    if L1<L2
        x, y = y, x
    end
    return x[1:L2]==y[1:L2] && iszero(x[L2+1:L1])
end

function ⊻(x::SVector{L,T}, y::SVector{L,T}) where {L,T}
    return SVector{L,T}(x[i]⊻y[i] for i=1:L)
end
function ⊻(x::MVector{L,T}, y::MVector{L,T}) where {L,T}
    return MVector{L,T}(x[i]⊻y[i] for i=1:L)
end

function xor!(x::MVector{L1,T}, y::Union{SVector{L2,T},MVector{L2,T}}) where {L1,L2,T}
    for i in 1:min(L1,L2)
        x[i] ⊻= y[i]
    end
end

function leftshift!(x::MVector{L,T}, shift::Int) where {L,T}
    if shift==0 return
    elseif shift<0 throw(ArgumentError("Shift must be nonnegative."))
    end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    upperbits = shift % blocksize
    lowerbits = blocksize - upperbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    for block in L:-1:2+blockshift
        x[block] = (x[block-blockshift]&lowermask)<<upperbits
        x[block] ⊻= (x[block-blockshift-1]>>lowerbits)&uppermask
    end
    x[1+blockshift] = (x[1+blockshift]&lowermask)<<upperbits
end

#writes x ⊻ (y<<shift) into x
#not guaranteed to produce the same result as doing the shift and xor separately
#(due to overflow of behaviour of y when shifted left)
function shiftedxor!(x::MVector{L1,T}, y::Union{SVector{L2,T},MVector{L2,T}}, shift::Int)::Nothing where {L1,L2,T}
    if shift<0 throw(ArgumentError("Cannot shift by a negative amount")) end
    if shift==0
        xor!(x, y)
        return
    end
    if shift%bitsize(T)==0
        word_shiftedxor!(x, y, shift ÷ bitsize(T))
        return
    end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    upperbits = shift % blocksize

    #otherwise, each new block requires bits from two shifted blocks
    lowerbits = blocksize - upperbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    top = L2
    if L2+blockshift==L1
        top = L2-1
        x[L1] ⊻= (y[L2]&lowermask)<<upperbits
    elseif L2+blockshift>L1
        top = L1-blockshift-1
        x[L1] ⊻= (y[L2-blockshift]&lowermask)<<upperbits
    end
    for block in 1:top
        x[block+blockshift] ⊻= (y[block]&lowermask)<<upperbits
        x[block+blockshift+1] ⊻= (y[block]>>lowerbits)&uppermask
    end
end

#shift by an entire word at a time
function word_shiftedxor!(x::MVector{L1,T}, y::Union{SVector{L2,T},MVector{L2,T}}, wordshift::Int)::Nothing where {L1,L2,T}
    top = L2
    if L2+wordshift>L1
        top = L1-wordshift
    end
    for block in 1:top
        x[block+wordshift] ⊻= y[block]
    end
end

#returns a new vector that has "newL" words
#if newL>L, the extra words are zeros
#if newL<L, it simply chops off the superfluous words
function changelength(x::SVector{L,T}, newL::Int)::SVector{newL,T} where {L,T}
    if newL<1 throw(ArgumentError("New length must be positive")) end
    if L==newL return x
    elseif L>newL return SVector{newL,T}(x[1:newL])
    end

    value::MVector{newL,T} = zeros(MVector{newL,T})
    for i in 1:min(L,newL)
        value[i] = x[i]
    end
    return SVector(value)
end
function changelength(x::MVector{L,T}, newL::Int)::MVector{newL,T} where {L,T}
    if newL<1 throw(ArgumentError("New length must be positive")) end
    if L==newL return x
    elseif L>newL return MVector{newL,T}(x[1:newL])
    end

    value::MVector{newL,T} = zeros(MVector{newL,T})
    for i in 1:min(L,newL)
        value[i] = x[i]
    end
    return value
end

#sec1 v2, 2.3.9
#returns the number stored in a, as a BigInt
function convert(::Type{BigInt}, a::Union{SVector{L,T},MVector{L,T}})::BigInt where {L,T}
    b = BigInt(0)
    for i in 1:L
        b += BigInt(a[i])<<(8*sizeof(T)*(i-1))
    end
    return b
end

#sec1 v2, 2.3.5
#returns the number stored in a as a hex string
function convert(::Type{String}, a::Union{SVector{L,T},MVector{L,T}})::String where {L,T}
    M = ""
    for block in a
        Mi = string(block, base=16)
        Mi = "0"^(2*sizeof(T)-length(Mi)) *Mi
        M = Mi*M
    end
    return M
end

function convert(::Type{MVector{L,T}}, value::Integer)::MVector{L,T} where {L,T}
    x = zeros(MVector{L,T})
    i = 1
    while value>0 && i<L
        x[i] = value & (T(1)<<bitsize(T) -1)
        value >>>= bitsize(T)
        i += 1
    end
    return x
end
function convert(::Type{SVector{L,T}}, value::Integer)::SVector{L,T} where {L,T}
    return SVector(convert(MVector{L,T}, value))
end

#returns a random number of the given type
function random(::Type{MVector{L,T}})::MVector{L,T} where {L,T}
    return MVector{L,T}([rand(T) for i=1:L])
end
function random(::Type{SVector{L,T}})::SVector{L,T} where {L,T}
    return SVector{L,T}([rand(T) for i=1:L])
end

#retuns a random number of the given type, no longer than "bits" bits long
function random(::Type{MVector{L,T}}, bits)::MVector{L,T} where {L,T}
    value = zeros(MVector{L,T})
    blocksize = sizeof(T)*8
    zerowords = L - ceil(Int, bits / blocksize)
    for i in 1:(L-zerowords)
        value[i] = rand(T)
    end
    value[L-zerowords] &= (T(1)<<(bits%blocksize)) -1
    return value
end
function random(::Type{SVector{L,T}}, bits)::SVector{L,T} where {L,T}
    value = zeros(MVector{L,T})
    blocksize = sizeof(T)*8
    zerowords = L - ceil(Int, bits / blocksize)
    for i in 1:(L-zerowords)
        value[i] = rand(T)
    end
    value[L-zerowords] &= (T(1)<<(bits%blocksize)) -1
    return SVector(value)
end
