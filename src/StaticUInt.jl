#L is the number of blocks
#T is the type of each block
struct StaticUInt{L,T<:Unsigned}
    value::MVector{L,T}
    StaticUInt{L,T}(value::MVector{L,T}) where {L,T} = new(value)
end

#create a staticuint from an integer
function StaticUInt{L,T}(value::Integer)::StaticUInt{L,T} where {L,T}
    if value<0 throw(ArgumentError("Expected a nonnegative value.")) end

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
#create a staticuint from a hex string
function StaticUInt{L,T}(s::String)::StaticUInt{L,T} where {L,T}
    value = zeros(MVector{L,T})
    blocksize = sizeof(T)*2
    for i in 1:L
        range = max(1, length(s)-blocksize*i+1):(length(s)-blocksize*(i-1))
        value[i] = parse(T, s[range], base=16)
    end
    return StaticUInt{L,T}(value)
end

#deep copy of a staticuint
function copy(x::StaticUInt{L,T})::StaticUInt{L,T} where {L,T}
    return StaticUInt{L,T}(copy(x.value))
end

function isone(x::StaticUInt{L,T})::Bool where {L,T}
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

function one(::Type{StaticUInt{L,T}})::StaticUInt{L,T} where {L,T}
    value = zeros(MVector{L,T})
    value[1] = T(1)
    return StaticUInt{L,T}(value)
end

function iszero(x::StaticUInt)::Bool
    return iszero(x.value)
end

function zero(::Type{StaticUInt{L,T}})::StaticUInt{L,T} where {L,T}
    return StaticUInt{L,T}(zeros(MVector{L,T}))
end

#returns the bit (as a 0 or 1 of type T) in position i of the number stored in x
function getbit(x::StaticUInt{L,T}, i::Int)::Int where {L,T}
    blocksize = sizeof(T)*8
    bit = i%blocksize
    block = (i÷(8*sizeof(T))) +1
    if block>L return 0 end
    return x.value[block]>>bit & 1
end

function getbits(x::StaticUInt{L,T}, start::Int, len::Int)::Int where {L,T}
    blocksize = sizeof(T)*8
    bit = start%blocksize
    block = (start÷(8*sizeof(T))) +1
    if block>L return 0 end
    bitmask = (1<<len) -1
    return (x.value[block]>>bit) & bitmask
end

#updates x to have the bit at positon i flipped
function flipbit!(x::StaticUInt{L,T}, i::Integer) where {L,T}
    blocksize = sizeof(T)*8
    bit = i%blocksize
    block = (i÷(8*sizeof(T))) +1
    if block>L
        throw(ArgumentError("Cannot flip bit beyond the size of the number."))
    end
    x.value[block] ⊻= T(1)<<bit
end

#length in bits of the number stored in x
function bits(x::StaticUInt{L,T})::Int where {L,T}
    block = L
    while x.value[block]==0
        block -= 1
        if block==0 return 0 end
    end
    return sizeof(T)*8*block - leading_zeros(x.value[block])
end

#returns true if the numbers stored in x and y are the same
function ==(x::StaticUInt{L1,T}, y::StaticUInt{L2,T})::Bool where {L1,L2,T}
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

#returns the result of x xor y, in a staticuint{L,T}, where L = max(L1,L2)
function ⊻(x::StaticUInt{L1,T}, y::StaticUInt{L2,T})::StaticUInt{max(L1,L2),T} where {L1,L2,T}
    if L1<L2
        x, y = y, x
    end
    value = copy(x.value)
    for i in 1:L2
        value[i] ⊻= y.value[i]
    end
    return StaticUInt{L1,T}(value)
end

function xor!(x::StaticUInt{L1,T}, y::StaticUInt{L2,T}) where {L1,L2,T}
    for i in 1:min(L1,L2)
        x.value[i] ⊻= y.value[i]
    end
end

#returns the result of x shifted left by "shift" bits
function <<(x::StaticUInt{L,T}, shift::Int)::StaticUInt{L,T} where {L,T}
    if shift<0 x>>(-shift)
    elseif shift==0 return x
    end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    upperbits = shift % blocksize
    lowerbits = blocksize - upperbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    value = zeros(MVector{L,T})
    value[L] = (x.value[L-blockshift]&lowermask)<<upperbits
    for block in 1:(L-blockshift-1)
        value[block+blockshift] ⊻= (x.value[block]&lowermask)<<upperbits
        value[block+blockshift+1] ⊻= (x.value[block]>>lowerbits)&uppermask
    end

    return StaticUInt{L,T}(value)
end

function leftshift!(x::StaticUInt{L,T}, shift::Int) where {L,T}
    if shift==0 return
    elseif shift<0 throw(ArgumentError("Shift must be nonnegative.")) end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    upperbits = shift % blocksize
    lowerbits = blocksize - upperbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    for block in L:-1:2+blockshift
        x.value[block] = (x.value[block-blockshift]&lowermask)<<upperbits
        x.value[block] ⊻= (x.value[block-blockshift-1]>>lowerbits)&uppermask
    end
    x.value[1+blockshift] = (x.value[1+blockshift]&lowermask)<<upperbits
end

#returns x ⊻ (y<<shift) in a staticuint{L1,T}
#not guaranteed to produce the same result as doing the shift and xor separately
#(due to overflow of behaviour of y when shifted left)
function shiftedxor(x::StaticUInt{L1,T}, y::StaticUInt{L2,T}, shift::Int)::StaticUInt{L1,T} where {L1,L2,T}
    if shift<0 throw(ArgumentError("Cannot shift by a negative amount"))
    elseif shift==0 return x ⊻ y
    end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    value = copy(x.value)

    #if the shift is a whole number of blocks, the calculation is much simpler
    if shift%blocksize==0
        for i in 1:min(L1-blockshift, L2)
            value[i+blockshift] ⊻= y.value[i]
        end
        return StaticUInt{L1,T}(value)
    end

    #otherwise, each new block requires bits from two shifted blocks
    upperbits = shift % blocksize
    lowerbits = blocksize - upperbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    top = L2
    if L2+blockshift==L1
        top = L2-1
        value[L1] ⊻= (y.value[L2]&lowermask)<<upperbits
    elseif L2+blockshift>L1
        top = L1-blockshift-1
    end
    for block in 1:top
        value[block+blockshift] ⊻= (y.value[block]&lowermask)<<upperbits
        value[block+blockshift+1] ⊻= (y.value[block]>>lowerbits)&uppermask
    end

    return StaticUInt{L1,T}(value)
end

function shiftedxor!(x::StaticUInt{L1,T}, y::StaticUInt{L2,T}, shift::Int)::Nothing where {L1,L2,T}
    if shift<0 throw(ArgumentError("Cannot shift by a negative amount")) end

    #if shift is zero, this is just xor that doesn't perform a copy
    if shift==0
        for i in 1:min(L1, L2)
            x.value[i] ⊻= y.value[i]
        end
        return
    end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    upperbits = shift % blocksize

    #if the shift is a whole number of blocks
    #the calculation is simpler because you don't need bitmasks
    if upperbits==0
        for i in 1:min(L1-blockshift, L2)
            x.value[i+blockshift] ⊻= y.value[i]
        end
        return
    end

    #otherwise, each new block requires bits from two shifted blocks
    lowerbits = blocksize - upperbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    top = L2
    if L2+blockshift==L1
        top = L2-1
        x.value[L1] ⊻= (y.value[L2]&lowermask)<<upperbits
    elseif L2+blockshift>L1
        top = L1-blockshift-1
    end
    for block in 1:top
        x.value[block+blockshift] ⊻= (y.value[block]&lowermask)<<upperbits
        x.value[block+blockshift+1] ⊻= (y.value[block]>>lowerbits)&uppermask
    end
end

#returns the result of right shifting the number stored in x by "shift" bits
function >>(x::StaticUInt{L,T}, shift::Int)::StaticUInt{L,T} where {L,T}
    if shift<0 x<<(-shift)
    elseif shift==0 return x
    end

    blocksize = 8*sizeof(T)
    blockshift = shift ÷ blocksize
    lowerbits = shift % blocksize
    upperbits = blocksize - lowerbits
    lowermask = (T(1)<<lowerbits)-1
    uppermask = (T(1)<<upperbits)-1

    value = zeros(MVector{L,T})
    value[1] = (x.value[1+blockshift]>>lowerbits)&uppermask
    for block in (2+blockshift):L
        value[block-blockshift] ⊻= (x.value[block]>>lowerbits)&uppermask
        value[block-blockshift-1] ⊻= (x.value[block]&lowermask)<<upperbits
    end
    return StaticUInt{L,T}(value)
end

#returns the length of x in words
function length(x::StaticUInt{L,T})::Int where {L,T}
    return L
end

#returns a new staticuint that has "newL" words
#if newL>L, the extra words are zeros
#if newL<L, it simply chops off the superfluous words
function changelength(x::StaticUInt{L,T}, newL::Int)::StaticUInt{newL,T} where {L,T}
    if newL<1 throw(ArgumentError("New length must be positive")) end
    if L==newL return copy(x) end

    value::MVector{newL,T} = zeros(MVector{newL,T})
    for i in 1:min(L,newL)
        value[i] = x.value[i]
    end
    return StaticUInt{newL,T}(value)
end

#sec1 v2, 2.3.9
#returns the number stored in a, as a BigInt
function convert(::Type{BigInt}, a::StaticUInt{L,T})::BigInt where {L,T}
    b = BigInt(0)
    for i in 1:L
        b += BigInt(a.value[i])<<(8*sizeof(T)*(i-1))
    end
    return b
end

#sec1 v2, 2.3.5
#returns the number stored in a as a hex string
function convert(::Type{String}, a::StaticUInt{L,T})::String where {L,T}
    M = ""
    for block in a.value
        Mi = string(block, base=16)
        Mi = "0"^(2*sizeof(T)-length(Mi)) *Mi
        M = Mi*M
    end
    return M
end

#returns a random number of the given type
function random(::Type{StaticUInt{L,T}})::StaticUInt{L,T} where {L,T}
    value = MVector{L,T}([rand(T) for i=1:L])
    return StaticUInt{L,T}(value)
end

#retuns a random number of the given type, no longer than "bits" bits long
function random(::Type{StaticUInt{L,T}}, bits)::StaticUInt{L,T} where {L,T}
    value = zeros(MVector{L,T})
    blocksize = sizeof(T)*8
    zerowords = L - ceil(Int, bits / blocksize)
    for i in 1:(L-zerowords)
        value[i] = rand(T)
    end
    value[L-zerowords] &= (T(1)<<(bits%blocksize)) -1
    return StaticUInt{L,T}(value)
end
