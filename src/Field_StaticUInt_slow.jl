#D is the degree of the reduction polynomial
#R is the reduction polynomial without the x^D term
struct FieldPoint{D,R}
    value::StaticUInt
    FieldPoint{D,R}(value::Integer) where {D,R} = new(StaticUInt{ceil(Int,D/64),UInt64}(value))
    FieldPoint{D,R}(value::StaticUInt) where {D,R} = new(value)
end

#sec1v2 2.3.6
#convert a hex string to a field element
function FieldPoint{D,R}(s::String) where {D,R}
    s = replace(s, " " => "")
    if length(s)!=ceil(D / 8)*2 throw(ArgumentError("Octet string is of the incorrect length for this field.")) end
    value = StaticUInt{ceil(Int,D/64),UInt64}(s)
    return FieldPoint{D,R}(value)
end

function ==(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    return a.value==b.value
end

"""
    +(x1::FieldPoint{D,R}, x2::FieldPoint{D,R}) where {D,R}
Addition of elements x1, x2 in the field represented by D and R.
"""
function +(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    return FieldPoint{D,R}(a.value ⊻ b.value)
end

function -(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    return a+b
end

function -(a::FieldPoint{D,R}) where {D,R}
    return a
end

function reduce(a::FieldPoint{D,R}) where {D,R}
    #b will should always be such that a ≡ b (mod R)
    #the loop will modify it until it reaches the smallest value that makes that true
    b = copy(a.value)
    r = typeof(a.value)(R)

    #iterate over the excess bits of a, left to right
    for i in (64*length(b)-1):-1:D
        if getbit(b, i)==1
            b = flipbit(b, i)
            b ⊻= r << (i-D)
        end
    end

    #remove excess blocks from b
    b = changelength(b, ceil(Int,D/64))
    return FieldPoint{D,R}(b)
end

#right to left, shift and add
function *(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    if a.value==b.value return square(a) end

    c = zero(StaticUInt{length(a.value)*2,UInt64})
    shiftedb = changelength(b.value, 2*length(a.value))

    for i in 0:(D-1)
        if getbit(a.value,i)==1
            c ⊻= shiftedb
        end
        shiftedb <<= 1
    end

    return reduce(FieldPoint{D,R}(c))
end

#add a zero between every digit of the original
function square(a::FieldPoint{D,R}) where {D,R}
    b = zero(StaticUInt{ceil(Int, D/32),UInt64})
    for i in 0:(D-1)
        if getbit(a.value,i)==1
            b = flipbit(b, i*2)
        end
    end

    return reduce(FieldPoint{D,R}(b))
end

#uses a version of egcd to invert a
#Algorithm 2.48, Guide to Elliptic Curve Cryptography
function inv(a::FieldPoint{D,R}) where {D,R}
    if iszero(a.value) throw(DivideError()) end

    u = a.value
    v = flipbit(typeof(a.value)(R), D)
    g1 = one(typeof(a.value))
    g2 = zero(typeof(a.value))

    while !isone(u)
        j = bits(u) - bits(v)
        if j<0
            u, v = v, u
            g1, g2 = g2, g1
            j = -j
        end
        u ⊻= v << j
        g1 ⊻= g2 << j
    end
    return FieldPoint{D,R}(g1)
end

function /(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    return a * inv(b)
end

#right to left, square and multiply method
function ^(a::FieldPoint{D,R}, b::Integer) where {D,R}
    c = one(typeof(a))
    squaring = a

    while b>BigInt(0)
        if b & BigInt(1) == BigInt(1)
            c *= squaring
        end
        squaring *= squaring
        b >>>= 1
    end

    return c
end

function random(::Type{FieldPoint{D,R}}) where {D,R}
    return FieldPoint{D,R}(random(StaticUInt{ceil(Int,D/64),UInt64}, D-1))
end

function iszero(a::FieldPoint)
    return iszero(a.value)
end

function zero(::Type{FieldPoint{D,R}}) where {D,R}
    return FieldPoint{D,R}(zero(StaticUInt{ceil(Int,D/64),UInt64}))
end

function isone(a::FieldPoint)
    return isone(a.value)
end

function one(::Type{FieldPoint{D,R}}) where {D,R}
    return FieldPoint{D,R}(one(StaticUInt{ceil(Int,D/64),UInt64}))
end

function convert(::Type{BigInt}, a::FieldPoint)
    return convert(BigInt, a.value)
end

#sec2 v2 (and v1), table 3:
FieldPoint113 = FieldPoint{113, Int128(512+1)} #v1 only
FieldPoint131 = FieldPoint{131, Int128(256+8+4+1)} #v1 only
FieldPoint163 = FieldPoint{163, Int128(128+64+8+1)}
FieldPoint193 = FieldPoint{193, (Int128(1)<<15) + Int128(1)} #v1 only
FieldPoint233 = FieldPoint{233, (Int128(1)<<74) + Int128(1)}
FieldPoint239 = FieldPoint{239, (Int128(1)<<36) + Int128(1)}
FieldPoint283 = FieldPoint{283, (Int128(1)<<12) + Int128(128+32+1)}
FieldPoint409 = FieldPoint{409, (Int128(1)<<87) + Int128(1)}
FieldPoint571 = FieldPoint{571, (Int128(1)<<10) + Int128(32+4+1)}
