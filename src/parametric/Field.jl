import Base: +, -, *, /, ^, ==, repr, inv, sqrt, iszero

#D is the degree of the reduction polynomial
#R is the reduction polynomial without the x^D term
struct FieldPoint{D,R}
    value::BigInt
    FieldPoint{D,R}(value::Integer) where {D,R} = new(convert(BigInt, value))
end

#sec1v2 2.3.6
#convert a hex string to a field element
function FieldPoint{D,R}(s::String) where {D,R}
    s = replace(s, " " => "")
    if length(s)!=ceil(D / 8)*2 throw(ArgumentError("Octet string is of the incorrect length for this field.")) end
    value = parse(BigInt, s, base=16)
    return FieldPoint{D,R}(value)
end

function repr(a::FieldPoint)
    return repr(a.value)
end

function ==(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    return a.value==b.value
end

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
    if a.value<(BigInt(1)<<D) return a end

    #k is the number of bits of a that need to be removed from a
    k = bits(a.value) - D

    #shift the reduction polynomial left
    #the loop will slowly shift it back down again
    t = BigInt(R) << k

    #b will eventually be such that a ≡ b (mod R)
    b = a.value

    for i in (k+D):-1:D
        if b & (BigInt(1) << i) != BigInt(0)
            b ⊻= t + (BigInt(1)<<i)
        end
        t >>>= 1
    end

    return FieldPoint{D,R}(b)
end

#right to left, shift and add
function *(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    if a.value==b.value return square(a) end

    c = BigInt(0)
    shiftedb = b.value

    for i in 0:(D-1)
        if a.value & (BigInt(1)<<i) != BigInt(0)
            c ⊻= shiftedb
        end
        shiftedb <<= 1
    end

    return reduce(FieldPoint{D,R}(c))
end

#add a zero between every digit of the original
function square(a::FieldPoint{D,R}) where {D,R}
    b = BigInt(0)
    counter = BigInt(1)
    for i in 0:(D-1)
        if a.value & counter != BigInt(0)
            b += BigInt(1) << (i*2)
        end
        counter <<= 1
    end

    return reduce(FieldPoint{D,R}(b))
end

#length of a number in bits
function bits(a::Integer)
    i = 0
    while a > (BigInt(1)<<i)
        i += 1
    end
    if a == (BigInt(1)<<i)
        return i+1
    else
        return i
    end
    #return floor(Int, log2(a)) +1 #log is an approximation
    #TODO find a constant time method for this?
end

#uses a version of egcd to invert a
#Algorithm 2.48, Guide to Elliptic Curve Cryptography
function inv(a::FieldPoint{D,R}) where {D,R}
    if a.value==0 throw(DivideError()) end

    u = a.value
    v = BigInt(R) + (BigInt(1)<<D)
    g1 = BigInt(1)
    g2 = BigInt(0)

    while u!=BigInt(1)
        j = bits(u) - bits(v)
        if j<0
            (u, v) = (v, u)
            (g1, g2) = (g2, g1)
            j = -j
        end
        u ⊻= v << j
        g1 ⊻= g2 << j
    end
    return reduce(FieldPoint{D,R}(g1))
end

function /(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    return a * inv(b)
end

#right to left, square and multiply method
function ^(a::FieldPoint{D,R}, b::Integer) where {D,R}
    c = FieldPoint{D,R}(1)
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

function sqrt(a::FieldPoint{D,R}) where {D,R}
    return a^(BigInt(1)<<(D-1))
end

function random(::Type{FieldPoint{D,R}}) where {D,R}
    range = BigInt(0):((BigInt(1)<<D)-BigInt(1))
    return FieldPoint{D,R}(rand(range))
end

function iszero(a::FieldPoint)
    return a.value==0
end

#sec2 v2, table 3:
FieldPoint163 = FieldPoint{163, Int128(128+64+8+1)}
FieldPoint233 = FieldPoint{233, (Int128(1)<<74) + Int128(1)}
FieldPoint239 = FieldPoint{239, (Int128(1)<<36) + Int128(1)}
FieldPoint283 = FieldPoint{283, (Int128(1)<<12) + Int128(128+32+1)}
FieldPoint409 = FieldPoint{409, (Int128(1)<<87) + Int128(1)}
FieldPoint571 = FieldPoint{571, (Int128(1)<<10) + Int128(32+4+1)}
