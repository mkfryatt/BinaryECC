import Base: +, -, *, /, ^, ==, repr, inv, sqrt, iszero, convert

struct FieldMismatchException <: Exception end

struct Field
    degree::Int16 #the degree of the reduction polynomial
    reduction::BigInt #the reduction polynomial
    Field(d::Integer, r::Integer) = new(convert(Int16, d), convert(BigInt, r))
end

function ==(a::Field, b::Field)
    return a.degree==b.degree && a.reduction==b.reduction
end

struct FieldPoint
    value::BigInt
    field::Field
    FieldPoint(x::Integer, field::Field) = reduce(new(convert(BigInt, x), field))
end

#sec1v2 2.3.6
#convert a hex string to a field element
function FieldPoint(s::String, f::Field)
    s = replace(s, " " => "")
    if length(s)!=ceil(f.degree / 8)*2 throw(ArgumentError("Octet string is of the incorrect length for this field.")) end
    value = parse(BigInt, s, base=16)
    return FieldPoint(value, f)
end

function repr(a::FieldPoint)
    return repr(a)
end

function ==(a::FieldPoint, b::FieldPoint)
    return a.value==b.value && a.field==b.field
end

function +(a::FieldPoint, b::FieldPoint)
    if a.field!=b.field throw(FieldMismatchException()) end
    return FieldPoint(a.value ⊻ b.value, a.field)
end

function -(a::FieldPoint, b::FieldPoint)
    return a+b
end

function -(a::FieldPoint)
    return a
end

function reduce(a::FieldPoint)
    if a.value<(BigInt(1)<<a.field.degree) return a end

    #k is the number of bits of a that need to be removed from a
    k = bits(a.value) - a.field.degree

    #shift the reduction polynomial left
    #the loop will slowly shift it back down again
    t = a.field.reduction << k

    #b will eventually be such that a ≡ b (mod reduction polynomial)
    b = a.value

    #left to right, get rid of the extra bits of b
    for i in (k+a.field.degree):-1:a.field.degree
        if b & (BigInt(1)<<i) != BigInt(0)
            b ⊻= t + (BigInt(1)<<i)
        end
        t >>>= 1
    end

    return FieldPoint(b, a.field)
end

#right to left, shift and add
function *(a::FieldPoint, b::FieldPoint)
    if a.field!=b.field throw(FieldMismatchException()) end
    if a.value==b.value return square(a) end

    c = BigInt(0)
    shiftedb = b.value

    for i in 0:(a.field.degree-1)
        if a.value & (BigInt(1)<<i) != BigInt(0)
            c ⊻= shiftedb
        end
        shiftedb <<= 1
    end

    return reduce(FieldPoint(c, a.field))
end

#number of bits in the binary representation of this number
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
end

#uses a version of egcd to invert a
#Algorithm 2.48, Guide to Elliptic Curve Cryptography
function inv(a::FieldPoint)
    if a.value==0 throw(DivideError()) end

    u = a.value
    v = a.field.reduction + (BigInt(1)<<a.field.degree)
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
    return reduce(FieldPoint(g1, a.field)) #TODO is this reduce needed?
end

function /(a::FieldPoint, b::FieldPoint)
    if a.field!=b.field throw(FieldMismatchException()) end
    return a * inv(b)
end

#add a zero between every digit of the original
function square(a::FieldPoint)
    b = BigInt(0)
    counter = BigInt(1)
    for i in 0:(a.field.degree-1)
        if a.value & counter != BigInt(0)
            b += BigInt(1) << (i*2)
        end
        counter <<= 1
    end

    return reduce(FieldPoint(b, a.field))
end

#right to left, square and multiply method
function ^(a::FieldPoint, b::Integer)
    c = FieldPoint(1, a.field)
    squaring = a

    while b>BigInt(0)
        if b & BigInt(1) == BigInt(1)
            c *= squaring
        end
        squaring = square(squaring)
        b >>>= 1
    end

    return c
end

function sqrt(a::FieldPoint)
    return a^(BigInt(1)<<(a.field.degree-1))
end

#return a random element of the specified field
function random(f::Field)
    range = BigInt(0):((BigInt(1)<<f.degree)-BigInt(1))
    return FieldPoint(rand(range), f)
end

function iszero(a::FieldPoint)
    return a.value==0
end

#sec1 v2, 2.3.9
function convert(::Type{BigInt}, a::FieldPoint)
    return a.value
end

#sec2 v2 (and v1), table 3:
const FIELD113 = Field(113, BigInt(512+1)) #v1 only
const FIELD131 = Field(131, BigInt(256+8+4+1)) #v1 only
const FIELD163 = Field(163, BigInt(128+64+8+1))
const FIELD193 = Field(193, (BigInt(1)<<15) + BigInt(1)) #v1 only
const FIELD233 = Field(233, (BigInt(1)<<74) + BigInt(1))
const FIELD239A = Field(239, (BigInt(1)<<36) + BigInt(1))
const FIELD239B = Field(239, (BigInt(1)<<158) + BigInt(1))
const FIELD283 = Field(283, (BigInt(1)<<12) + BigInt(128+32+1))
const FIELD409 = Field(409, (BigInt(1)<<87) + BigInt(1))
const FIELD571 = Field(571, (BigInt(1)<<10) + BigInt(32+4+1))
