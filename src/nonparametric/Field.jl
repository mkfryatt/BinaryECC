import Base: +, -, *, /, ^, ==, repr, inv, sqrt, iszero
using Base: log2, floor, rand, ceil

struct FieldMismatchException <: Exception end

struct Field
    order::Int16 #the degree of the polynomial
    reduction::BigInt #the reduction polynomial without the most significant bit
    Field(m::Integer, f::Integer) = new(convert(Int16, m), convert(BigInt, f))
end

function repr(f::Field)
    return "Degree: "*repr(f.order)*"\nReduction: "*repr(f.reduction)
end

function ==(a::Field, b::Field)
    return a.order==b.order && a.reduction==b.reduction
end

struct FieldPoint
    value::BigInt
    field::Field
    FieldPoint(value::Integer, field::Field) = new(convert(BigInt, value), field)
end

#sec1v2 2.3.6
function FieldPoint(s::String, f::Field)
    s = replace(s, " " => "")
    if length(s)!=ceil(f.order / 8)*2 throw(ArgumentError("Octet string is of the incorrect length for this field.")) end
    value = parse(BigInt, s, base=16)
    return FieldPoint(value, f)
end

function repr(a::FieldPoint)
    acc = ""
    i = 0
    val = a.value
    while val > 0
        if val & BigInt(1) > BigInt(0)
            if i==0
                acc = "1"
            elseif acc==""
                acc = "x^"*repr(i)
            else
                acc = "x^"*repr(i)*" + "*acc
            end
        end
        i += 1
        val >>>= BigInt(1)
    end
    return acc
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
    if a.value<(BigInt(1)<<a.field.order) return a end

    #k = order(a) - a.field.order
    #i.e. k is the number of bits of a that need to be reduced
    k = bits(a.value) - a.field.order

    #reduction(x) = x^order + t(x)
    #now shift t left, and the loop will slowly shift it back down again
    t = a.field.reduction << k

    result = a.value

    for i in k:-1:0
        if result & (BigInt(1) << (a.field.order+i)) != BigInt(0)
            result ⊻= t + (BigInt(1)<<(a.field.order+i))
        end
        t >>>= 1
    end

    return FieldPoint(result, a.field)
end

#right to left, shift and add
function *(a::FieldPoint, b::FieldPoint)
    if a.field!=b.field throw(FieldMismatchException()) end

    if a.value & BigInt(1) != BigInt(0)
        c = b.value
    else
        c = BigInt(0)
    end

    for i in 1:(a.field.order-1)
        if a.value & (BigInt(1)<<i) != BigInt(0)
            temp = reduce(FieldPoint(b.value << i, b.field))
            c ⊻= temp.value
        end
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
    #return floor(Int, log2(a)) +1 #log is an approximation
end

function inv(a::FieldPoint)
    if a.value==0 throw(DivideError()) end

    u = a.value
    v = a.field.reduction + (BigInt(1)<<a.field.order)
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

#square and multiply method
function ^(a::FieldPoint, b::Integer)
    a = reduce(a)

    result = FieldPoint(1, a.field)
    squaring = a

    while b>BigInt(0)
        if b & BigInt(1) == BigInt(1)
            result *= squaring
        end
        squaring *= squaring
        b >>>= 1
    end

    return result
end

function sqrt(a::FieldPoint)
    return a^(BigInt(1)<<(a.field.order-1))
end

function random(f::Field)
    range = BigInt(0):((BigInt(1)<<f.order)-BigInt(1))
    return FieldPoint(rand(range), f)
end

function iszero(a::FieldPoint)
    return a.value==0
end

#sec2 v2, table 3:
const FIELD163 = Field(163, BigInt(128+64+8+1))
const FIELD233 = Field(233, (BigInt(1)<<74) + BigInt(1))
const FIELD239A = Field(239,  (BigInt(1)<<36) + BigInt(1))
const FIELD239B = Field(239, (BigInt(1)<<158) + BigInt(1))
const FIELD283 = Field(283, (BigInt(1)<<12) + BigInt(128+32+1))
const FIELD409 = Field(409, (BigInt(1)<<87) + BigInt(1))
const FIELD571 = Field(571, (BigInt(1)<<10) + BigInt(32+4+1))
