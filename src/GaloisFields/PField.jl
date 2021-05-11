struct PFieldEltMismatchException <: Exception end

"""
    PFieldElt
Represents an element of a prime order field, with the named field `value::BigInt`
holding the element itself, and `p::BigInt` holding the field order.

Supports all standard arithmetic operations, such as `==`, `+`, `-`, `*`, `/`, `inv`, `^`,
`isone`, `iszero`, `one`, `zero`, `isvalid`.
"""
struct PFieldElt
    value::BigInt
    p::BigInt
    PFieldElt(value::Union{BFieldElt,StaticUInt,Integer}, p::Union{BFieldElt,StaticUInt,Integer}) =
        new(convert(BigInt, value)%p, convert(BigInt, p))
end

function repr(x::PFieldElt)
    return repr(x.value)
end

function ==(x::PFieldElt, y::PFieldElt)
    if x.p!=y.p throw(PFieldEltMismatchException()) end
    return x.value==y.value
end

function +(x::PFieldElt, y::PFieldElt)
    if x.p!=y.p throw(PFieldEltMismatchException()) end
    return PFieldElt(x.value+y.value, x.p)
end

function -(x::PFieldElt)
    return PFieldElt(x.p-x.value, x.p)
end

function -(x::PFieldElt, y::PFieldElt)
    if x.p!=y.p throw(PFieldEltMismatchException()) end
    return x + (-y)
end

function *(x::PFieldElt, y::PFieldElt)
    if x.p!=y.p throw(PFieldEltMismatchException()) end
    return PFieldElt(x.value*y.value, x.p)
end

function inv(x::PFieldElt)
    if x.value==0 throw(DivideError()) end
    t, new_t = BigInt(0), BigInt(1)
    r, new_r = x.p, x.value

    while new_r!=0
        quotient = r ÷ new_r
        t, new_t = new_t, t - quotient*new_t
        r, new_r = new_r, r - quotient*new_r
    end

    if t<0
        t = t+x.p
    end

    return PFieldElt(t, x.p)
end

function /(x::PFieldElt, y::PFieldElt)
    if x.p!=y.p throw(PFieldEltMismatchException()) end
    return x * inv(y)
end

function zero(::Type{PFieldElt}, p::Integer)
    return PFieldElt(0, p)
end
function zero(x::PFieldElt)
    return PFieldElt(0, x.p)
end

function one(::Type{PFieldElt}, p::Integer)
    return PFieldElt(1, p)
end
function one(x::PFieldElt)
    return PFieldElt(1, x.p)
end

function iszero(x::PFieldElt)
    return x.value==0
end

function isone(x::PFieldElt)
    return x.value==1
end

#if log2(n)>=8hashlen, returns the digest as a BigInt
#otherwise, returns the leftmost log2(n) bits as a BigInt
function from_digest(::Type{PFieldElt}, H, n)
    hashlen = 32
    len = ceil(Int, log2(n))
    e = BigInt(0)

    if len < 8*hashlen
        bytes = len ÷ 8
        extra_bits = len % 8
        for i in 1:bytes
            e += BigInt(H[i])<<(8*(bytes-i)+extra_bits)
        end
        extras = H[bytes+1]
        extras >>>= 8-extra_bits
        extras &= 2^extra_bits -1
        e += extras

    else
        for i in 0:31
            e += BigInt(H[32-i])<<(8*i)
        end
    end
    return PFieldElt(e, n)
end

function isvalid(x::PFieldElt)
    return x.value<x.p
end

function random(::Type{PFieldElt}, p::Integer)
    return PFieldElt(rand(1:p), p)
end

function ^(x::PFieldElt, y::Integer)
    z = one(PFieldElt, x.p)
    squaring = x
    while y>0
        if y%2==1
            z *= squaring
        end
        squaring *= squaring
        y ÷= 2
    end
    return z
end
