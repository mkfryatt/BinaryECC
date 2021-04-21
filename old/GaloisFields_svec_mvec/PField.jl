struct PFieldPointMismatchException <: Exception end

struct PFieldPoint
    value::BigInt
    p::BigInt
    PFieldPoint(value::Union{BFieldElt,SVector,MVector,Integer}, p::Union{BFieldElt,SVector,MVector,Integer}) =
        new(convert(BigInt, value)%p, convert(BigInt, p))
end

function ==(x::PFieldPoint, y::PFieldPoint)
    if x.p!=y.p throw(PFieldPointMismatchException()) end
    return x.value==y.value
end

function +(x::PFieldPoint, y::PFieldPoint)
    if x.p!=y.p throw(PFieldPointMismatchException()) end
    return PFieldPoint(x.value+y.value, x.p)
end

function -(x::PFieldPoint)
    return PFieldPoint(x.p-x.value, x.p)
end

function -(x::PFieldPoint, y::PFieldPoint)
    if x.p!=y.p throw(PFieldPointMismatchException()) end
    return x + (-y)
end

function *(x::PFieldPoint, y::PFieldPoint)
    if x.p!=y.p throw(PFieldPointMismatchException()) end
    return PFieldPoint(x.value*y.value, x.p)
end

function inv(x::PFieldPoint)
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

    return PFieldPoint(t, x.p)
end

function /(x::PFieldPoint, y::PFieldPoint)
    if x.p!=y.p throw(PFieldPointMismatchException()) end
    return x * inv(y)
end

function zero(::Type{PFieldPoint}, p::Integer)
    return PFieldPoint(0, p)
end

function one(::Type{PFieldPoint}, p::Integer)
    return PFieldPoint(1, p)
end

function iszero(x::PFieldPoint)
    return x.value==0
end

function isone(x::PFieldPoint)
    return x.value==1
end

#if log2(n)>=8hashlen, returns the digest as a BigInt
#otherwise, returns the leftmost log2(n) bits as a BigInt
function from_digest(::Type{PFieldPoint}, H, n)
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
    return PFieldPoint(e, n)
end

function isvalid(x::PFieldPoint)
    return x.value<x.p
end

function random(::Type{PFieldPoint}, p::Integer)
    return PFieldPoint(rand(1:p), p)
end

function ^(x::PFieldPoint, y::Integer)
    z = one(PFieldPoint, x.p)
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
