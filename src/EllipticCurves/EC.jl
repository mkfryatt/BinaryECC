"""
    ECMismatchException <: Exception
Indicates that an operation has been attempted
on several points that are not on the same curve.
"""
struct ECMismatchException <: Exception end


"""
    AbstractECPoint
Abstract type for points on an elliptic curve.
"""
abstract type AbstractECPoint end

"""
    EC{D,R}
Represents a non-supersingular elliptic curve over the
field given by D and R.

Contains fields ``a`` and ``b``, where:

``y^2 + xy = x^3 + ax^2 + b``
"""
struct EC{D,R}
    a::BFieldPoint{D,R}
    b::BFieldPoint{D,R}
end

"""
    repr(ec::EC)
Returns a string representation of an elliptic curve equation,
 "``y^2 + xy = x^3 + ax^2 + b``".
"""
function repr(ec::EC)
    return "E: y² + xy = x³ + "*repr(ec.a.value)*"x² + "*repr(ec.b.value)
end

"""
    ==(ec1::EC{D,R}, ec2::EC{D,R}) where {D,R}
Two elliptic curves are equal if they have the
same ``a`` and ``b`` values, and defined over the same field.
"""
function ==(ec1::EC{D,R}, ec2::EC{D,R}) where {D,R}
    return ec1.a==ec2.a && ec1.b==ec2.b
end

"""
    -(p1::AbstractECPoint, p2::AbstractECPoint)
Returns ``p_1-p_2``.
"""
function -(p1::AbstractECPoint, p2::AbstractECPoint)
    return p1 + (-p2)
end

function *(n::Integer, p::AbstractECPoint)
    return p*n
end

function *(n::PFieldPoint, p::AbstractECPoint)
    return p*n.value
end

function *(p::AbstractECPoint, n::PFieldPoint)
    return p*n.value
end

"""
    mult_mont_general(p::AbstractECPoint, n::Integer)
Performs ``p \\cdot n`` with a fixed sequence of curve and field operations.
More resistant to timing attacks than the standard double and add algorithm.
"""
function mult_mont_general(p::AbstractECPoint, n::Integer)
    R0 = p
    R1 = double(p)
    for i in (bits(n)-2):-1:0
        bit = (n>>>i)&1
        if bit==0
            R1, R0 = R0+R1, double(R0)
        else
            R0, R1 = R0+R1, double(R1)
        end
    end
    return R0
end

#number of bits in the binary representation of this number
function bits(a::T) where T<:Integer
    return sizeof(T)*8 - leading_zeros(a)
end
function bits(a::BigInt)
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

#Guide to ECC, algorithm 3.30
#computes the non adjacent form of a positive integer
function naf(k::T) where T<:Integer
    i = 0
    adds = T(0)
    subs = T(0)
    while k>=1
        if k%2==1
            ki = 2 - (k%4)
            k -= ki
            if ki==1 adds ⊻= T(1)<<i
            else subs ⊻= T(1)<<i
            end
        end
        k = k ÷ 2
        i += 1
    end

    return (adds, subs, i)
end

#Guide to ECC, algorithm 3.35
#Computing the width-w NAF of a positive integer
function naf(k::Integer, w::Int)
    if k<=0 throw(ArgumentError("Expected a positive value.")) end
    i = 1
    windowsize = 1<<w
    representation = zeros(MVector{bits(k)+1,Int8})
    while k>0
        if k%2==1
            digit = k % windowsize
            if digit>=(1<<(w-1)) digit -= windowsize end
            representation[i] = digit
            k -= digit
        end
        k >>>= 1
        i += 1
    end
    return representation
end
