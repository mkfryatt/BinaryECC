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
Represents an elliptic curve over the
field given by D and R.

Contains fields ``a`` and ``b``, where:

``y^2 + xy = x^3 + ax^2 + b``
"""
struct EC{D,R}
    a::FieldPoint{D,R}
    b::FieldPoint{D,R}
end

"""
    repr(ec::EC)
Returns a string representation of an elliptic curve equation, "``y^2 + xy = x^3 + ax^2 + b``".
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

"""
    montmul(p::AbstractECPoint, n::Integer)
Performs ``p \\cdot n`` with a fixed sequence of curve and field operations.
More resistant to timing attacks than the standard "double and add" algorithm.
"""
function montmul(p::AbstractECPoint, n::Integer)
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
    i = 0
    while a > (T(1)<<i)
        i += 1
    end
    if a == (T(1)<<i)
        return i+1
    else
        return i
    end
end
