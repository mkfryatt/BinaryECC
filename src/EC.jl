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
