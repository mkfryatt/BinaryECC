"""
    ECPointAffine{D,R} <: AbstractECPoint
Represents a point on an elliptic curve over the field represented by D and R.
Contains fields ``x``, ``y``, and the elliptic field ("ec") that it is on.

``E: y^2 +  xy = x^3 + ax^2 + b``
"""
struct ECPointAffine{D,R} <: AbstractECPoint
    x::FieldPoint{D,R}
    y::FieldPoint{D,R}
    ec::EC{D,R}
end

"""
    ECPointAffine(ec::EC{D,R}) where {D,R}
Returns an object representing the point at infinity on the given curve.
"""
function ECPointAffine(ec::EC{D,R}) where {D,R}
    return ECPointAffine{D,R}(FieldPoint{D,R}(0), FieldPoint{D,R}(0), ec)
end

"""
    ECPointAffine(s::String, ec::EC{D,R}) where {D,R}
Convert a hex string to a point on the given elliptic curve
using the procedure in SEC 2 (version 2), section 2.3.4.
"""
function ECPointAffine(s::String, ec::EC{D,R}) where {D,R}
    s = replace(s, " " => "")

    #point is id
    if s=="00" return ECPointAffine(ec) end

    #input string is not of the uncompressed format specified by sec1v2
    if length(s) != 4*ceil(Int16, D / 8)+2
        throw(ArgumentError("Octet string is of the wrong length for this curve."))
    end
    if s[1:2]!="04"
        throw(ArgumentError("Octet string must start with '04'."))
    end

    #TODO handle compressed format

    x = FieldPoint{D,R}(s[3:4+2*floor(Int16, D / 8)])
    y = FieldPoint{D,R}(s[5+2*floor(Int16, D / 8):6+4*floor(Int16, D / 8)])

    return ECPointAffine(x, y, ec)
end

function repr(p::ECPointAffine)
    return "("*repr(p.x)*", "*repr(p.y)*")"
end

"""
    ==(p1::ECPointAffine{D,R}, p2::ECPointAffine{D,R}) where {D,R}
Two points are equal iff they have the same ``x`` and ``y`` coordinate,
 and are on the same elliptic curve.
"""
function ==(p1::ECPointAffine{D,R}, p2::ECPointAffine{D,R}) where {D,R}
    return p1.x==p2.x && p1.y==p2.y && iszero(p1)==iszero(p2) && p1.ec==p2.ec
end

"""
    +(p1::ECPointAffine{D,R}, p2::ECPointAffine{D,R}) where {D,R}
Returns ``p_1+p_2``.

If the points are not on the same curve, this will throw an ECMismatchException.
"""
function +(p1::ECPointAffine{D,R}, p2::ECPointAffine{D,R}) where {D,R}
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
    if p1==-p2 return ECPointAffine(p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 8
    #Mults: 2
    #Sqrs: 1
    #Invs: 1
    x1x2 = p1.x + p2.x
    lambda = (p1.y+p2.y) / x1x2
    x3 = lambda^2 + lambda + x1x2 + p1.ec.a
    y3 = lambda*(p1.x+x3) + x3 + p1.y
    return ECPointAffine(x3, y3, p1.ec)
end

"""
    -(p::ECPointAffine) where {D,R}
Returns additive inverse of the given point, ``-p``.
"""
function -(p::ECPointAffine) where {D,R}
    if iszero(p) return p end
    return ECPointAffine(p.x, p.x+p.y, p.ec)
end

function double(p::ECPointAffine{D,R}) where {D,R}
    if iszero(p) return p end
    if p==-p return ECPointAffine(p.ec) end

    #Adds: 5
    #Mults: 2
    #Sqrs: 2
    #Invs: 1
    lambda = p.x + (p.y / p.x)
    x_new = lambda^2 + lambda + p.ec.a
    y_new = p.x^2 + lambda*x_new + x_new
    return ECPointAffine(x_new, y_new, p.ec)
end

"""
    *(p::ECPointAffine, n::Integer) where {D,R}
Returns the result of the scalar multiplication ``p \\cdot n``, using a double and add method.
"""
function *(p::ECPointAffine, n::Integer) where {D,R}
    if n<0 return (-p)*(-n) end
    if iszero(p) return p end
    if n==0 return ECPointAffine(p.ec) end
    if n==1 return p end

    result = ECPointAffine(p.ec)
    doubling = p
    while n>0
        if n&1==1
            result += doubling
        end
        doubling = double(doubling)
        n >>= 1
    end
    return result
end

"""
    isvalid(p::ECPointAffine)
Returns true if ``p`` is a point on the elliptic curve that it is associated with.
"""
function isvalid(p::ECPointAffine)
    return iszero(p) || (p.y^2 + p.x*p.y == p.x^3 + p.ec.a*p.x^2 + p.ec.b)
end

"""
    iszero(p::ECPointAffine)
Returns true if ``p = \\mathcal{O}``, i.e it is the point at infinity.
"""
function iszero(p::ECPointAffine)
    return iszero(p.x) && iszero(p.y)
end
