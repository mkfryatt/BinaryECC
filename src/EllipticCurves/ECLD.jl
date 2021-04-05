"""
    ECPointLD{D,R} <: AbstractECPoint{D,R}
Represents a point on an elliptic curve over the field represented by D and R.
Contains fields ``x``, ``y``, ``z``, and the elliptic field ("ec") that it is on.

``E: y^2 +  xyz = x^3z + ax^2z^2 + bz^4``

Each (affine) point ``(x, y)`` is represented by a set of Lopez-Dahab points,
``\\{(\\lambda x, \\lambda^2 y, \\lambda) : \\lambda \\in K^* \\}``
(where ``K^*`` is the binary field that the curve is based on).
"""
struct ECPointLD{D,R} <: AbstractECPoint{D,R}
    x::BFieldPoint{D,R}
    y::BFieldPoint{D,R}
    z::BFieldPoint{D,R}
    ec::EC{D,R}
end

"""
    repr(p::ECPointLD)
Returns a string representation of an elliptic curve point, "``(x, y, z)``".
"""
function repr(p::ECPointLD)
    return "("*repr(p.x)*", "*repr(p.y)*", "*repr(p.z)*")"
end

function ==(p1::ECPointLD{D,R}, p2::ECPointLD{D,R}) where {D,R}
    return iszero(p1)==iszero(p2) && p1.ec==p2.ec && p1.x*p2.z==p2.x*p1.z && p1.y*p2.z^2==p2.y*p1.z^2
end

function +(p1::ECPointLD{D,R}, p2::ECPointLD{D,R}) where {D,R}
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
    if p1==-p2 return zero(ECPointLD{D,R}, p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 8
    #Mults: 16
    #Sqrs: 4
    #Invs: 0

    A = p1.y*p2.z^2 + p2.y*p1.z^2
    C = p1.x*p2.z + p2.x*p1.z
    E = p2.z*C
    B = p1.z*E

    z3 = B^2

    x3 = A*(A+B) + B*(C^2 + p1.ec.a*B)

    t = x3*p1.z

    y3 = A*(p1.x*z3 + t) + (t*p1.z + p1.y*z3)*E
    y3 *= E

    return ECPointLD{D,R}(x3, y3, z3, p1.ec)
end

function -(p::ECPointLD{D,R}) where {D,R}
    if iszero(p) return p end
    return ECPointLD{D,R}(p.x, p.x+p.y, p.z, p.ec)
end

function double(p::ECPointLD{D,R}) where {D,R}
    if iszero(p) return p end
    if p==-p return zero(ECPointLD{D,R}, p.ec) end

    #Adds: 4
    #Mults: 6
    #Sqrs: 3
    #Invs: 0

    x_2 = p.x^2
    A = x_2 + p.y
    B = p.x*p.z

    z_new = B^2

    AB = A+B

    x_new = A*AB + p.ec.a*z_new

    y_new = x_2^2 * z_new + x_new*B*AB

    return ECPointLD{D,R}(x_new, y_new, z_new, p.ec)
end

function *(p::ECPointLD{D,R}, n::Integer) where {D,R}
    if iszero(p) return p end
    if n==0 return zero(ECPointLD{D,R}, p.ec) end
    if n==1 return p end

    result = zero(ECPointLD{D,R}, p.ec)
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

function isvalid(p::ECPointLD)
    return iszero(p) || (p.y^2 + p.x*p.y*p.z == p.x^3*p.z + p.ec.a*p.x^2*p.z^2 + p.ec.b*p.z^4)
end

"""
    iszero(p::ECPointLD)
Returns true if ``p = \\mathcal{O}``, i.e it is the point at infinity.
"""
function iszero(p::ECPointLD)
    return iszero(p.y) && iszero(p.z)
end

"""
    zero(::Type{ECPointLD{D,R}}, ec::EC{D,R}) where {D,R}
Returns an object representing the point at infinity on the given curve.
"""
function zero(::Type{ECPointLD{D,R}}, ec::EC{D,R}) where {D,R}
    return ECPointLD{D,R}(BFieldPoint{D,R}(0), BFieldPoint{D,R}(0), BFieldPoint{D,R}(0), ec)
end

"""
    zero(::Type{ECPointLD}, ec::EC{D,R}) where {D,R}
Returns an object representing the point at infinity on the given curve.
"""
function zero(::Type{ECPointLD}, ec::EC{D,R}) where {D,R}
    return ECPointLD{D,R}(BFieldPoint{D,R}(0), BFieldPoint{D,R}(0), BFieldPoint{D,R}(0), ec)
end
