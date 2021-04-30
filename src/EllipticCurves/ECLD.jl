"""
    ECPointLD{B} <: AbstractECPoint{B}
Represents a point on an elliptic curve over the binary `B`.
Contains fields `x::B`, `y::B`, `z::B` and `ec::EC{B}`.

``E: y^2 +  xyz = x^3z + ax^2z^2 + bz^4``

Each (affine) point ``(x, y)`` is represented by a set of Lopez-Dahab points,
``\\{(\\lambda x, \\lambda^2 y, \\lambda) : \\lambda \\in K^* \\}``
(where ``K^*`` is the binary field that the curve is based on).
"""
struct ECPointLD{B} <: AbstractECPoint{B}
    x::B
    y::B
    z::B
    ec::EC{B}
end

"""
    repr(p::ECPointLD)
Returns a string representation of an elliptic curve point, "``(x, y, z)``".
"""
function repr(p::ECPointLD)::String
    return "("*repr(p.x)*", "*repr(p.y)*", "*repr(p.z)*")"
end

function ==(p1::ECPointLD{B}, p2::ECPointLD{B})::Bool where B
    return iszero(p1)==iszero(p2) && p1.ec==p2.ec && p1.x*p2.z==p2.x*p1.z && p1.y*square(p2.z)==p2.y*square(p1.z)
end

function +(p1::ECPointLD{BF}, p2::ECPointLD{BF})::ECPointLD{BF} where BF
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
    if p1==-p2 return zero(ECPointLD{BF}, p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 8
    #Mults: 16
    #Sqrs: 4
    #Invs: 0

    A = p1.y*square(p2.z) + p2.y*square(p1.z)
    C = p1.x*p2.z + p2.x*p1.z
    E = p2.z*C
    B = p1.z*E

    z3 = square(B)

    x3 = A*(A+B) + B*(square(C) + p1.ec.a*B)

    t = x3*p1.z

    y3 = A*(p1.x*z3 + t) + (t*p1.z + p1.y*z3)*E
    y3 *= E

    return ECPointLD{BF}(x3, y3, z3, p1.ec)
end

function -(p::ECPointLD{B})::ECPointLD{B} where B
    if iszero(p) return p end
    return ECPointLD{B}(p.x, p.x+p.y, p.z, p.ec)
end

function double(p::ECPointLD{BF})::ECPointLD{BF} where BF
    if iszero(p) return p end
    if p==-p return zero(ECPointLD{BF}, p.ec) end

    #Adds: 4
    #Mults: 6
    #Sqrs: 3
    #Invs: 0

    x_2 = square(p.x)
    A = x_2 + p.y
    B = p.x*p.z

    z_new = square(B)

    AB = A+B

    x_new = A*AB + p.ec.a*z_new

    y_new = square(x_2) * z_new + x_new*B*AB

    return ECPointLD{BF}(x_new, y_new, z_new, p.ec)
end

function *(p::ECPointLD{B}, n::Integer)::ECPointLD{B} where B
    if iszero(p) return p end
    if n==0 return zero(ECPointLD{B}, p.ec) end
    if n==1 return p end

    result = zero(ECPointLD{B}, p.ec)
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

function isvalid(p::ECPointLD)::Bool
    return iszero(p) || (square(p.y) + p.x*p.y*p.z == p.x^3*p.z + p.ec.a*square(p.x*p.z) + p.ec.b*p.z^4)
end

"""
    iszero(p::ECPointLD)
Returns true if ``p = \\mathcal{O}``, i.e it is the point at infinity.
"""
function iszero(p::ECPointLD)::Bool
    return iszero(p.y) && iszero(p.z)
end

"""
    zero(::Type{ECPointLD}, ec::EC{B}) where B
Returns an object representing the point at infinity on the given curve.
"""
function zero(::Type{ECPointLD}, ec::EC{B})::ECPointLD{B} where B
    return ECPointLD{B}(B(0), B(0), B(0), ec)
end

"""
    zero(::Type{ECPointLD{B}}, ec::EC{B}) where B
Returns an object representing the point at infinity on the given curve.
"""
function zero(::Type{ECPointLD{B}}, ec::EC{B})::ECPointLD{B} where B
    return ECPointLD{B}(B(0), B(0), B(0), ec)
end
