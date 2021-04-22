"""
    ECPointJacobian{B} <: AbstractECPoint{B}
Represents a point on an elliptic curve over the field represented by D and R.
Contains fields ``x``, ``y``, ``z``, and the elliptic field ("ec") that it is on.

``E: y^2 +  xyz = x^3 + ax^2z^2 + bz^6``

Each point ``(x, y)`` on the curve is represented by a set of equivalent Jacobian points,
``\\{(\\lambda^2 x, \\lambda^3 y, \\lambda) : \\lambda \\in K^* \\}``
(where ``K^*`` is the binary field that the curve is based on).
"""
struct ECPointJacobian{B} <: AbstractECPoint{B}
    x::B
    y::B
    z::B
    ec::EC{B}
end

"""
    repr(p::ECPointJacobian)
Returns a string representation of an elliptic curve point, "``(x, y, z)``".
"""
function repr(p::ECPointJacobian)::String
    return "("*repr(p.x)*", "*repr(p.y)*", "*repr(p.z)*")"
end

function ==(p1::ECPointJacobian{B}, p2::ECPointJacobian{B})::Bool where B
    z1_2 = square(p1.z)
    z2_2 = square(p2.z)
    return iszero(p1)==iszero(p2) && p1.ec==p2.ec && p1.x*z2_2==p2.x*z1_2 && p1.y*p2.z*z2_2==p2.y*p1.z*z1_2
end

function +(p1::ECPointJacobian{BF}, p2::ECPointJacobian{BF})::ECPointJacobian{BF} where BF
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
    if p1==-p2 return zero(ECPointJacobian{BF}, p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 8
    #Mults: 20
    #Sqrs: 6
    #Invs: 0

    z1_2 = square(p1.z)
    z2_2 = square(p2.z)
    A = p1.y*p2.z*z2_2 + p2.y*p1.z*z1_2
    C = p1.x*z2_2 + p2.x*z1_2
    B = p1.z * p2.z * C

    z3 = B * z1_2

    x3 = A*(A+B) + C^3 + p1.ec.a*square(B)
    x3 *= square(z1_2)

    z3_2 = square(z3)
    t = x3*z1_2
    y3 = A*(p1.x*z3_2 + t)
    y3 += C*p2.z*(t*p1.z + p1.y*z3_2)

    return ECPointJacobian{BF}(x3, y3, z3, p1.ec)
end

function -(p::ECPointJacobian{B})::ECPointJacobian{B} where B
    if iszero(p) return p end
    return ECPointJacobian{B}(p.x, p.x+p.y, p.z, p.ec)
end

function double(p::ECPointJacobian{BF})::ECPointJacobian{BF} where BF
    if iszero(p) return p end
    if p==-p return zero(ECPointJacobian{BF}, p.ec) end

    #Adds: 5
    #Mults: 10
    #Sqrs: 6
    #Invs: 0

    z_2 = square(p.z)
    A = square(p.x) + p.y*p.z
    B = p.x*z_2

    z_4 = square(z_2)

    z_new = B*z_4

    x_new = A*(A+B)
    x_new += p.ec.a * square(B)
    x_new *= square(z_4)

    y_new = B * square(p.x*z_new)
    y_new += (A+B)*x_new*z_4

    return ECPointJacobian{BF}(x_new, y_new, z_new, p.ec)
end

function *(p::ECPointJacobian{B}, n::Integer)::ECPointJacobian{B} where B
    if iszero(p) return p end
    if n==0 return zero(ECPointJacobian{B}, p.ec) end
    if n==1 return p end

    result = zero(ECPointJacobian{B}, p.ec)
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

function isvalid(p::ECPointJacobian)::Bool
    return iszero(p) || (square(p.y) + p.x*p.y*p.z == p.x^3 + p.ec.a*square(p.x*p.z) + p.ec.b*p.z^6)
end

"""
    iszero(p::ECPointJacobian)
Returns true if ``p = \\mathcal{O}``, i.e it is the point at infinity.
"""
function iszero(p::ECPointJacobian)::Bool
    return iszero(p.z)
end

"""
    zero(::Type{ECPointJacobian}, ec::EC{B}) where B
Returns an object representing the point at infinity on the given curve.
"""
function zero(::Type{ECPointJacobian}, ec::EC{B})::ECPointJacobian{B} where B
    return ECPointJacobian{B}(B(0), B(0), B(0), ec)
end
function zero(::Type{ECPointJacobian{B}}, ec::EC{B})::ECPointJacobian{B} where B
    return ECPointJacobian{B}(B(0), B(0), B(0), ec)
end
