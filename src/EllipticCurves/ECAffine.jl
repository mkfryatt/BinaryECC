"""
    ECPointAffine{D,R,T} <: AbstractECPoint{D,R,T}
Represents a point on an elliptic curve over the field represented by D and R.
Contains fields ``x``, ``y``, and the elliptic field ("ec") that it is on.

``E: y^2 +  xy = x^3 + ax^2 + b``
"""
struct ECPointAffine{D,R,T} <: AbstractECPoint{D,R,T}
    x::BFieldPoint{D,R,T}
    y::BFieldPoint{D,R,T}
    ec::EC{D,R,T}
end

"""
    ECPointAffine(s::String, ec::EC{D,R,T}) where {D,R,T}
Convert a hex string to a point on the given elliptic curve
using the procedure in SEC 2 (version 2), section 2.3.4.
"""
function ECPointAffine(s::String, ec::EC{D,R,T})::ECPointAffine{D,R,T} where {D,R,T}
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

    x = BFieldPoint{D,R,T}(s[3:4+2*floor(Int16, D / 8)])
    y = BFieldPoint{D,R,T}(s[5+2*floor(Int16, D / 8):6+4*floor(Int16, D / 8)])

    return ECPointAffine(x, y, ec)
end

"""
    repr(p::ECPointAffine)
Returns a string representation of an elliptic curve point, "``(x, y)``".
"""
function repr(p::ECPointAffine)::String
    return "("*repr(p.x)*", "*repr(p.y)*")"
end

"""
    ==(p1::ECPointAffine{D,R,T}, p2::ECPointAffine{D,R,T}) where {D,R,T}
Two points are equal iff they have the same ``x`` and ``y`` coordinate,
 and are on the same elliptic curve.
"""
function ==(p1::ECPointAffine{D,R,T}, p2::ECPointAffine{D,R,T})::Bool where {D,R,T}
    return p1.x==p2.x && p1.y==p2.y && iszero(p1)==iszero(p2) && p1.ec==p2.ec
end

"""
    +(p1::ECPointAffine{D,R,T}, p2::ECPointAffine{D,R,T}) where {D,R,T}
Returns ``p_1+p_2``.

If the points are not on the same curve, this will throw an ECMismatchException.
"""
function +(p1::ECPointAffine{D,R,T}, p2::ECPointAffine{D,R,T})::ECPointAffine{D,R,T} where {D,R,T}
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
    if p1==-p2 return zero(ECPointAffine{D,R,T}, p1.ec) end
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
    -(p::ECPointAffine{D,R,T}) where {D,R,T}
Returns additive inverse of the given point, ``-p``.
"""
function -(p::ECPointAffine{D,R,T})::ECPointAffine{D,R,T} where {D,R,T}
    if iszero(p) return p end
    return ECPointAffine(p.x, p.x+p.y, p.ec)
end

function double(p::ECPointAffine{D,R,T}) where {D,R,T}
    return double_threaded(p)
end

function double_standard(p::ECPointAffine{D,R,T}) where {D,R,T}
    if iszero(p) return p end
    if p==-p return zero(ECPointAffine{D,R,T}, p.ec) end

    #Adds: 5
    #Mults: 2
    #Sqrs: 2
    #Invs: 1
    lambda = p.x + (p.y / p.x)
    x_new = lambda^2 + lambda + p.ec.a
    y_new = p.x^2 + lambda*x_new + x_new
    return ECPointAffine(x_new, y_new, p.ec)
end

function double_threaded(p::ECPointAffine{D,R,T}) where {D,R,T}
    if iszero(p) return p end
    if p==-p return zero(ECPointAffine{D,R,T}, p.ec) end

    lambda = p.x + (p.y / p.x)
    x_new_task = Threads.@spawn $lambda^2 + $lambda + $(p.ec.a)
    y_new = (p.x)^2
    x_new::BFieldPoint{D,R,T} = fetch(x_new_task)
    y_new += lambda*x_new + x_new
    return ECPointAffine(x_new, y_new, p.ec)
end

"""
    mult_mont_affine(p::ECPointAffine{D,R,T}, n::Integer) where {D,R,T}
Returns ``p \\cdot n``.

More resistant to timing attacks than the standard double and add algorithm.

Fast Multiplication on Elliptic Curves over ``GF(2^m)`` without Precomputation,
Algorithm 2A: Montgomery Scalar Multiplication.
"""
function mult_mont_affine(p::ECPointAffine{D,R,T}, n::I)::ECPointAffine{D,R,T} where {D,R,T} where I<:Integer
    if n==0 || iszero(p) return p end
    if n<0 return mont_pow_ladder(-p, -n) end

    b = p.ec.b
    x = p.x

    x1 = x
    x2 = x^2
    x2 += b/x2
    v1, v2 = I(1), I(2)
    for i in (bits(n)-2):-1:0
        if x1==x2
            #ie p1==-p2
            if i==0
                return zero(ECPointAffine{D,R,T}, p.ec)
            else
                #think of a better way to recover here
                #because this repeats a lot of the arithmetic that has already been done
                return find_point(x1, x2, p) + mont_pow_ladder(p, n-v1)
            end
        end
        t = x1 / (x1+x2)
        if (n>>>i)&1==1
            x1 = x + t^2 + t
            x2 = x2^2
            x2 += b/x2
            v1, v2 = v1+v2, 2*v2

        else
            x1 = x1^2
            x1 += b/x1
            x2 = x + t^2 + t
            v1, v2 = 2*v1, v1+v2
        end
    end
    return find_point(x1, x2, p)
end

#needed for montgomery's powering ladder
#finds y1 given that (x1, y1) + p == (x2, y2)
function find_point(x1::BFieldPoint{D,R,T}, x2::BFieldPoint{D,R,T}, p::ECPointAffine{D,R,T})::ECPointAffine{D,R,T} where {D,R,T}
    r1 = x1+p.x
    r2 = x2+p.x
    y1 = r1*r2 + p.x^2 + p.y
    y1 *= r1/p.x
    y1 += p.y
    return ECPointAffine{D,R,T}(x1, y1, p.ec)
end

"""
    isvalid(p::ECPointAffine)
Returns true if ``p`` is a point on the elliptic curve that it is associated with.
"""
function isvalid(p::ECPointAffine)::Bool
    return iszero(p) || (p.y^2 + p.x*p.y == p.x^3 + p.ec.a*p.x^2 + p.ec.b)
end

"""
    iszero(p::ECPointAffine)
Returns true if ``p = \\mathcal{O}``, i.e it is the point at infinity.
"""
function iszero(p::ECPointAffine)::Bool
    return iszero(p.x) && iszero(p.y)
end

"""
    zero(::Type{ECPointAffine, ec::EC{D,R,T}) where {D,R,T}
Returns an object representing the point at infinity on the given curve.
"""
function zero(::Type{ECPointAffine}, ec::EC{D,R,T})::ECPointAffine{D,R,T} where {D,R,T}
    return ECPointAffine{D,R,T}(BFieldPoint{D,R,T}(0), BFieldPoint{D,R,T}(0), ec)
end

function zero(::Type{ECPointAffine{D,R,T}}, ec::EC{D,R,T})::ECPointAffine{D,R,T} where {D,R,T}
    return ECPointAffine{D,R,T}(BFieldPoint{D,R,T}(0), BFieldPoint{D,R,T}(0), ec)
end
