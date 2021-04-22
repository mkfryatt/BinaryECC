"""
    ECPointAffine{B} <: AbstractECPoint{B}
Represents a point on an elliptic curve over the field represented by D and R.
Contains fields ``x``, ``y``, and the elliptic field ("ec") that it is on.

``E: y^2 +  xy = x^3 + ax^2 + b``
"""
struct ECPointAffine{B} <: AbstractECPoint{B}
    x::B
    y::B
    ec::EC{B}
end

"""
    ECPointAffine(s::String, ec::EC{B}) where B
Convert a hex string to a point on the given elliptic curve
using the procedure in SEC 2 (version 2), section 2.3.4.
"""
function ECPointAffine(s::String, ec::EC{BFieldPoint{D,R,T,L}})::ECPointAffine{BFieldPoint{D,R,T,L}} where {D,R,T,L}
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

    x = BFieldPoint{D,R,T,L}(s[3:4+2*floor(Int16, D / 8)])
    y = BFieldPoint{D,R,T,L}(s[5+2*floor(Int16, D / 8):6+4*floor(Int16, D / 8)])

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
    ==(p1::ECPointAffine{B}, p2::ECPointAffine{B}) where B
Two points are equal iff they have the same ``x`` and ``y`` coordinate,
 and are on the same elliptic curve.
"""
function ==(p1::ECPointAffine{B}, p2::ECPointAffine{B})::Bool where B
    return p1.x==p2.x && p1.y==p2.y && iszero(p1)==iszero(p2) && p1.ec==p2.ec
end

"""
    +(p1::ECPointAffine{B}, p2::ECPointAffine{B}) where B
Returns ``p_1+p_2``.

If the points are not on the same curve, this will throw an ECMismatchException.
"""
function +(p1::ECPointAffine{B}, p2::ECPointAffine{B})::ECPointAffine{B} where B
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
    if p1==-p2 return zero(ECPointAffine{B}, p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 8
    #Mults: 2
    #Sqrs: 1
    #Invs: 1
    x1x2 = p1.x + p2.x
    lambda = (p1.y+p2.y) / x1x2
    x3 = square(lambda) + lambda + x1x2 + p1.ec.a
    y3 = lambda*(p1.x+x3) + x3 + p1.y
    return ECPointAffine(x3, y3, p1.ec)
end

"""
    -(p::ECPointAffine{B}) where B
Returns additive inverse of the given point, ``-p``.
"""
function -(p::ECPointAffine{B})::ECPointAffine{B} where B
    if iszero(p) return p end
    return ECPointAffine(p.x, p.x+p.y, p.ec)
end

function double(p::ECPointAffine{B}) where B
    return double_standard(p)
end

function double_standard(p::ECPointAffine{B}) where B
    if iszero(p) return p end
    if p==-p return zero(ECPointAffine{B}, p.ec) end

    #Adds: 5
    #Mults: 2
    #Sqrs: 2
    #Invs: 1
    lambda = p.x + (p.y / p.x)
    x_new = square(lambda) + lambda + p.ec.a
    y_new = square(p.x) + lambda*x_new + x_new
    return ECPointAffine(x_new, y_new, p.ec)
end

function double_threaded(p::ECPointAffine{B}) where B
    if iszero(p) return p end
    if p==-p return zero(ECPointAffine{B}, p.ec) end

    lambda = p.x + (p.y / p.x)
    x_new_task = Threads.@spawn square($lambda) + $lambda + $(p.ec.a)
    y_new = square(p.x)
    x_new::B = fetch(x_new_task)
    y_new += lambda*x_new + x_new
    return ECPointAffine(x_new, y_new, p.ec)
end

"""
    mult_mont_affine(p::ECPointAffine{B}, n::Integer) where B
Returns ``p \\cdot n``.

More resistant to timing attacks than the standard double and add algorithm.

Fast Multiplication on Elliptic Curves over ``GF(2^m)`` without Precomputation,
Algorithm 2A: Montgomery Scalar Multiplication.
"""
function mult_mont_affine(p::ECPointAffine{B}, n::I)::ECPointAffine{B} where B where I<:Integer
    if n==0 || iszero(p) return p end
    if n<0 return mont_pow_ladder(-p, -n) end

    b = p.ec.b
    x = p.x

    x1 = x
    x2 = square(x)
    x2 += b/x2
    v1, v2 = I(1), I(2)
    for i in (bits(n)-2):-1:0
        if x1==x2
            #ie p1==-p2
            if i==0
                return zero(ECPointAffine{B}, p.ec)
            else
                #think of a better way to recover here
                #because this repeats a lot of the arithmetic that has already been done
                return find_point(x1, x2, p) + mont_pow_ladder(p, n-v1)
            end
        end
        t = x1 / (x1+x2)
        if (n>>>i)&1==1
            x1 = x + square(t) + t
            x2 = square(x2)
            x2 += b/x2
            v1, v2 = v1+v2, 2*v2

        else
            x1 = square(x1)
            x1 += b/x1
            x2 = x + square(t) + t
            v1, v2 = 2*v1, v1+v2
        end
    end
    return find_point(x1, x2, p)
end

#needed for montgomery's powering ladder
#finds y1 given that (x1, y1) + p == (x2, y2)
function find_point(x1::B, x2::B, p::ECPointAffine{B})::ECPointAffine{B} where B
    r1 = x1+p.x
    r2 = x2+p.x
    y1 = r1*r2 + square(p.x) + p.y
    y1 *= r1/p.x
    y1 += p.y
    return ECPointAffine{B}(x1, y1, p.ec)
end

"""
    isvalid(p::ECPointAffine)
Returns true if ``p`` is a point on the elliptic curve that it is associated with.
"""
function isvalid(p::ECPointAffine)::Bool
    return iszero(p) || (square(p.y) + p.x*p.y == p.x^3 + p.ec.a*square(p.x) + p.ec.b)
end

"""
    iszero(p::ECPointAffine)
Returns true if ``p = \\mathcal{O}``, i.e it is the point at infinity.
"""
function iszero(p::ECPointAffine)::Bool
    return iszero(p.x) && iszero(p.y)
end

"""
    zero(::Type{ECPointAffine, ec::EC{B}) where B
Returns an object representing the point at infinity on the given curve.
"""
function zero(::Type{ECPointAffine}, ec::EC{B})::ECPointAffine{B} where B
    return ECPointAffine{B}(B(0), B(0), ec)
end

function zero(::Type{ECPointAffine{B}}, ec::EC{B})::ECPointAffine{B} where B
    return ECPointAffine{B}(B(0), B(0), ec)
end
