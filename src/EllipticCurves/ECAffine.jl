"""
    ECPointAffine{D,R} <: AbstractECPoint{D,R}
Represents a point on an elliptic curve over the field represented by D and R.
Contains fields ``x``, ``y``, and the elliptic field ("ec") that it is on.

``E: y^2 +  xy = x^3 + ax^2 + b``
"""
struct ECPointAffine{D,R} <: AbstractECPoint{D,R}
    x::BFieldPoint{D,R}
    y::BFieldPoint{D,R}
    ec::EC{D,R}
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

    x = BFieldPoint{D,R}(s[3:4+2*floor(Int16, D / 8)])
    y = BFieldPoint{D,R}(s[5+2*floor(Int16, D / 8):6+4*floor(Int16, D / 8)])

    return ECPointAffine(x, y, ec)
end

"""
    repr(p::ECPointAffine)
Returns a string representation of an elliptic curve point, "``(x, y)``".
"""
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
    if p1==-p2 return zero(ECPointAffine{D,R}, p1.ec) end
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
    if p==-p return zero(ECPointAffine{D,R}, p.ec) end

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
    *(p::ECPointAffine{D,R}, n::Integer) where {D,R}
Returns the result of the scalar multiplication ``p \\cdot n``, using a double and add method.
"""
function *(p::ECPointAffine{D,R}, n::Integer) where {D,R}
    if n<0 return (-p)*(-n) end
    if iszero(p) return p end
    if n==0 return zero(ECPointAffine{D,R}, p.ec) end
    if n==1 return p end

    result = zero(ECPointAffine{D,R}, p.ec)
    while n>0
        if n&1==1
            result += p
        end
        p = double(p)
        n >>= 1
    end
    return result
end

"""
    mult_mont_affine(p::ECPointAffine{D,R}, n::Integer) where {D,R}
Returns ``p \\cdot n``.

More resistant to timing attacks than the standard double and add algorithm.

Fast Multiplication on Elliptic Curves over ``GF(2^m)`` without Precomputation,
Algorithm 2A: Montgomery Scalar Multiplication.
"""
function mult_mont_affine(p::ECPointAffine{D,R}, n::T) where {D,R} where T<:Integer
    if n==0 || iszero(p) return p end
    if n<0 return mont_pow_ladder(-p, -n) end

    b = p.ec.b
    x = p.x

    x1 = x
    x2 = x^2
    x2 += b/x2
    v1, v2 = T(1), T(2)
    for i in (bits(n)-2):-1:0
        if x1==x2
            #ie p1==-p2
            if i==0
                return zero(ECPointAffine{D,R}, p.ec)
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
function find_point(x1::BFieldPoint{D,R}, x2::BFieldPoint{D,R}, p::ECPointAffine{D,R}) where {D,R}
    r1 = x1+p.x
    r2 = x2+p.x
    y1 = r1*r2 + p.x^2 + p.y
    y1 *= r1/p.x
    y1 += p.y
    return ECPointAffine{D,R}(x1, y1, p.ec)
end

"""
    mult_naf(p::ECPointAffine{D,R}, n::Integer) where {D,R}
Returns ``p \\cdot n``.

Uses the binary NAF multiplication method described in Guide to Elliptic Curve Cryptography,
algorithm 3.31.
"""
function mult_naf(p::ECPointAffine{D,R}, n::Integer) where {D,R}
    (adds, subs, l) = naf(n)
    Q = zero(ECPointAffine{D,R}, p.ec)
    for i in (l-1):-1:0
        Q = double(Q)
        if (adds>>i)&1==1
            Q += p
        elseif (subs>>i)&1==1
            Q -= p
        end
    end
    return Q
end

#precompute array of scalar mults of P
function precompute(P::ECPointAffine{D,R}, n, step) where {D,R}
    precomp::Array{ECPointAffine{D,R},1} = []
    Pn = P*step
    append!(precomp, [P])
    for i in 2:n
        append!(precomp, [precomp[i-1]+Pn])
    end
    return precomp
end

#Guide to ECC, algorithm 3.36
#Window NAF method for point multiplication
function mult_naf(P::ECPointAffine{D,R}, k::Integer, w::Int) where {D,R}
    if w==1 return mult_naf(P, k) end
    naf_k = naf(k, w)

    precomp = precompute(P, 1<<(w-2), 2)

    Q = zero(ECPointAffine{D,R}, P.ec)
    for i in length(naf_k):-1:1
        Q = double(Q)
        if naf_k[i]>0
            Q += precomp[naf_k[i]>>1+1]
        elseif naf_k[i]<0
            Q -= precomp[(-naf_k[i])>>1+1]
        end
    end
    return Q
end


#Guide to ECC, algorithm 3.38
#Window NAF method for point multiplication
function mult_naf_window(P::ECPointAffine{D,R}, k::Integer, w::Int) where {D,R}
    (adds, subs, l) = naf(k)

    #make this size  more accurate
    precomp = precompute(P, 1<<(w-1), 2)

    Q = zero(ECPointAffine{D,R}, P.ec)
    i = l-1
    while i>=0
        if (adds>>i)&1==0 && (subs>>i)&1==0
            t, u = 1, 0
        else
            t = w
            while t>0
                if (adds>>(i-t+1))&1==1 || (subs>>(i-t+1))&1==1 break end
                t -= 1
            end
            u = (adds>>>(i-t+1)) & (1<<t -1)
            u -= (subs>>>(i-t+1)) & (1<<t -1)
        end

        for j in 1:t Q = double(Q) end

        if u>0 Q += precomp[u÷2+1]
        elseif u<0 Q -= precomp[(-u)÷2+1]
        end

        i -= t
    end
    return Q
end

#windowed scalar mult, left to right
function mult_window(P::ECPointAffine{D,R}, k::Integer, w::Int) where {D,R}
    l = bits(k)
    precomp = precompute(P, 1<<w -1, 1)
    Q = zero(ECPointAffine{D,R}, P.ec)
    for i in (l-(l%w)):-w:0
        for j in 1:w Q = double(Q) end
        n = (k>>i)&(1<<w -1)
        if n>0 Q += precomp[n] end
    end
    return Q
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

"""
    zero(::Type{ECPointAffine{D,R}}, ec::EC{D,R}) where {D,R}
Returns an object representing the point at infinity on the given curve.
"""
function zero(::Type{ECPointAffine{D,R}}, ec::EC{D,R}) where {D,R}
    return ECPointAffine{D,R}(BFieldPoint{D,R}(0), BFieldPoint{D,R}(0), ec)
end

"""
    zero(::Type{ECPointAffine, ec::EC{D,R}) where {D,R}
Returns an object representing the point at infinity on the given curve.
"""
function zero(::Type{ECPointAffine}, ec::EC{D,R}) where {D,R}
    return ECPointAffine{D,R}(BFieldPoint{D,R}(0), BFieldPoint{D,R}(0), ec)
end
