"""
    ECMismatchException <: Exception
Indicates that an operation has been attempted
on several points that are not on the same curve.
"""
struct ECMismatchException <: Exception end


"""
    AbstractECPoint{D,R,T}
Abstract type for points on an elliptic curve.
"""
abstract type AbstractECPoint{D,R,T} end

"""
    EC{D,R,T}
Represents a non-supersingular elliptic curve over the
field given by D and R.

Contains fields ``a`` and ``b``, where:

``y^2 + xy = x^3 + ax^2 + b``
"""
struct EC{D,R,T}
    a::BFieldPoint{D,R,T}
    b::BFieldPoint{D,R,T}
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
    ==(ec1::EC{D,R,T}, ec2::EC{D,R,T}) where {D,R,T}
Two elliptic curves are equal if they have the
same ``a`` and ``b`` values, and defined over the same field.
"""
function ==(ec1::EC{D,R,T}, ec2::EC{D,R,T}) where {D,R,T}
    return ec1.a==ec2.a && ec1.b==ec2.b
end

"""
    -(p1::AbstractECPoint, p2::AbstractECPoint)
Returns ``p_1-p_2``.
"""
function -(p1::AbstractECPoint{D,R,T}, p2::AbstractECPoint{D,R,T}) where {D,R,T}
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

"""
    *(p::AbstractECPoint{D,R,T}, n::Integer) where {D,R,T}
Returns the result of the scalar multiplication ``p \\cdot n``, using a double and add method.
"""
function *(P::AbstractECPoint{D,R,T}, k::Integer)::AbstractECPoint{D,R,T} where {D,R,T}
    return mult_window(P, k, 4)
end

function mult_standard(P::AbstractECPoint{D,R,T}, k::Integer)::AbstractECPoint{D,R,T} where {D,R,T}
    if k<0 return (-P)*(-k) end
    if iszero(P) return P end

    Q = zero(typeof(P), P.ec)
    while k>0
        if k&1==1
            Q += P
        end
        P = double(P)
        k >>= 1
    end
    return Q
end

function mult_threaded(P::AbstractECPoint{D,R,T}, k::Integer)::AbstractECPoint{D,R,T} where {D,R,T}
    if k<0 return (-P)*(-k) end
    if iszero(P) return P end

    len = ceil(Int, D/2)
    k1 = k & (1<<len -1)
    k2 = k & (1<<len -1)<<len

    Q1 = Threads.@spawn $P*$k1
    Q2 = Threads.@spawn $P*$k2
    return fetch(Q1)+fetch(Q2)
end

#windowed scalar mult, left to right
function mult_window(P::AbstractECPoint{D,R,T}, k::Integer, w::Int=1)::AbstractECPoint{D,R,T} where {D,R,T}
    if w==1 return mult_standard(P, k) end
    if k<0 return (-P)*(-k) end
    if iszero(P) return P end
    if k==0 return zero(typeof(P), P.ec) end
    if k==1 return P end

    l = bits(k)
    precomp = precompute(P, 1<<w -1, 1)
    Q = zero(typeof(P), P.ec)
    for i in (l-(l%w)):-w:0
        for j in 1:w Q = double(Q) end
        n = (k>>i)&(1<<w -1)
        if n>0 Q += precomp[n] end
    end
    return Q
end

"""
    mult_naf(p::AbstractECPoint{D,R,T}, n::Integer) where {D,R,T}
Returns ``p \\cdot n``.

Uses the binary NAF multiplication method described in Guide to Elliptic Curve Cryptography,
algorithm 3.31.
"""
function mult_bnaf(P::AbstractECPoint{D,R,T}, k::Integer)::AbstractECPoint{D,R,T} where {D,R,T}
    if k<0 return (-P)*(-k) end
    if iszero(P) return P end

    Q = zero(typeof(P), P.ec)
    while k>0
        if k&1==1
            t = 2 - (k%4)
            k -= t
            if t==1 Q+=P
            else Q-=P
            end
        end
        P = double(P)
        k >>= 1
    end
    return Q
end

#Guide to ECC, algorithm 3.38
#Window NAF method for point multiplication
function mult_bnaf_window(P::AbstractECPoint{D,R,T}, k::Integer, w::Int=1)::AbstractECPoint{D,R,T} where {D,R,T}
    if w==1 return mult_bnaf(P, k) end
    if k<0 return (-P)*(-k) end
    if iszero(P) return P end
    if k==0 return zero(typeof(P), P.ec) end
    if k==1 return P end

    (adds, subs, l) = naf(k)

    n = (2^w - (-1)^w) ÷ 3
    precomp = precompute(P, n, 2)

    Q = zero(typeof(P), P.ec)
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

#precompute array of scalar mults of P
function precompute(P::AbstractECPoint{D,R,T}, n::Int, step::Int) where {D,R,T}
    precomp::Array{typeof(P),1} = []
    Pn = P*step
    append!(precomp, [P])
    for i in 2:n
        append!(precomp, [precomp[i-1]+Pn])
    end
    return precomp
end

#Guide to ECC, algorithm 3.36
#Window NAF method for point multiplication
function mult_wnaf(P::AbstractECPoint{D,R,T}, k::Integer, w::Int)::AbstractECPoint{D,R,T} where {D,R,T}
    if w==1 return mult_bnaf(P, k) end

    if k<0 return (-P)*(-k) end
    if iszero(P) return P end
    if k==0 return zero(typeof(P), P.ec) end
    if k==1 return P end

    naf_k = naf(k, w)

    precomp = precompute(P, 1<<(w-2), 2)

    Q = zero(typeof(P), P.ec)
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
