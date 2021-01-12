import Base: isvalid
using SHA: sha256

struct ECKeyPair{D,R}
    d::BigInt
    Q::ECPointAffine{D,R}
end

struct ECDSASignature
    r::BigInt
    s::BigInt
end

#sec1 v2, 3.2.1
#Elliptic Curve Key Pair Generation Primitive
function generate_keypair(T::CurveDomainParams{D,R}) where {D,R}
    d = rand(1:T.n)
    Q = T.G*d
    return ECKeyPair{D,R}(d, Q)
end

#sec1 v2, 3.2.2.1
#Elliptic Curve Public Key Validation Primitive
function isvalid(T::CurveDomainParams{D,R}, Q::ECPointAffine{D,R}) where {D,R}
    #1
    if iszero(Q) return false end

    #2 is for curves over prime fields

    #3
    if T.G.ec!=Q.ec return false end
    if !isvalid(Q) return false end

    #4
    if T.h!=1 && !iszero(T.n*Q) return false end
    #if h==1, then nQ==0 is already implied by step 3

    #5
    return true
end

#sec1 v2, 4.1.3
#Signing Operation
function ecdsa_sign(T::CurveDomainParams{D,R}, U::ECKeyPair{D,R}, M::String) where {D,R}
    #loops until it chooses ephemeral key pair that results in nonzero r,s
    r, s = 0, 0
    while true
        #1
        ephemeral = generate_keypair(T)

        #2
        xbar = convert(BigInt, ephemeral.Q.x)

        #3
        r = xbar % T.n
        if r==0 continue end

        #4
        H = sha256(M)

        #5
        e = digest_to_int(H, T.n)

        #6
        s = (invert(ephemeral.d, T.n)*(e+ r*U.d)) % T.n
        if s==0 continue end
        break
    end

    #7
    return ECDSASignature(r,s)
end

#sec1 v2, 4.1.4
#Verifying Operation
function ecdsa_verify(T::CurveDomainParams{D,R}, Q::ECPointAffine{D,R}, sig::ECDSASignature, M::String) where {D,R}
    #1
    if sig.r>=T.n || sig.s>=T.n return false end

    #2
    H = sha256(M)

    #3
    e = digest_to_int(H, T.n)

    #4
    s_inv = invert(sig.s, T.n)
    u1 = (e*s_inv) % T.n
    u2 = (sig.r*s_inv) % T.n

    #5
    R1 = u1*T.G + u2*Q
    if iszero(R1) return false end

    #6
    xbar = convert(BigInt, R1.x)

    #7
    v = xbar % T.n

    #8
    return v==sig.r
end

#sec1 v2, 3.3.1
#Elliptic Curve Diffie Hellman Primitive
function ecdh_calculate(T::CurveDomainParams{D,R}, dU::BigInt, QV::ECPointAffine{D,R}) where {D,R}
    #check that QV is associated with T and valid
    if !isvalid(T, QV) return nothing end

    #1
    P = dU*QV

    #2
    if iszero(P) return nothing end

    #3
    return P.x
end

#sec1 v2, 6.1.2
#ECDH Key Deployment Operation
#from the perspective of U
function ecdh_deployment1(T::CurveDomainParams)
    #1
    return generatekeypair(T)
end
function ecdh_deployment2(T::CurveDomainParams{D,R}, QV::ECPointAffine{D,R}) where {D,R}
    #3
    #QV has been received
    #4
    return isvalid(T, QV)
end

#sec1 v2, 6.1.3
#ECDH Key Agreement Operation
function ecdh_agreement(T::CurveDomainParams{D,R}, ukey::ECKeyPair{D,R}, QV::ECPointAffine{D,R}) where {D,R}
    #1
    z = ecdh_calculate(T, ukey.d, QV)

    #2
    z = convert(String, z)

    #3
    #K = use kdf to generate keying data of length keydatalen

    #4
    return K
end

#if log2(n)>=8hashlen, returns the digest as a BigInt
#otherwise, returns the leftmost log2(n) bits as a BigInt
function digest_to_int(H, n)
    hashlen = 32
    len = ceil(Int, log2(n))
    e = BigInt(0)

    if len < 8*hashlen
        bytes = len ÷ 8
        extra_bits = len % 8
        for i in 1:bytes
            e += BigInt(H[i])<<(8*(bytes-i)+extra_bits)
        end
        extras = H[bytes+1]
        extras >>>= 8-extra_bits
        extras &= 2^extra_bits -1
        e += extras

    else
        for i in 0:31
            e += BigInt(H[32-i])<<(8*i)
        end
    end
    return e
end

#finds t such that xt ≡ 1 (mod p)
function invert(x::Integer, p::Integer)
    if x==0 throw(DivideError()) end
    t, new_t = 0, 1
    r, new_r = p, x

    while new_r!=0
        quotient = r ÷ new_r
        t, new_t = new_t, t - quotient*new_t
        r, new_r = new_r, r - quotient*new_r
    end

    if t<0
        t = t+p
    end

    return t
end
