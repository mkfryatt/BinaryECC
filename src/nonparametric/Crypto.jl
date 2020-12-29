import Base: isvalid
using SHA: sha256

struct ECKeyPair
    d::BigInt
    Q::ECPointAffine
end

struct ECDSASignature
    r::BigInt
    s::BigInt
end

#sec1 v2, 3.2.1
#Elliptic Curve Key Pair Generation Primitive
function generate_keypair(T::CurveDomainParams)
    d = rand(1::T.n)
    Q = T.G*d
    return ECKeyPair(d, Q)
end

#sec1 v2, 3.2.2.1
#Elliptic Curve Public Key Validation Primitive
function isvalid(Q::ECPointAffine, T::CurveDomainParams)
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
function ecdsa_sign(T::CurveDomainParams, U::ECKeyPair, M::String)
    #loops until it chooses ephemeral key pair that results in nonzero r,s
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
function ecdsa_verify(T::CurveDomainParams, Q::ECPointAffine, M::String, sig::ECDSASignature)
    #1
    if sig.r>=T.n || sig.s>=T.n return false end

    #2
    H = sha256(M)

    #3
    e = digest_to_int(H, T.n)

    #4
    s_inv = invert(s, T.n)
    u1 = (e*s_inv) % T.n
    u2 = (r*s_inv) % T.n

    #5
    R = u1*G + u2*Q
    if iszero(R) return false end

    #6
    xbar = convert(BigInt, R.x)

    #7
    v = xbar % T.n

    #8
    return v==r
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
