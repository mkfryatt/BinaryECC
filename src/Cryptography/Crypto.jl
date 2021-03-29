"""
    ECKeyPair{D,R}
Represents an elliptic curve key pair (described in SEC 1, version 2, 3.2) with
 fields ``d`` and ``Q`` (where ``Q = d \\cdot G``, and ``G`` is the generator of
  the curve domain paramters used to generate this key pair).
"""
struct ECKeyPair{D,R}
    d::PFieldPoint
    Q::ECPointAffine{D,R}
end


"""
    ECDSASignature
Represents a signature produced by ECDSA (Elliptic Curve DSA), with the fields
 ``r`` and ``s`` (both integers).
"""
struct ECDSASignature
    r::PFieldPoint
    s::PFieldPoint
end

function ==(ukey::ECKeyPair{D,R}, vkey::ECKeyPair{D,R}) where {D,R}
    return ukey.d==vkey.d && ukey.Q==vkey.Q
end

function ==(sig1::ECDSASignature, sig2::ECDSASignature)
    return sig1.r==sig2.r && sig1.s==sig2.s
end

#sec1 v2, 3.2.1
#Elliptic Curve Key Pair Generation Primitive
"""
    generate_keypair(T::CurveDomainParams{D,R}) where {D,R}
Gnerates a new random ECKeyPair associated with T, as described in SEC 1 (version 2)
 3.2.1.
"""
function generate_keypair(T::CurveDomainParams{D,R}) where {D,R}
    d = random(PFieldPoint, T.n)
    Q = T.G*d
    return ECKeyPair{D,R}(d, Q)
end

#sec1 v2, 3.2.2.1
#Elliptic Curve Public Key Validation Primitive
"""
    isvalid(T::CurveDomainParams{D,R}, Q::ECPointAffine{D,R}) where {D,R}
Returns true if ``Q`` is a valid public key associated with the curve domain
 parameters ``T``, using the procedure in SEC 1 (version 2) 3.2.2.1, and false
 otherwise.
"""
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
"""
    ecdsa_sign(T::CurveDomainParams{D,R}, U::ECKeyPair{D,R}, M::String) where {D,R}
Creates an ECDSASignature using the key pair ``U`` (associated with the curve
 domain parameters ``T``) for the message ``M`` (a string).

This follows the signing  procedure described in SEC 1 (version 2) 4.1.3.
"""
function ecdsa_sign(T::CurveDomainParams{D,R}, U::ECKeyPair{D,R}, M::String) where {D,R}
    #loops until it chooses ephemeral key pair that results in nonzero r,s
    r, s = zero(PFieldPoint,T.n), zero(PFieldPoint,T.n)
    while true
        #1
        ephemeral = generate_keypair(T)

        #2
        r = PFieldPoint(ephemeral.Q.x, T.n)

        #3
        if iszero(r) continue end

        #4
        H = sha256(M)

        #5
        e = from_digest(PFieldPoint, H, T.n)

        #6
        s = (e + r*U.d) / ephemeral.d
        if iszero(s) continue end
        break
    end

    #7
    return ECDSASignature(r,s)
end

#sec1 v2, 4.1.4
#Verifying Operation
"""
    ecdsa_verify(T::CurveDomainParams{D,R}, Q::ECPointAffine{D,R}, sig::ECDSASignature, M::String) where {D,R}
Returns true if ``\\textit{sig}`` is valid signature for message ``M`` and
 public key ``Q`` (associated with curve domain parameters ``T``), following the
 verifying operation described in SEC 1 (version 2) 4.1.4, and false otherwise.
"""
function ecdsa_verify(T::CurveDomainParams{D,R}, Q::ECPointAffine{D,R}, sig::ECDSASignature, M::String) where {D,R}
    #1
    if !isvalid(sig.r) || !isvalid(sig.s) return false end

    #2
    H = sha256(M)

    #3
    e = from_digest(PFieldPoint, H, T.n)

    #4
    s_inv = inv(sig.s)
    u1 = e * s_inv
    u2 = sig.r * s_inv

    #5
    R1 = u1*T.G + u2*Q
    if iszero(R1) return false end

    #6
    #7
    v = PFieldPoint(R1.x, T.n)

    #8
    return v==sig.r
end

#sec1 v2, 3.3.1
#Elliptic Curve Diffie Hellman Primitive
"""
    ecdh_calculate(T::CurveDomainParams{D,R}, dU::PFieldPoint, QV::ECPointAffine{D,R}) where {D,R}
Calculates the shared secret value for entity "U"'s private key
 (``\\textit{dU}``) and entity "V"'s public key (``\\textit{QV}``), which are
 associated with curve domain parameters ``T``.

This follows the procedure described in SEC 1 (version 2) 3.3.1.
"""
function ecdh_calculate(T::CurveDomainParams{D,R}, dU::PFieldPoint, QV::ECPointAffine{D,R}) where {D,R}
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
"""
    ecdh_deployment1(T::CurveDomainParams)
Performs the first stage of the ECDH deployment operation (described in SEC 1,
 version 2, 6.1.2) from the perspective of entity "U".
"""
function ecdh_deployment1(T::CurveDomainParams)
    #1
    return generate_keypair(T)
end
"""
    ecdh_deployment2(T::CurveDomainParams{D,R}, QV::ECPointAffine{D,R}) where {D,R}
Performs the second stage of the ECDH deployment operation (described in SEC 1,
 version 2, 6.1.2) from the perspective of entity "U", using entity "V"'s
 public key (``\\textit{QV}``).
"""
function ecdh_deployment2(T::CurveDomainParams{D,R}, QV::ECPointAffine{D,R}) where {D,R}
    #3
    #QV has been received
    #4
    return isvalid(T, QV)
end

#sec1 v2, 6.1.3
#ECDH Key Agreement Operation
#TODO finish
"""
    ecdh_agreement(T::CurveDomainParams{D,R}, ukey::ECKeyPair{D,R}, QV::ECPointAffine{D,R}) where {D,R}
This performs the ECDH key agreement operation as described in SEC 1 (version 2) 6.1.3.

It is performed from the perspective of entity "U", using their ECKeyPair ``\\textit{ukey}`` and
the public key of entity "V" (``\\textit{QV}``).
"""
function ecdh_agreement(T::CurveDomainParams{D,R}, ukey::ECKeyPair{D,R}, QV::ECPointAffine{D,R}) where {D,R}
    #1
    z = ecdh_calculate(T, ukey.d, QV)
    return z

    #2
    #z = convert(String, z)

    #3
    #K = use kdf to generate keying data of length keydatalen

    #4
    return K
end