# Cryptographic Primitives

## Curve Domain Parameters
```@docs
CurveDomainParams{D,R}

isvalid(T::CurveDomainParams{D,R}, t::Int) where {D,R}
```

## Elliptic Curve Key Pairs
```@docs
ECKeyPair{D,R}

generate_keypair(T::CurveDomainParams{D,R}) where {D,R}

isvalid(T::CurveDomainParams{D,R}, Q::ECPointAffine{D,R}) where {D,R}
```

## ECDSA
```@docs
ECDSASignature

ecdsa_sign(T::CurveDomainParams{D,R}, U::ECKeyPair{D,R}, M::String) where {D,R}

ecdsa_verify(T::CurveDomainParams{D,R}, Q::ECPointAffine{D,R}, sig::ECDSASignature, M::String) where {D,R}
```

## ECDH
```@docs
ecdh_calculate(T::CurveDomainParams{D,R}, dU::BigInt, QV::ECPointAffine{D,R}) where {D,R}

ecdh_deployment1(T::CurveDomainParams)

ecdh_deployment2(T::CurveDomainParams{D,R}, QV::ECPointAffine{D,R}) where {D,R}

ecdh_agreement(T::CurveDomainParams{D,R}, ukey::ECKeyPair{D,R}, QV::ECPointAffine{D,R}) where {D,R}
```
