# Cryptographic Primitives

## Prime Order Fields
```@docs
PFieldElt
```

## Curve Domain Parameters
```@docs
CurveDomainParams{B}

isvalid(T::CurveDomainParams{B}, t::Int) where B
```

## Elliptic Curve Key Pairs
```@docs
ECKeyPair{B}

generate_keypair(T::CurveDomainParams{B}) where B

isvalid(T::CurveDomainParams{B}, Q::ECPointAffine{B}) where B
```

## ECDSA
```@docs
ECDSASignature

ecdsa_sign(T::CurveDomainParams{B}, U::ECKeyPair{B}, M::String) where B

ecdsa_verify(T::CurveDomainParams{B}, Q::ECPointAffine{B}, sig::ECDSASignature, M::String) where B
```

## ECDH
```@docs
ecdh_calculate(T::CurveDomainParams{B}, dU::PFieldElt, QV::ECPointAffine{B}) where B

ecdh_deployment1(T::CurveDomainParams)

ecdh_deployment2(T::CurveDomainParams{B}, QV::ECPointAffine{B}) where B

ecdh_agreement(T::CurveDomainParams{B}, ukey::ECKeyPair{B}, QV::ECPointAffine{B}) where B
```
