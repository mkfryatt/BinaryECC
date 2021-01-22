module BinaryECC

import Base:
    +, -, *, /, ^, ==, ⊻, >>, <<,
    repr, inv, isvalid,
    iszero, isone, zero, one,
    convert, length, copy

using SHA: sha256
using StaticArrays: MVector

include("StaticUInt.jl")
include("Field_StaticUInt.jl")
include("Field_fastreduce.jl")
#include("Field_BigInt.jl")
include("EC.jl")
include("ECAffine.jl")
include("ECLD.jl")
include("ECJacobian.jl")
include("ECProjective.jl")
include("ECMix.jl")
include("CurveDomainParams.jl")
include("Crypto.jl")

export
    FieldPoint,
    EC,
    AbstractECPoint,
    ECPointAffine,
    ECPointLD,
    ECPointJacobian,
    ECPointProjective,
    CurveDomainParams,
    ECMismatchException,
    ECKeyPair,
    ECDSASignature,
    StaticUInt,

    +,
    -,
    *,
    montmul,
    /,
    inv,
    ^,
    ==,

    repr,
    convert,
    isvalid,
    iszero,
    isone,
    zero,
    one,
    random,
    comb_mult,
    noreduce_mult,

    generate_keypair,
    ecdsa_sign,
    ecdsa_verify,
    ecdh_calculate,
    ecdh_deployment1,
    ecdh_deployment2,
    ecdh_agreement,

    SECT163K1,
    SECT163R1,
    SECT163R2,
    SECT233K1,
    SECT233R1,
    SECT283K1,
    SECT283R1,
    SECT409K1,
    SECT409R1,
    SECT571K1,
    SECT571R1,

    FieldPoint113,
    FieldPoint131,
    FieldPoint163,
    FieldPoint193,
    FieldPoint233,
    FieldPoint239,
    FieldPoint283,
    FieldPoint409,
    FieldPoint571

end
