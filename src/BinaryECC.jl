module BinaryECC

import Base:
    +, -, *, /, ^, ==, ⊻, >>, <<,
    repr, inv, isvalid,
    iszero, isone, zero, one,
    convert, length, copy, sqrt,
    reduce

using SHA: sha256
using StaticArrays: MVector

bitsize(T) = 8*sizeof(T)

include("GaloisFields/StaticUInt.jl")
include("GaloisFields/BField.jl")
include("GaloisFields/BField_fastreduce.jl")
include("GaloisFields/PField.jl")

include("EllipticCurves/EC.jl")
include("EllipticCurves/ECAffine.jl")
include("EllipticCurves/ECLD.jl")
include("EllipticCurves/ECJacobian.jl")
include("EllipticCurves/ECMix.jl")

include("Cryptography/CurveDomainParams.jl")
include("Cryptography/Crypto.jl")

export
    BFieldPoint,
    PFieldPoint,
    PFieldPointMismatchException,
    StaticUInt,

    EC,
    AbstractECPoint,
    ECPointAffine,
    ECPointLD,
    ECPointJacobian,
    ECMismatchException,

    CurveDomainParams,
    ECKeyPair,
    ECDSASignature,

    +,
    -,
    *,
    /,
    inv,
    ^,
    ==,
    sqrt,
    repr,
    convert,
    isvalid,
    iszero,
    isone,
    zero,
    one,
    reduce,
    random,

    #elliptic curve routines
    double,
    mult_standard,
    mult_window,
    mult_threaded,
    mult_threaded_window,
    mult_threaded_ownreduce,
    mult_bnaf,
    mult_wnaf,
    mult_bnaf_window,
    mult_bnaf_window_test,
    mult_mont_general,
    mult_mont_affine,
    naf,

    #binary field routines
    mult_shiftandadd,
    mult_shiftandadd_window,
    mult_threaded,
    mult_ownreduce,
    mult_comb_rtl,
    mult_comb_ltr,
    mult_comb_window,
    square,
    square_standard,
    square_window,

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

    BFieldPoint113,
    BFieldPoint131,
    BFieldPoint163,
    BFieldPoint193,
    BFieldPoint233,
    BFieldPoint239,
    BFieldPoint283,
    BFieldPoint409,
    BFieldPoint571

end
