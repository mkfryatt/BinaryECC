#T is the tuple of curve domain params
#Defined by sec1v2, 3.1.2
#T = (m, f(x), a, b, G, n, h)
#m: D
#f(x): R + x^m
#a: G.ec.a
#b: G.ec.b
#G: G
#n: n
#h: h

"""
    CurveDomainParams{B}
Represents the elliptic curve domain parameters for elliptic curve groups defined
over binary field `B`,  as described in SEC 1 (version 2), 3.1.2.

It contains three fields:
- `G::ECPointAffine{B}`, a generating point, in affine coordinates
- `n::BigInt`, the order of `G` (i.e. the smallest `n` satisfying ``G \\cdot n = \\mathcal{O}``)
- `h::BigInt`, the cofactor, ``h = \\#E(\\mathbb{F}_{2^m}) / n``

The other elements of the septuple described in 3.1.2 are accessible through
the fields of `G`.

Several standard curves domain parameters (taken from SEC 2, section 3) can be created
by calling the following functions with a word type `T`:
- `SECT163K1(T::Type{U}) where U<:Unsigned`
- `SECT163R1(T::Type{U}) where U<:Unsigned`
- `SECT233K1(T::Type{U}) where U<:Unsigned`
- `SECT233R1(T::Type{U}) where U<:Unsigned`
- `SECT283K1(T::Type{U}) where U<:Unsigned`
- `SECT283R1(T::Type{U}) where U<:Unsigned`
- `SECT409K1(T::Type{U}) where U<:Unsigned`
- `SECT409R1(T::Type{U}) where U<:Unsigned`
- `SECT571K1(T::Type{U}) where U<:Unsigned`
- `SECT571R1(T::Type{U}) where U<:Unsigned`
"""
struct CurveDomainParams{B}
    G::ECPointAffine{B} #generator point
    n::BigInt #order of G, ie nG = O
    h::Int #cofactor
    CurveDomainParams(G::ECPointAffine{B}, n::Integer, h::Integer) where B =
        new{B}(G, convert(BigInt, n), h)

    CurveDomainParams{B}(G::ECPointAffine{B}, n::Integer, h::Integer) where B =
        new{B}(G, convert(BigInt, n), h)
end

#sec1 v2, 3.1.2.2.1
#Elliptic Curve Domain Parameters over F_{2^m} Validation Primitive
"""
    isvalid(T::CurveDomainParams{B}, t::Int) where B
Returns true if the curve domain parameters ``T`` meet the security level ``t``,
 using the procedure in SEC 1 (version 2) 3.1.2.2.1, and false otherwise.

Note: does not currently perform step 6 (checking that ``n`` is prime).
"""
function isvalid(T::CurveDomainParams{B}, t::Int) where B
    #1
    levels = [80, 112, 128, 192, 256]
    i=1
    while levels[i]<t
        i+=1
        if i>5 return false end #t has been set too high
    end
    tprime = levels[i]

    if !(D in [163, 233, 239, 283, 409, 571]) return false end
    if D<=2*t || D>=2*tprime return false end

    #2
    reductions = Dict([
        (163, (Int128(128+64+8+1)),
        (233, (Int128(1)<<74) + Int128(1)),
        (239, Int128(1)<<36) + Int128(1)),
        (483, (Int128(1)<<12) + Int128(128+32+1)),
        (409, (Int128(1)<<87) + Int128(1)),
        (571, (Int128(1)<<10) + Int128(32+4+1))
    ])
    if R!=reductions[D] return false end

    #3
    #check that a,b,x,y are all in the field
    #this should already be true due to checks in the constructor for G

    #4
    if iszero(T.G.ec.b) return false end

    #5
    if !isvalid(T.G) return false end

    #6
    #TODO Check that n is prime

    #7
    if T.h> 2^(t/8) return false end
    if T.h != floor(Int, ((sqrt(2^D)+1)^2/T.n)) return false end

    #8
    if !iszero(T.n*T.G) return false end

    #9
    if T.n*T.h == 2^m return false end
    for b in 1:(200*D)
        if (2^b % T.n)==0 return false end
    end

    return true
end

function random(T::CurveDomainParams{B})::ECPointAffine{B} where {B}
    return mult_memo(T.G, rand(1:T.n))
end

#sec2 v2 curve parameters:
SECT163K1(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "0402FE 13C0537B BC11ACAA 07D793DE 4E6D5E5C 94EEE802 89070FB0 5D38FF58 321F2E80 0536D538 CCDAA3D9",
        (EC(B163(T)(1), B163(T)(1)))
    ),
    parse(BigInt, "04 00000000 00000000 00020108 A2E0CC0D 99F8A5EF", base=16),
    2
)

SECT163R1(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "040369 979697AB 43897789 56678956 7F787A78 76A65400 435EDB42 EFAFB298 9D51FEFC E3C80988 F41FF883",
        (EC(
            B163(T)("07 B6882CAA EFA84F95 54FF8428 BD88E246 D2782AE2"),
            B163(T)("07 13612DCD DCB40AAB 946BDA29 CA91F73A F958AFD9")
        ))
    ),
    parse(BigInt, "03 FFFFFFFF FFFFFFFF FFFF48AA B689C29C A710279B", base=16),
    2
)

SECT163R2(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "0403F0 EBA16286 A2D57EA0 991168D4 994637E8 343E3600 D51FBC6C 71A0094F A2CDD545 B11C5C0C 797324F1",
        (EC(
            B163(T)(1),
            B163(T)("02 0A601907 B8C953CA 1481EB10 512F7874 4A3205FD")
        ))
    ),
    parse(BigInt, "04 00000000 00000000 000292FE 77E70C12 A4234C33", base=16),
    2
)

SECT233K1(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "04 017232BA 853A7E73 1AF129F2 2FF41495 63A419C2 6BF50A4C 9D6EEFAD 612601DB 537DECE8 19B7F70F 555A67C4 27A8CD9B F18AEB9B 56E0C11056FAE6A3",
        (EC(
            B233(T)(0),
            B233(T)(1)
        ))
    ),
    parse(BigInt, "80 00000000 00000000 00000000 00069D5B B915BCD4 6EFB1AD5 F173ABDF", base=16),
    4
)

SECT233R1(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "04 00FAC9DF CBAC8313 BB2139F1 BB755FEF 65BC391F 8B36F8F8 EB7371FD 558B0100 6A08A419 03350678 E58528BE BF8A0BEF F867A7CA 36716F7E 01F81052",
        (EC(
            B233(T)(1),
            B233(T)("0066 647EDE6C 332C7F8C 0923BB58 213B333B 20E9CE42 81FE115F 7D8F90AD")
        ))
    ),
    parse(BigInt, "0100 00000000 00000000 00000000 0013E974 E72F8A69 22031D26 03CFE0D7", base=16),
    2
)

SECT283K1(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "04 0503213F 78CA4488 3F1A3B81 62F188E5 53CD265F 23C1567A 16876913 B0C2AC24 58492836 01CCDA38 0F1C9E31 8D90F95D 07E5426F E87E45C0 E8184698 E4596236 4E341161 77DD2259",
        (EC(
            B283(T)(0),
            B283(T)(1)
        ))
    ),
    parse(BigInt, "01FFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFE9AE 2ED07577 265DFF7F 94451E06 1E163C61", base=16),
    4
)

SECT283R1(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "04 05F93925 8DB7DD90 E1934F8C 70B0DFEC 2EED25B8 557EAC9C 80E2E198 F8CDBECD 86B12053 03676854 FE24141C B98FE6D4 B20D02B4 516FF702 350EDDB0 826779C8 13F0DF45 BE8112F4",
        (EC(
            B283(T)(1),
            B283(T)("027B680A C8B8596D A5A4AF8A 19A0303F CA97FD76 45309FA2 A581485A F6263E31 3B79A2F5")
        ))
    ),
    parse(BigInt, "03FFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFEF90 399660FC 938A9016 5B042A7C EFADB307", base=16),
    2
)

SECT409K1(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "04 0060F05F 658F49C1 AD3AB189 0F718421 0EFD0987 E307C84C 27ACCFB8 F9F67CC2 C460189E B5AAAA62 EE222EB1 B35540CF E9023746 01E36905 0B7C4E42 ACBA1DAC BF04299C 3460782F 918EA427 E6325165 E9EA10E3 DA5F6C42 E9C55215 AA9CA27A 5863EC48 D8E0286B",
        (EC(
            B409(T)(0),
            B409(T)(1)
        ))
    ),
    parse(BigInt, "7FFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFE5F 83B2D4EA 20400EC4 557D5ED3 E3E7CA5B 4B5C83B8 E01E5FCF", base=16),
    4
)

SECT409R1(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "04 015D4860 D088DDB3 496B0C60 64756260 441CDE4A F1771D4D B01FFE5B 34E59703 DC255A86 8A118051 5603AEAB 60794E54 BB7996A7 0061B1CF AB6BE5F3 2BBFA783 24ED106A 7636B9C5 A7BD198D 0158AA4F 5488D08F 38514F1F DF4B4F40 D2181B36 81C364BA 0273C706",
        (EC(
            B409(T)(1),
            B409(T)("0021A5C2 C8EE9FEB 5C4B9A75 3B7B476B 7FD6422E F1F3DD67 4761FA99 D6AC27C8 A9A197B2 72822F6C D57A55AA 4F50AE31 7B13545F")
        ))
    ),
    parse(BigInt, "01000000 00000000 00000000 00000000 00000000 00000000 000001E2 AAD6A612 F33307BE 5FA47C3C 9E052F83 8164CD37 D9A21173", base=16),
    2
)

SECT571K1(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "04 026EB7A8 59923FBC 82189631 F8103FE4 AC9CA297 0012D5D4 60248048 01841CA4 43709584 93B205E6 47DA304D B4CEB08C BBD1BA39 494776FB 988B4717 4DCA88C7 E2945283 A01C8972 0349DC80 7F4FBF37 4F4AEADE 3BCA9531 4DD58CEC 9F307A54 FFC61EFC 006D8A2C 9D4979C0 AC44AEA7 4FBEBBB9 F772AEDC B620B01A 7BA7AF1B 320430C8 591984F6 01CD4C14 3EF1C7A3",
        (EC(
            B571(T)(0),
            B571(T)(1)
        ))
    ),
    parse(BigInt, "02000000 00000000 00000000 00000000 00000000 00000000 00000000 00000000 00000000 131850E1 F19A63E4 B391A8DB 917F4138 B630D84B E5D63938 1E91DEB4 5CFE778F 637C1001", base=16),
    4
)

SECT571R1(T::Type{U}) where U<:Unsigned = CurveDomainParams(
    ECPointAffine(
        "04 0303001D 34B85629 6C16C0D4 0D3CD775 0A93D1D2 955FA80A A5F40FC8 DB7B2ABD BDE53950 F4C0D293 CDD711A3 5B67FB14 99AE6003 8614F139 4ABFA3B4 C850D927 E1E7769C 8EEC2D19 037BF273 42DA639B 6DCCFFFE B73D69D7 8C6C27A6 009CBBCA 1980F853 3921E8A6 84423E43 BAB08A57 6291AF8F 461BB2A8 B3531D2F 0485C19B 16E2F151 6E23DD3C 1A4827AF 1B8AC15B",
        (EC(
            B571(T)(1),
            B571(T)("02F40E7E 2221F295 DE297117 B7F3D62F 5C6A97FF CB8CEFF1 CD6BA8CE 4A9A18AD 84FFABBD 8EFA5933 2BE7AD67 56A66E29 4AFD185A 78FF12AA 520E4DE7 39BACA0C 7FFEFF7F 2955727A")
        ))
    ),
    parse(BigInt, "03FFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF E661CE18 FF559873 08059B18 6823851E C7DD9CA1 161DE93D 5174D66E 8382E9BB 2FE84E47", base=16),
    2
)
