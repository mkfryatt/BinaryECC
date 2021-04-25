using Test

function testcurve(curve, verbose=true)
    if verbose println("Testing $curve:") end
    tests = readlines("./test/$curve.txt")
    G = curve(UInt).G

    for i in 1:4:length(tests)
        k = parse(BigInt, tests[i][5:length(tests[i])], base=10)
        expectedx = typeof(G.ec.a)(tests[i+1][5:length(tests[i+1])])
        expectedy = typeof(G.ec.a)(tests[i+2][5:length(tests[i+2])])

        if verbose print("G*$k") end
        Gk = G*k
        if Gk.x==expectedx && Gk.y==expectedy
            if verbose println(" is correct.") end
        else
            if verbose println(" is incorrect.") end
            return false
        end
    end
    return true
end

function testecdsa(verbose=true)
    #https://www.ietf.org/rfc/rfc6979.txt
    if verbose println("Testing ECDSA") end
    tests = readlines("./test/ecdsa.txt")
    tests = [split(t, " = ")[2] for t in tests]

    i = 1
    while i<length(tests)
        if verbose println("Curve: $(tests[i])") end

        T = eval(Meta.parse(tests[i]))(UInt)
        privkey = parse(BigInt, tests[i+1], base=16)
        U = ECKeyPair(PFieldElt(privkey, T.n), T.G*privkey)

        Ux = convert(BigInt, U.Q.x)
        Uy = convert(BigInt, U.Q.y)
        if parse(BigInt, tests[i+2], base=16) != Ux ||
            parse(BigInt, tests[i+3], base=16) != Uy
            if verbose println("Failed to compute public key") end
            return false
        end

        i += 4
        for j in 1:2
            msg = string(tests[i])
            k = string(tests[i+1])
            sig = ecdsa_sign_deterministic(T, U, msg, k)
            r = parse(BigInt, tests[i+2], base=16)
            s = parse(BigInt, tests[i+3], base=16)
            if r!=sig.r.value || s!=sig.s.value
                if verbose println("Failed to compute signature") end
                return false
            end
            i += 4
        end
    end
    return true
end

@testset "ECDSA testvectors" begin
    @test testecdsa(false)
    msg1 = "message 1"
    msg2 = "message 2"
    sect163k1 = SECT163K1(UInt)
    ukey = generate_keypair(sect163k1)
    vkey = generate_keypair(sect163k1)
    @test ukey != vkey
    sig1_real = ecdsa_sign(sect163k1, ukey, msg1)
    sig1_fake = ecdsa_sign(sect163k1, vkey, msg1)
    @test sig1_real!=sig1_fake
    @test ecdsa_verify(sect163k1, ukey.Q, sig1_real, msg1)
    @test !ecdsa_verify(sect163k1, ukey.Q, sig1_fake, msg1)
    @test !ecdsa_verify(sect163k1, ukey.Q, sig1_fake, msg2)
end

@testset "Binary Field" begin
    f = B163(UInt)
    x = random(f)
    y = random(f)
    @test x*y == y*x #test commutativity for mult
    @test x+y == y+x #test commutativity for add
    @test iszero(x+x) #test doubling
    @test x^5 == x*x*x*x*x #test exponentiation for an odd power
    @test x^6 == x*x*x*x*x*x #test exponentiation for an even power
    @test x*inv(y) == x/y #test division is inv and mult
    @test isone(x/x) #test division
end

@testset "Prime Field" begin
    sect163k1 = SECT163K1(UInt)
    x = random(PFieldElt, sect163k1.n)
    y = random(PFieldElt, sect163k1.n)
    @test x*y == y*x #test commutativity for mult
    @test x+y == y+x #test commutativity for add
    @test x^5 == x*x*x*x*x #test exponentiation for an odd power
    @test x^6 == x*x*x*x*x*x #test exponentiation for an even power
    @test x*inv(y) == x/y #test division is inv and mult
    @test isone(x/x) #test division
end

@testset "Elliptic Curve" begin
    sect163k1 = SECT163K1(UInt)
    G = sect163k1.G
    sect163r2 = SECT163R2(UInt)
    G_other = sect163r2.G
    O = zero(typeof(G), G.ec)
    @test iszero(O) #test additive identity is correct
    @test G-G == O #test subtraction
    @test G*5 == G+G+G+G+G #test scalar mult
    G2 = G+G
    @test G2+G == G+G2 #test addition of different points
    @test_throws ECMismatchException G+G_other #test different curves
end

@testset "Elliptic Curve Testvectors" begin
    @test testcurve(SECT163K1, false) #test the smallest curve
    @test testcurve(SECT571R1, false) #test the largest curve
end

@testset "Jacobian Coordinates" begin
    sect163k1 = SECT163K1(UInt)
    G = sect163k1.G
    G_J = convert(ECPointJacobian, G)
    @test G == G_J #test the conversion worked
    @test G_J*2 == G*2 #test doubling
    @test G_J*3 == G*3 #test addition
end

@testset "Lopez-Dahab Coordinates" begin
    sect163k1 = SECT163K1(UInt)
    G = sect163k1.G
    G_LD = convert(ECPointLD, G)
    @test G == G_LD #test the conversion worked
    @test G_LD*2 == G*2 #test doubling
    @test G_LD*3 == G*3 #test addition
end

@testset "ECDH" begin
    sect163k1 = SECT163K1(UInt)
    u1 = ecdh_deployment1(sect163k1)
    v1 = ecdh_deployment1(sect163k1)
    @test u1!=v1
    @test ecdh_deployment2(sect163k1, v1.Q)
    @test ecdh_deployment2(sect163k1, u1.Q)
    @test ecdh_agreement(sect163k1, u1, v1.Q)==ecdh_agreement(sect163k1, v1, u1.Q)
end
