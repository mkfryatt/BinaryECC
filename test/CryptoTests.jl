using Test

#107 tests

function testecdsa(verbose=true)
    #https://www.ietf.org/rfc/rfc6979.txt
    if verbose println("Testing ECDSA") end
    tests = readlines("ecdsa.txt")
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
            sig = Main.BinaryECC.ecdsa_sign_deterministic(T, U, msg, k)
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

@testset "Prime Field" begin
    T = SECT163K1(UInt)
    n = T.n
    x = random(PFieldElt, n)
    y = random(PFieldElt, n)
    @test x*y == y*x
    @test x+y == y+x
    @test x^5 == x*x*x*x*x
    @test x^6 == x*x*x*x*x*x
    @test x*inv(y) == x/y
    @test isone(x/x)
    @test iszero(x-x)
end

@testset "ECDSA Test Vector" begin
    @test testecdsa(false)
end

for curve in [SECT163K1, SECT163R1, SECT163R2, SECT233K1, SECT233R1, SECT283K1,
    SECT283R1, SECT409K1, SECT409R1, SECT571K1, SECT571R1]
    @testset "ECDSA General: curve $curve" begin
        msg1 = "message 1"
        msg2 = "message 2"
        T = curve(UInt)
        ukey = generate_keypair(T)
        vkey = generate_keypair(T)
        @test ukey != vkey
        sig1_real = ecdsa_sign(T, ukey, msg1)
        sig1_fake = ecdsa_sign(T, vkey, msg1)
        @test sig1_real!=sig1_fake
        @test ecdsa_verify(T, ukey.Q, sig1_real, msg1)
        @test !ecdsa_verify(T, ukey.Q, sig1_fake, msg1)
        @test !ecdsa_verify(T, ukey.Q, sig1_fake, msg2)
    end

    @testset "ECDH General: curve $curve" begin
        T = curve(UInt)
        u1 = ecdh_deployment1(T)
        v1 = ecdh_deployment1(T)
        @test u1!=v1
        @test ecdh_deployment2(T, v1.Q)
        @test ecdh_deployment2(T, u1.Q)
        @test ecdh_agreement(T, u1, v1.Q)==ecdh_agreement(T, v1, u1.Q)
    end
end
