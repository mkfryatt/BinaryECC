include("..\\testvectors\\test.jl")

using Test

@testset "Binary Field" begin
    x = random(BFieldPoint163)
    y = random(BFieldPoint163)
    @test x*y == y*x #test commutativity for mult
    @test x+y == y+x #test commutativity for add
    @test iszero(x+x) #test doubling
    @test x^5 == x*x*x*x*x #test exponentiation for an odd power
    @test x^6 == x*x*x*x*x*x #test exponentiation for an even power
    @test x*inv(y) == x/y #test division is inv and mult
    @test isone(x/x) #test division
end

@testset "Prime Field" begin
    x = random(PFieldPoint, SECT163K1.n)
    y = random(PFieldPoint, SECT163K1.n)
    @test x*y == y*x #test commutativity for mult
    @test x+y == y+x #test commutativity for add
    @test x^5 == x*x*x*x*x #test exponentiation for an odd power
    @test x^6 == x*x*x*x*x*x #test exponentiation for an even power
    @test x*inv(y) == x/y #test division is inv and mult
    @test isone(x/x) #test division
end

@testset "Elliptic Curve" begin
    G = SECT163K1.G
    G_other = SECT163R2.G
    O = zero(typeof(G), G.ec)
    @test iszero(O) #test additive identity is correct
    @test G-G == O #test subtraction
    @test G*5 == G+G+G+G+G #test scalar mult
    G2 = G+G
    @test G2+G == G+G2 #test addition of different points
    @test_throws ECMismatchException G+G_other #test different curves
end

@testset "Elliptic Curve Testvectors" begin
    @test testcurve("SECT163K1", false) #test the smallest curve
    @test testcurve("SECT571R1", false) #test the largest curve
end

@testset "Jacobian Coordinates" begin
    G = SECT163K1.G
    G_J = convert(ECPointJacobian, G)
    @test G == G_J #test the conversion worked
    @test G_J*2 == G*2 #test doubling
    @test G_J*3 == G*3 #test addition
end

@testset "Lopez-Dahab Coordinates" begin
    G = SECT163K1.G
    G_LD = convert(ECPointLD, G)
    @test G == G_LD #test the conversion worked
    @test G_LD*2 == G*2 #test doubling
    @test G_LD*3 == G*3 #test addition
end

@testset "ECDSA" begin
    msg1 = "message 1"
    msg2 = "message 2"
    ukey = generate_keypair(SECT163K1)
    vkey = generate_keypair(SECT163K1)
    @test ukey != vkey
    sig1_real = ecdsa_sign(SECT163K1, ukey, msg1)
    sig1_fake = ecdsa_sign(SECT163K1, vkey, msg1)
    @test sig1_real!=sig1_fake
    @test ecdsa_verify(SECT163K1, ukey.Q, sig1_real, msg1)
    @test !ecdsa_verify(SECT163K1, ukey.Q, sig1_fake, msg1)
    @test !ecdsa_verify(SECT163K1, ukey.Q, sig1_fake, msg2)
end

@testset "ECDH" begin
    u1 = ecdh_deployment1(SECT163K1)
    v1 = ecdh_deployment1(SECT163K1)
    @test u1!=v1
    @test ecdh_deployment2(SECT163K1, v1.Q)
    @test ecdh_deployment2(SECT163K1, u1.Q)
    @test ecdh_agreement(SECT163K1, u1, v1.Q)==ecdh_agreement(SECT163K1, v1, u1.Q)
end
