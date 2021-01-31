include("..\\testvectors\\test.jl")

using Test

@testset "Field" begin
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
    @test testcurve("SECT163K1") #test the smallest curve
    @test testcurve("SECT571R1") #test the largest curve
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
