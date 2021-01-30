include("..\\testvectors\\test.jl")

using Test

@testset "Field" begin
    x = random(FieldPoint163)
    y = random(FieldPoint163)
    @test x*y == y*x
    @test x+y == y+x
    @test iszero(x+x)
    @test x^5 == x*x*x*x*x
    @test x^6 == x*x*x*x*x*x
    @test x*inv(y) == x/y
    @test isone(x/x)
end

@testset "Elliptic Curve" begin
    G = SECT163K1.G
    G_other = SECT163R2.G
    O = zero(typeof(G), G.ec)
    @test iszero(O)
    @test G-G == O
    @test G*5 == G+G+G+G+G
    G2 = G+G
    @test G2+G == G+G2
    @test_throws ECMismatchException G+G_other
end

@testset "Elliptic Curve Testvectors" begin
    @test testcurve("SECT163K1")
    @test testcurve("SECT571R1")
end
