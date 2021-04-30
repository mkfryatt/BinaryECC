using Test

#155 tests

for size in [UInt8, UInt16, UInt32, UInt64, UInt128]
    @testset "Wordsize $size" begin
        f = B163(size)
        x = random(f)
        y = random(f)
        @test x*y == y*x
        @test x+y == y+x
        @test iszero(x+x)
        @test iszero(y-y)
        @test x^6 == x*x*x*x*x*x
        @test x^5 == x*x*x*x*x
        @test square(x) == x^2 == x*x
        @test x*inv(y) == x/y
        @test isone(x/x)
        @test square(sqrt(x)) == x
    end
end

for field in [B113, B131, B163, B193, B233, B239, B283, B409, B571]
    @testset "Field $field" begin
        f = field(UInt)
        x = random(f)
        y = random(f)
        @test x*y == y*x
        @test x+y == y+x
        @test iszero(x+x)
        @test iszero(y-y)
        @test x^6 == x*x*x*x*x*x
        @test x^5 == x*x*x*x*x
        @test square(x) == x^2 == x*x
        @test x*inv(y) == x/y
        @test isone(x/x)
        @test square(sqrt(x)) == x
    end
end

@testset "Alternative multiplication" begin
    field = B163(UInt)
    x = random(field)
    y = random(field)
    z = x*y
    @test z == mult_shiftandadd(x, y)
    @test z == mult_threaded(x, y)
    @test z == mult_ownreduce(x, y)
    @test z == mult_comb_rtl(x, y)
    @test z == mult_comb_ltr(x, y)
    for w in [2,4,8]
        @test z == mult_shiftandadd_window(x, y, w)
        @test z == mult_comb_window(x, y, w)
    end
end

@testset "Alternative squaring" begin
    field = B163(UInt)
    x = random(field)
    x2 = square(x)
    @test x2 == square_standard(x)
    for w in [2,4,8]
        @test x2 == square_window(x, w)
    end
end
