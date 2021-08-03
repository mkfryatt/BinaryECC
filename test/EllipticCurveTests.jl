using Test

#89 tests

function testcurve(curve, verbose=true)
    if verbose println("Testing $curve:") end
    tests = readlines("$curve.txt")
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

for curve in [SECT163K1, SECT163R2, SECT233K1, SECT233R1, SECT283K1,
    SECT283R1, SECT409K1, SECT409R1, SECT571K1, SECT571R1]
    @testset "Curve $curve" begin
        T = curve(UInt)
        G = T.G
        n = T.n
        @test iszero(G*n)
        @test iszero(G-G)
        @test G+G == double(G)
        testcurve(curve, false)
    end
end

@testset "Jacobian Coordinates" begin
    T = SECT163K1(UInt)
    G = T.G
    GJ = convert(ECPointJacobian, G)
    n = T.n
    @test iszero(GJ*n)
    @test iszero(GJ-GJ)
    @test GJ+GJ == double(GJ)
    x = rand(1:n)
    @test G*n == GJ*n
end

@testset "Lopez-Dahab Coordinates" begin
    T = SECT163K1(UInt)
    G = T.G
    GLD = convert(ECPointLD, G)
    n = T.n
    @test iszero(GLD*n)
    @test iszero(GLD-GLD)
    @test GLD+GLD == double(GLD)
    x = rand(1:n)
    @test G*n == GLD*n
end

@testset "Alternative doubling" begin
    T = SECT163K1(UInt)
    G = T.G
    G2 = double(G)
    @test G2 == double_threaded(G)
    @test G2 == double_memo(G)
    @test G2 ==  double_standard(G)
end

@testset "Alternative scalar multiplication" begin
    T = SECT163K1(UInt)
    G = T.G
    n = T.n
    for x in [1, rand(1:n), n]
        Q = G*x
        @test Q == mult_standard_ltr(G, x)
        @test Q == mult_standard_rtl(G, x)
        @test Q == mult_bnaf_threaded(G, x)
        @test Q == mult_bnaf(G, x)
        @test Q == mult_memo(G, x)
        for w in [2,4,8]
            @test Q == mult_window(G, x, w)
            @test Q == mult_bnaf_window(G, x, w)
            @test Q == mult_wnaf(G, x, w)
        end
    end
end

@testset "Montgomery's powering ladder" begin
    T = SECT163K1(UInt)
    G = T.G
    n = T.n
    for x in [1, rand(1:n), n]
        Q = G*x
        @test Q == mult_mont_general(G, x)
        @test Q == mult_mont_affine(G, x)
    end
end
