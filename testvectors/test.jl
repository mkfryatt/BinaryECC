function testcurve(curve::String, verbose=true)
    if verbose println("Testing ", curve, ":") end
    tests = readlines("./testvectors/$curve.txt")
    G = eval(Symbol(curve))(UInt).G

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
