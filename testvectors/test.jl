function testcurve(curve::String)
    println("Testing ", curve, ":")
    tests = readlines(".\\testvectors\\"*curve*".txt")
    G = eval(Symbol(curve)).G

    for i in 1:4:length(tests)
        k = parse(BigInt, tests[i][5:length(tests[i])], base=10)
        expectedx = typeof(G.ec.a)(tests[i+1][5:length(tests[i+1])])
        expectedy = typeof(G.ec.a)(tests[i+2][5:length(tests[i+2])])

        print("G*", k)
        Gk = G*k
        if Gk.x==expectedx && Gk.y==expectedy
            println(" is correct.")
        else
            println(" is incorrect.")
            return false
        end
    end
    return true
end
