function testcurve(curve::String)
    println("Testing ", curve, ":")
    tests = readlines("D:\\OneDrive - University Of Cambridge\\diss_code\\BinaryECC\\testvectors\\"*curve*".txt")
    G = eval(Symbol(curve)).G

    for i in 1:4:length(tests)
        k = parse(BigInt, tests[i][5:length(tests[i])], base=10)

        expectedx, expectedy = 0,0

        try
            expectedx = typeof(G.ec.a)(tests[i+1][5:length(tests[i+1])])
            expectedy = typeof(G.ec.a)(tests[i+2][5:length(tests[i+2])])
        catch e
            println("k: ", k)
            println("ERROR: ", e)
            if e isa ArgumentError
                continue
            else
                break
            end
        end

        print("G*", k)
        Gk = G*k
        if Gk.x==expectedx && Gk.y==expectedy
            println(" is correct.")
        else
            println(" is incorrect.")
            break
        end
    end
end
