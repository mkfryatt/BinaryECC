#Guide to ECC, Algorithm 2.41
#altered for a 64bit word size and several different fields

function reduce(a::FieldPoint113)
    b = copy(a.value)
    for i in length(b):-1:3
        t = b.value[i]
        b.value[i-1] ⊻= (t>>>40) ⊻ (t>>>49)
        b.value[i-2] ⊻= (t<<24) ⊻ (t<<15)
    end
    t = b.value[2]>>>49
    b.value[1] ⊻= (t<<9) ⊻ t
    b.value[2] &= (UInt64(1)<<49)-1
    return FieldPoint113(changelength(b, 2))
end

function reduce(a::FieldPoint131)
    b = copy(a.value)
    for i in length(b):-1:4
        t = b.value[i]
        b.value[i-1] ⊻= t>>>59
        b.value[i-2] ⊻= (t<<5) ⊻ t ⊻ (t>>>1) ⊻ (t>>>3)
        b.value[i-3] ⊻= (t<<63) ⊻ (t<<61)
    end
    t = b.value[3]>>>3
    b.value[2] ⊻= t>>>56
    b.value[1] ⊻= (t<<8) ⊻ (t<<3) ⊻ (t<<2) ⊻ t
    b.value[3] &= (UInt64(1)<<3)-1
    return FieldPoint131(changelength(b, 3))
end

function reduce(a::FieldPoint163)
    b = copy(a.value)
    for i in length(b):-1:4 #reduce b[i]z^{64(i-1)} mod reduction polynomial
        t = b.value[i]
        b.value[i-2] ⊻= (t>>>28) ⊻ (t>>>29) ⊻ (t>>>32) ⊻ (t>>>35)
        b.value[i-3] ⊻= (t<<36) ⊻ (t<<35) ⊻ (t<<32) ⊻ (t<<29)
    end
    #reduce the extra bits of b (from word 3)
    t = b.value[3]>>>35 #extract bits 35-63 from b[3], i.e. bits 163-191 of b
    b.value[1] ⊻= (t<<7) ⊻ (t<<6) ⊻ (t<<3) ⊻ t
    b.value[3] &= (UInt64(1)<<35)-1 #clear the reduced bits of b[3]
    return FieldPoint163(changelength(b, 3))
end

function reduce(a::FieldPoint193)
    b = copy(a.value)
    for i in length(b):-1:5
        t = b.value[i]
        b.value[i-2] ⊻= t>>>50
        b.value[i-3] ⊻= (t<<14) ⊻ (t>>>1)
        b.value[i-4] ⊻= t<<63
    end
    t = b.value[4]>>>1
    b.value[1] ⊻= t ⊻ (t<<15)
    b.value[2] ⊻= t>>>49
    b.value[4] &= 1
    return FieldPoint193(changelength(b, 4))
end

function reduce(a::FieldPoint233)
    b = copy(a.value)
    for i in length(b):-1:5
        t = b.value[i]
        b.value[i-2] ⊻= t>>>31
        b.value[i-3] ⊻= (t<<33) ⊻ (t>>>41)
        b.value[i-4] ⊻= t<<23
    end
    t = b.value[4]>>>41
    b.value[1] ⊻= t
    b.value[2] ⊻= t<<10
    b.value[4] &= (UInt64(1)<<41)-1
    return FieldPoint233(changelength(b, 4))
end

function reduce(a::FieldPoint239)
    b = copy(a.value)
    for i in length(b):-1:5
        t = b.value[i]
        b.value[i-3] ⊻= (t>>>11) ⊻ (t>>>47)
        b.value[i-4] ⊻= (t<<53) ⊻ (t<<17)
    end
    t = b.value[4]>>>47
    b.value[1] ⊻= t ⊻ (t<<36)
    b.value[4] &= (UInt64(1)<<47)-1
    return FieldPoint239(changelength(b, 4))
end

function reduce(a::FieldPoint283)
    b = copy(a.value)
    for i in length(b):-1:6
        t = b.value[i]
        b.value[i-4] ⊻= (t>>>15) ⊻ (t>>>20) ⊻ (t>>>22) ⊻ (t>>>27)
        b.value[i-5] ⊻= (t<<49) ⊻ (t<<44) ⊻ (t<<42) ⊻ (t<<37)
    end
    t = b.value[5]>>>27
    b.value[1] ⊻= t ⊻ (t<<5) ⊻ (t<<7) ⊻ (t<<12)
    b.value[5] &= (UInt64(1)<<27)-1
    return FieldPoint283(changelength(b, 5))
end

function reduce(a::FieldPoint409)
    b = copy(a.value)
    for i in length(b):-1:8
        t = b.value[i]
        b.value[i-5] ⊻= t>>>2
        b.value[i-6] ⊻= (t<<62) ⊻ (t>>>25)
        b.value[i-7] ⊻= t<<39
    end
    t = b.value[7]>>>25
    b.value[1] ⊻= t
    b.value[2] ⊻= t<<23
    b.value[7] &= (UInt64(1)<<25)-1
    return FieldPoint409(changelength(b, 7))
end

function reduce(a::FieldPoint571)
    b = copy(a.value)
    for i in length(b):-1:10
        t = b.value[i]
        b.value[i-8] ⊻= (t>>>49) ⊻ (t>>>54) ⊻ (t>>>57) ⊻ (t>>>59)
        b.value[i-9] ⊻= (t<<15) ⊻ (t<<10) ⊻ (t<<7) ⊻ (t<<5)
    end
    t = b.value[9]>>>59
    b.value[1] ⊻= t ⊻ (t<<2) ⊻ (t<<5) ⊻ (t<<10)
    b.value[9] &= (UInt64(1)<<59)-1
    return FieldPoint571(changelength(b, 9))
end
