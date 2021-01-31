#Guide to ECC, Algorithm 2.41
function reduce(a::BFieldPoint163)::BFieldPoint163
    b::StaticUInt{length(a.value),UInt32} = copy(a.value)
    for i in length(b):-1:7
        t::UInt32 = b.value[i]
        b.value[i-6] ⊻= t<<29
        b.value[i-5] ⊻= (t<<4) ⊻ (t<<3) ⊻ t ⊻ (t>>>3)
        b.value[i-4] ⊻= (t>>28) ⊻ (t>>29)
    end

    t = b.value[6]>>>3
    b.value[1] ⊻= (t<<7) ⊻ (t<<6) ⊻ (t<<3) ⊻ t
    b.value[2] ⊻= (t>>>25) ⊻ (t>>>26)
    b.value[6] &= (UInt32(1)<<3)-1
    return BFieldPoint163(changelength(b, 6))
end

function reduce(a::BFieldPoint233)::BFieldPoint233
    b::StaticUInt{length(a.value),UInt32} = copy(a.value)
    for i in length(b):-1:9
        t::UInt32 = b.value[i]
        b.value[i-8] ⊻= t<<23
        b.value[i-7] ⊻= t>>>9
        b.value[i-5] ⊻= t<<1
        b.value[i-4] ⊻= t>>>31
    end

    t = b.value[8]>>>9
    b.value[1] ⊻= t
    b.value[3] ⊻= t<<10
    b.value[4] ⊻= t>>>22
    b.value[8] &= (UInt32(1)<<9)-1
    return BFieldPoint233(changelength(b, 8))
end

function reduce(a::BFieldPoint283)::BFieldPoint283
    b::StaticUInt{length(a.value),UInt32} = copy(a.value)
    for i in length(b):-1:10
        t::UInt32 = b.value[i]
        b.value[i-9] ⊻= (t<<5) ⊻ (t<<10) ⊻ (t<<12) ⊻ (t<<17)
        b.value[i-8] ⊻= (t>>>27) ⊻ (t>>>22) ⊻ (t>>>20) ⊻ (t>>>15)
    end

    t = b.value[9]>>>27
    b.value[1] ⊻= t ⊻ (t<<5) ⊻ (t<<7) ⊻ (t<<12)
    b.value[9] &= (UInt32(1)<<27)-1
    return BFieldPoint283(changelength(b, 9))
end

function reduce(a::BFieldPoint409)::BFieldPoint409
    b::StaticUInt{length(a.value),UInt32} = copy(a.value)
    for i in length(b):-1:14
        t::UInt32 = b.value[i]
        b.value[i-13] ⊻= t<<7
        b.value[i-12] ⊻= t>>>25
        b.value[i-11] ⊻= t<<30
        b.value[i-10] ⊻= t>>>2
    end

    t = b.value[13]>>>25
    b.value[1] ⊻= t
    b.value[3] ⊻= t<<23
    b.value[13] &= (UInt32(1)<<25)-1
    return BFieldPoint409(changelength(b, 13))
end

function reduce(a::BFieldPoint571)::BFieldPoint571
    b::StaticUInt{length(a.value),UInt32} = copy(a.value)
    for i in length(b):-1:19
        t::UInt32 = b.value[i]
        b.value[i-18] ⊻= (t<<5) ⊻ (t<<7) ⊻ (t<<10) ⊻ (t<<15)
        b.value[i-17] ⊻= (t>>>27) ⊻ (t>>>25) ⊻ (t>>>22) ⊻ (t>>>17)
    end

    t = b.value[18]>>>27
    b.value[1] ⊻= t ⊻ (t<<2) ⊻ (t<<5) ⊻ (t<<10)
    b.value[18] &= (UInt32(1)<<27)-1
    return BFieldPoint571(changelength(b, 18))
end
