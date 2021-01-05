import Base: +, -, *, /, ^, ==, repr, inv, sqrt, iszero, convert
using StaticArrays

#D is the degree of the reduction polynomial
#R is the reduction polynomial without the x^D term
struct FieldPoint{D,R}
    value::MVector
    FieldPoint{D,R}(value::MVector) where {D,R} = new(value)
    FieldPoint{D,R}(value::Integer) where {D,R} = new(tovector(value, D))
end

#sec1v2 2.3.6
#convert a hex string to a field element
function FieldPoint{D,R}(s::String) where {D,R}
    s = replace(s, " " => "")
    if length(s)!=ceil(D / 8)*2 throw(ArgumentError("Octet string is of the incorrect length for this field.")) end
    value = parse(BigInt, s, base=16)
    return FieldPoint{D,R}(value)
end

#stores a number ..., b191, ..., b1, b0 as a vector of:
#b63, b62, ..., b2, b1, b0
#b127, b126, ..., b66, b65, b64
#b191, b190, ..., b130, b129, b128
#...
function tovector(value::Integer, D::Integer)
    numberofblocks = ceil(Int, D/64)
    valuevector = zeros(MVector{numberofblocks, UInt64})
    bitmask = UInt64(0) -1
    for i in 1:numberofblocks
        valuevector[i] = UInt64(value & bitmask)
        value >>= 64
    end
    return valuevector
end

function test(vec::MVector)
    if vec[1]!=1
        return false
    end
    for i in 2:length(vec)
        if vec[i]!=0
            return false
        end
    end
    return true
end

function getbit(vec::MVector, i::Integer)
    bit = i%64
    block = (i÷64) +1
    return vec[block]>>bit & 1
end

function flipbit!(vec::MVector, i::Integer)
    bit = i%64
    block = (i÷64) +1
    vec[block] ⊻= 1<<bit
end

function repr(a::FieldPoint)
    return repr(a.value)
end

function ==(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    return a.value==b.value
end

#assumes that a and b are both in the given field
function +(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    c = FieldPoint{D,R}(copy(a.value))
    for i in 1:length(c.value)
        c.value[i] ⊻= b.value[i]
    end
    return c
end

function -(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    return a+b
end

function -(a::FieldPoint{D,R}) where {D,R}
    return a
end

function reduce(a::FieldPoint{D,R}) where {D,R}
    #b will should always be such that a ≡ b (mod R)
    #the loop will modify it until it reaches the smallest value that makes that true
    bvec = copy(a.value)

    #iterate over the excess bits of a, left to right
    for i in (64*length(a.value)-1):-1:D
        if getbit(bvec, i)==1
            flipbit!(bvec, i)
            #b ⊻= R<<(i-D)
            startblock = 1 + ((i-D)÷64)
            lowerbits = 64 - ((i-D) % 64)
            lowermask = (1<<lowerbits)-1
            middlemask = UInt64(0)-1
            uppermask = (1<<(64-lowerbits))-1

            bvec[startblock] ⊻= (R & lowermask)<<(64-lowerbits)
            bvec[startblock+1] ⊻= (R>>lowerbits) & middlemask
            bvec[startblock+2] ⊻= (R>>(64+lowerbits)) & uppermask
        end
    end

    #make a new vector without the spare blocks in it
    shortenedbvec = zeros(MVector{ceil(Int, D/64), UInt64})
    for block in 1:length(shortenedbvec)
        shortenedbvec[block] = bvec[block]
    end

    return FieldPoint{D,R}(shortenedbvec)
end

#right to left, shift and add
function *(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    if a.value==b.value return square(a) end

    c = zeros(MVector{2*length(a.value), UInt64})

    for i in 0:(D-1)
        if getbit(a.value, i)==1
            #c += b<<i
            startblock = 1+(i÷64)
            upperbits = i%64
            lowerbits = 64-upperbits
            lowermask = (1<<lowerbits)-1
            uppermask = (1<<upperbits)-1

            for block in 1:length(b.value)
                c[startblock+block-1] ⊻= (b.value[block] & lowermask)<<upperbits
                c[startblock+block] ⊻= (b.value[block]>>lowerbits) & uppermask
            end
        end
    end

    return reduce(FieldPoint{D,R}(c))
end

#add a zero between every digit of the original
function square(a::FieldPoint{D,R}) where {D,R}
    b = zeros(MVector{2*length(a.value), UInt64})
    for i in 0:(D-1)
        if getbit(a.value, i)==1
            flipbit!(b, 2*i)
        end
    end
    return reduce(FieldPoint{D,R}(b))
end

#length of a number in bits
function bits(a::MVector)
    block = length(a)
    while a[block]==0
        block -= 1
        if block==0
            return 0
        end
    end

    i = 0
    for i in 0:64
        if a[block] == (UInt64(1)<<i)
            return i+1 + 64*(block-1)
        elseif a[block] < (UInt64(1)<<i)
            return i + 64*(block-1)
        end
    end
    return 64*block
end

#uses a version of egcd to invert a
#Algorithm 2.48, Guide to Elliptic Curve Cryptography
function inv(a::FieldPoint{D,R}) where {D,R}
    if a.value==0 throw(DivideError()) end

    u = copy(a.value)
    v = tovector(R, D)
    flipbit!(v, D)
    g1 = zeros(typeof(a.value))
    flipbit!(g1, 0)
    g2 = zeros(typeof(a.value))

    while !test(u)
        j = bits(u) - bits(v)
        if j<0
            (u, v) = (v, u)
            (g1, g2) = (g2, g1)
            j = -j
        end

        #u ⊻= v << j
        #g1 ⊻= g2 << j
        upperbits = j%64
        lowerbits = 64-upperbits
        lowermask = (1<<lowerbits) -1
        uppermask = (UInt64(0)-1) ⊻ lowermask
        startblock = (j÷64+1)

        u[startblock] ⊻= (v[1]&lowermask)<<upperbits
        g1[startblock] ⊻= (g2[1]&lowermask)<<upperbits

        vblock = 2
        for ublock in (startblock+1):length(u)
            u[ublock] ⊻= (v[vblock-1]&uppermask)>>lowerbits
            g1[ublock] ⊻= (g2[vblock-1]&uppermask)>>lowerbits
            u[ublock] ⊻= (v[vblock]&lowermask)<<upperbits
            g1[ublock] ⊻= (g2[vblock]&lowermask)<<upperbits
            vblock += 1
        end
    end
    return FieldPoint{D,R}(g1)
end

function /(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}
    return a * inv(b)
end

#right to left, square and multiply method
function ^(a::FieldPoint{D,R}, b::Integer) where {D,R}
    c = FieldPoint{D,R}(1)
    squaring = a

    while b>BigInt(0)
        if b & BigInt(1) == BigInt(1)
            c *= squaring
        end
        squaring *= squaring
        b >>>= 1
    end

    return c
end

function random(::Type{FieldPoint{D,R}}) where {D,R}
    range = BigInt(0):((BigInt(1)<<D)-BigInt(1))
    return FieldPoint{D,R}(rand(range))
end

function iszero(a::FieldPoint)
    for block in a.value
        if block!=0
            return false
        end
    end
    return true
end

#sec1 v2, 2.3.9
function convert(::Type{BigInt}, a::FieldPoint)
    b = BigInt(0)
    for i in 1:length(a.value)
        b += BigInt(a.value[i])<<(64*(i-1))
    end
    return b
end

#sec2 v2 (and v1), table 3:
FieldPoint113 = FieldPoint{113, Int128(512+1)} #v1 only
FieldPoint131 = FieldPoint{131, Int128(256+8+4+1)} #v1 only
FieldPoint163 = FieldPoint{163, Int128(128+64+8+1)}
FieldPoint193 = FieldPoint{193, (Int128(1)<<15) + Int128(1)} #v1 only
FieldPoint233 = FieldPoint{233, (Int128(1)<<74) + Int128(1)}
FieldPoint239 = FieldPoint{239, (Int128(1)<<36) + Int128(1)}
FieldPoint283 = FieldPoint{283, (Int128(1)<<12) + Int128(128+32+1)}
FieldPoint409 = FieldPoint{409, (Int128(1)<<87) + Int128(1)}
FieldPoint571 = FieldPoint{571, (Int128(1)<<10) + Int128(32+4+1)}

function testmul(fs)
    for f1 in fs[1:10]
        for f2 in fs[11:20]
            f1*f2
        end
    end
end

function testsq(fs)
    for f in fs
        square(f)
    end
end

function testinv(fs)
    for f in fs
        inv(f)
    end
end
