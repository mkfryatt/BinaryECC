import Base: +, -, *, /, ^, ==, repr, inv, sqrt, iszero, isone, convert
using StaticArrays

struct FieldMismatchException <: Exception end

struct Field
    degree::UInt16 #the degree of the reduction polynomial
    reduction::UInt128 #the reduction polynomial
    Field(d::Integer, r::Integer) = new(convert(UInt16, d), convert(UInt128, r))
end

function ==(a::Field, b::Field)
    return a.degree==b.degree && a.reduction==b.reduction
end

struct FieldPoint
    value::MVector
    field::Field
    FieldPoint(x::MVector, field::Field) = new(x, field)
    FieldPoint(x::Integer, field::Field) = new(tovector(x, field.degree), field)
end

#sec1v2 2.3.6
#convert a hex string to a field element
function FieldPoint(s::String, f::Field)
    s = replace(s, " " => "")
    if length(s)!=ceil(f.degree / 8)*2 throw(ArgumentError("Octet string is of the incorrect length for this field.")) end
    value = parse(BigInt, s, base=16)
    return FieldPoint(value, f)
end

function repr(a::FieldPoint)
    return repr(a)
end

function ==(a::FieldPoint, b::FieldPoint)
    return a.value==b.value && a.field==b.field
end

function +(a::FieldPoint, b::FieldPoint)
    if a.field!=b.field throw(FieldMismatchException()) end
    c = FieldPoint(copy(a.value), a.field)
    for i in 1:length(c.value)
        c.value[i] ⊻= b.value[i]
    end
    return c
end

function -(a::FieldPoint, b::FieldPoint)
    return a+b
end

function -(a::FieldPoint)
    return a
end

function reduce(a::FieldPoint)
    #b will should always be such that a ≡ b (mod R)
    #the loop will modify it until it reaches the smallest value that makes that true
    bvec = copy(a.value)

    #iterate over the excess bits of a, left to right
    for i in (64*length(a.value)-1):-1:a.field.degree
        if getbit(bvec, i)==1
            flipbit!(bvec, i)
            #b ⊻= R<<(i-D)
            startblock = 1 + ((i-a.field.degree)÷64)
            lowerbits = 64 - ((i-a.field.degree) % 64)
            lowermask = (1<<lowerbits)-1
            middlemask = typemax(UInt64)
            uppermask = (1<<(64-lowerbits))-1

            bvec[startblock] ⊻= (a.field.reduction & lowermask)<<(64-lowerbits)
            bvec[startblock+1] ⊻= (a.field.reduction>>lowerbits) & middlemask
            bvec[startblock+2] ⊻= (a.field.reduction>>(64+lowerbits)) & uppermask
        end
    end

    #make a new vector without the spare blocks in it
    shortenedbvec = zeros(MVector{ceil(Int, a.field.degree/64), UInt64})
    for block in 1:length(shortenedbvec)
        shortenedbvec[block] = bvec[block]
    end

    return FieldPoint(shortenedbvec, a.field)
end

#right to left, shift and add
function *(a::FieldPoint, b::FieldPoint)
    if a.field!=b.field throw(FieldMismatchException()) end
    if a.value==b.value return square(a) end

    c = zeros(MVector{2*length(a.value), UInt64})

    for i in 0:(a.field.degree-1)
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

    return reduce(FieldPoint(c, a.field))
end

#add a zero between every digit of the original
function square(a::FieldPoint)
    b = zeros(MVector{2*length(a.value), UInt64})
    for i in 0:(a.field.degree-1)
        if getbit(a.value, i)==1
            flipbit!(b, 2*i)
        end
    end
    return reduce(FieldPoint(b, a.field))
end

#uses a version of egcd to invert a
#Algorithm 2.48, Guide to Elliptic Curve Cryptography
function inv(a::FieldPoint)
    if a.value==0 throw(DivideError()) end

    u = copy(a.value)
    v = tovector(a.field.reduction, a.field.degree)
    flipbit!(v, a.field.degree)
    g1 = zeros(typeof(a.value))
    flipbit!(g1, 0)
    g2 = zeros(typeof(a.value))

    while !isone(u)
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
        uppermask = typemax(UInt64) ⊻ lowermask
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
    return FieldPoint(g1, a.field)
end

function /(a::FieldPoint, b::FieldPoint)
    if a.field!=b.field throw(FieldMismatchException()) end
    return a * inv(b)
end

#right to left, square and multiply method
function ^(a::FieldPoint, b::Integer)
    c = FieldPoint(1, a.field)
    squaring = a

    while b>0
        if b&1 == 1
            c *= squaring
        end
        squaring = square(squaring)
        b >>>= 1
    end

    return c
end

#return a random element of the specified field
function random(f::Field)
    range = BigInt(0):((BigInt(1)<<f.degree)-BigInt(1))
    return FieldPoint(rand(range), f)
end

function iszero(a::FieldPoint)
    return iszero(a.value)
end

#sec1 v2, 2.3.9
function convert(::Type{BigInt}, a::FieldPoint)
    b = BigInt(0)
    for i in 1:length(a.value)
        b += BigInt(a.value[i])<<(64*(i-1))
    end
    return b
end

#stores a number ..., b191, ..., b1, b0 as a vector of:
#b63, b62, ..., b2, b1, b0
#b127, b126, ..., b66, b65, b64
#b191, b190, ..., b130, b129, b128
#...
#i.e. little endian on 64bit block basis (rather than byte by byte)
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

function isone(vec::MVector)
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

#number of bits in the binary representation of this number
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

#sec2 v2 (and v1), table 3:
const FIELD113 = Field(113, UInt128(512+1)) #v1 only
const FIELD131 = Field(131, UInt128(256+8+4+1)) #v1 only
const FIELD163 = Field(163, UInt128(128+64+8+1))
const FIELD193 = Field(193, (UInt128(1)<<15) + UInt128(1)) #v1 only
const FIELD233 = Field(233, (UInt128(1)<<74) + UInt128(1))
const FIELD239 = Field(239, (UInt128(1)<<36) + UInt128(1))
const FIELD283 = Field(283, (UInt128(1)<<12) + UInt128(128+32+1))
const FIELD409 = Field(409, (UInt128(1)<<87) + UInt128(1))
const FIELD571 = Field(571, (UInt128(1)<<10) + UInt128(32+4+1))
