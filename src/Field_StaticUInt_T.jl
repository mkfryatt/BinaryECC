#D is the degree of the reduction polynomial
#R is the reduction polynomial without the x^D term
"""
    FieldPoint{D,R,T}
Represents a point in the binary field which has order ``2^D`` and reduction polynomial

``x^D + x^{r_n} + \\cdots + x^{r_0}``

where ``R = r_n r_{n-1}\\ldots r_1 r_0`` in binary.

It stores the polynomial in an array of words of type T.

Types for points in the standard fields (taken from SEC 2, table 3)
 are available:
- FieldPoint113
- FieldPoint131
- FieldPoint163
- FieldPoint193
- FieldPoint233
- FieldPoint239
- FieldPoint283
- FieldPoint409
- FieldPoint571
"""
struct FieldPoint{D,R,T}
    value::StaticUInt
    FieldPoint{D,R,T}(value::Integer) where {D,R,T} = new(StaticUInt{ceil(Int,D/(8*sizeof(T))),T}(value))
    FieldPoint{D,R,T}(value::StaticUInt) where {D,R,T} = new(value)
end

"""
    FieldPoint{D,R,T}(s::String) where {D,R,T}
Using the procedure set out in SEC 1 (version 2) 2.3.6,
this converts a hex string to a field element.
"""
function FieldPoint{D,R,T}(s::String) where {D,R,T}
    s = replace(s, " " => "")
    if length(s)!=ceil(D / 8)*2 throw(ArgumentError("Octet string is of the incorrect length for this field.")) end
    value = StaticUInt{ceil(Int,D/(8*sizeof(T))),T}(s)
    return FieldPoint{D,R,T}(value)
end


"""
    ==(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T}) where {D,R,T}
Returns true if the points ``a`` and ``b`` from the same field are equal,
 and false otherwise.
"""
function ==(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T})::Bool where {D,R,T}
    return a.value==b.value
end

"""
    +(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T}) where {D,R,T}
Returns a new element (of the binary field represented by {D,R,T}) which is the result of ``a+b``.
"""
function +(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    return FieldPoint{D,R,T}(a.value ⊻ b.value)
end

"""
    -(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T}) where {D,R,T}
Returns a new element (of the binary field represented by {D,R,T}) which is the result of ``a-b``.
"""
function -(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    return a+b
end

"""
    -(a::FieldPoint{D,R,T}) where {D,R,T}
Returns a new element (of the binary field represented by {D,R,T}) which is the result of ``-a``.
"""
function -(a::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    return copy(a)
end

#returns the least element b, such that a ≡ b (mod R)
#note: this is the standard algorithm, but faster specialised versions of it are
#available for each of the standard fields (in Field_fastreduce.jl)
function reduce(a::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    #b will should always be such that a ≡ b (mod R)
    #the loop will modify it until it reaches the smallest value that makes that true
    b = copy(a.value)
    r = StaticUInt{ceil(Int,sizeof(R)/sizeof(T)),T}(R)

    #iterate over the excess bits of a, left to right
    for i in (8*sizeof(T)*length(b)-1):-1:D
        if getbit(b, i)==1
            flipbit!(b, i)
            shiftedxor!(b, r, i-D)
        end
    end

    #remove excess blocks from b
    b = changelength(b, ceil(Int,D/(8*sizeof(T))))
    return FieldPoint{D,R,T}(b)
end

#right to left, shift and add
"""
    *(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T}) where {D,R,T}
Returns a new element (of the binary field represented by {D,R,T}) which is the
 result of ``a \\cdot b``.
"""
function *(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    return right_to_left_mult(a, b)
end

#right to left, shift and add
function right_to_left_mult(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    if a.value==b.value return square(a) end

    #c needs to store a polynomial of degree 2D
    c = zero(StaticUInt{ceil(Int,D/(4*sizeof(T))),T})

    for i in 0:(D-1)
        if getbit(a.value,i)==1
            shiftedxor!(c, b.value, i)
        end
    end

    return reduce(FieldPoint{D,R,T}(c))
end

function threads_mult(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    if a.value==b.value return square(a) end

    #cs needs to store polynomials of degree 2D
    cs = [zero(StaticUInt{ceil(Int,D/(4*sizeof(T))),T}) for i=1:Threads.nthreads()]

    Threads.@threads for i in 0:(D-1)
        if getbit(a.value,i)==1
            t = Threads.threadid()
            cs[t] = shiftedxor(cs[t], b.value, i)
        end
    end

    #TODO use something like cumsum / reduce
    for i in 2:Threads.nthreads()
        cs[1] ⊻= cs[i]
    end

    return reduce(FieldPoint{D,R,T}(cs[1]))
end

#still uses shift and add
#but performs reduction itself, rather than calling a reduce function
#Guide to ECC, algorithm 2.33 (ish)
function noreduce_mult(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    if a.value==b.value return square(a) end

    L = ceil(Int,D/(8*sizeof(T)))
    shiftedb::StaticUInt{L,T} = b.value
    c = zero(StaticUInt{L,T})
    r = StaticUInt{ceil(Int,sizeof(R)/sizeof(T)),T}(R)

    for i in 0:(D-1)
        if getbit(a.value,i)==1
            c ⊻= shiftedb
        end
        shiftedb <<= 1
        if getbit(shiftedb, D)==1
            flipbit!(shiftedb, D)
            shiftedb ⊻= r
        end
    end

    return FieldPoint{D,R,T}(c)
end

#Guide to ECC, Algorithm 2.34, right to left comb method
function right_to_left_comb_mult(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    if a.value==b.value return square(a) end

    wordsize = 8*sizeof(T)
    L = ceil(D/wordsize)

    #c needs to store polynomials of degree 2D
    c = zero(StaticUInt{2*L,T})

    #b needs to store polynomials of degree D+wordsize
    bvalue = changelength(b.value, L+1)

    for k in 0:(wordsize-1)
        for j in 0:(L-1)
            if getbit(a.value, wordsize*j + k)==1
                shiftedxor!(c, bvalue, j*wordsize)
            end
        end
        if k!=(wordsize-1) bvalue <<= 1 end
    end

    return reduce(FieldPoint{D,R,T}(c))
end

#Guide to ECC, Algorithm 2.35, left to right comb method
function left_to_right_comb_mult(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    if a.value==b.value return square(a) end

    wordsize = 8*sizeof(T)
    L = ceil(Int,D/wordsize)
    c = zero(StaticUInt{2*L,T})

    for k in (wordsize-1):-1:0
        for j in 0:(L-1)
            if getbit(a.value, wordsize*j + k)==1
                shiftedxor!(c, b.value, j*wordsize)
            end
        end
        if k!=0 c <<= 1 end
    end

    return reduce(FieldPoint{D,R,T}(c))
end

#Guide to ECC, Algorithm 2.36, left to right comb method with windowing
function window_comb_mult(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T}, window::Int)::FieldPoint{D,R,T} where {D,R,T}

    wordsize = 8*sizeof(T)
    L = ceil(Int,D/wordsize)
    Bu::Array{StaticUInt{L+1,T},1} = [small_mult(b, u) for u=0:(1<<window -1)]
    c = zero(StaticUInt{2*L,T})

    for k in ((wordsize÷window)-1):-1:0
        for j in 0:(length(a.value)-1)
            u = getbits(a.value, window*k + wordsize*j, window)
            shiftedxor!(c, Bu[u+1], j*wordsize)
        end
        if k!=0 c <<= window end
    end

    return reduce(FieldPoint{D,R,T}(c))
end

#used for window_comb_mult
function small_mult(a::FieldPoint{D,R,T}, b::Int)::StaticUInt where {D,R,T}
    blen = 8*sizeof(b)
    maxlen = D + blen
    c = zero(StaticUInt{ceil(Int,maxlen/(8*sizeof(T))),T})

    for i in 0:(blen-1)
        if (b>>i)&1==1
            shiftedxor!(c, a.value, i)
        end
    end

    return c
end

#squares the given number, slightly faster than using the standard * algorithm
#adds a zero between every digit of the original
function square(a::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    b = zero(StaticUInt{ceil(Int,D/(4*sizeof(T))),T})
    for i in 0:(D-1)
        if getbit(a.value,i)==1
            flipbit!(b, i*2)
        end
    end

    return reduce(FieldPoint{D,R,T}(b))
end

#uses a version of egcd to invert a
#Algorithm 2.48, Guide to Elliptic Curve Cryptography
"""
    inv(a::FieldPoint{D,R,T}) where {D,R,T}
Returns a new element ``b`` such that ``a b ≡ 1 \\pmod{f_R(x)}``
 (where ``f_R(x)`` is the reduction polynomial for the field).
"""
function inv(a::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    if iszero(a.value) throw(DivideError()) end

    L = ceil(Int,D/(8*sizeof(T)))
    u = a.value
    v = StaticUInt{L,T}(R)
    flipbit!(v, D)
    g1 = one(StaticUInt{L,T})
    g2 = zero(StaticUInt{L,T})

    while !isone(u)
        j = bits(u) - bits(v)
        if j<0
            u, v = v, u
            g1, g2 = g2, g1
            j = -j
        end
        u = shiftedxor(u, v, j)
        g1 = shiftedxor(g1, g2, j)
    end
    return FieldPoint{D,R,T}(g1)
end

"""
    /(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T}) where {D,R,T}
Returns a new element (of the binary field represented by {D,R,T}) which is the
result of ``\\frac{a}{b}``.
"""
function /(a::FieldPoint{D,R,T}, b::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    return a * inv(b)
end

#right to left, square and multiply method
"""
    ^(a::FieldPoint{D,R,T}, b::Integer) where {D,R,T}
Returns a new element (of the binary field represented by {D,R,T}) which is the
result of ``a^b``.
"""
function ^(a::FieldPoint{D,R,T}, b::Integer)::FieldPoint{D,R,T} where {D,R,T}
    c = one(typeof(a))
    squaring = a

    while b>0
        if b & 1 == 1
            c *= squaring
        end
        squaring *= squaring
        b >>>= 1
    end

    return c
end

"""
    random(::Type{FieldPoint{D,R,T}}) where {D,R,T}
Returns a random element of the specified field.
"""
function random(::Type{FieldPoint{D,R,T}})::FieldPoint{D,R,T} where {D,R,T}
    return FieldPoint{D,R,T}(random(StaticUInt{ceil(Int,D/(8*sizeof(T))),T}, D-1))
end

"""
    iszero(a::FieldPoint)
Returns true if ``a`` is the zero element of the field represented by D and R,
 and false otherwise.
"""
function iszero(a::FieldPoint)::Bool
    return iszero(a.value)
end

"""
    zero(::Type{FieldPoint{D,R,T}}) where {D,R,T}
Returns the zero element of the specified field.
"""
function zero(::Type{FieldPoint{D,R,T}})::FieldPoint{D,R,T} where {D,R,T}
    return FieldPoint{D,R,T}(zero(StaticUInt{ceil(Int,D/(8*sizeof(T))),T}))
end

"""
    isone(a::FieldPoint)
Returns true if ``a`` is equal to one, and false otherwise.
"""
function isone(a::FieldPoint)::Bool
    return isone(a.value)
end

"""
    one(::Type{FieldPoint{D,R,T}}) where {D,R,T}
Returns element 1 of the specified field.
"""
function one(::Type{FieldPoint{D,R,T}})::FieldPoint{D,R,T} where {D,R,T}
    return FieldPoint{D,R,T}(one(StaticUInt{ceil(Int,D/(8*sizeof(T))),T}))
end

"""
    sqrt(a::FieldPoint{D,R,T}) where {D,R,T}
Returns ``b`` such that ``b^2 ≡ a \\pmod(R)``.
"""
function sqrt(a::FieldPoint{D,R,T})::FieldPoint{D,R,T} where {D,R,T}
    #a^{2^{D-1}}
    for i in 1:(D-1)
        a *= a
    end
    return a
end

"""
    convert(::Type{BigInt}, a::FieldPoint)
Converts the given field point to a number (of type BigInt), following the procedure
 set out in SEC 1 (version 2) 2.3.9.
"""
function convert(::Type{BigInt}, a::FieldPoint)::BigInt
    return convert(BigInt, a.value)
end

#sec2 v2 (and v1), table 3:
FieldPoint113 = FieldPoint{113, UInt16(512+1),UInt64} #v1 only
FieldPoint131 = FieldPoint{131, UInt16(256+8+4+1),UInt64} #v1 only
FieldPoint163 = FieldPoint{163, UInt16(128+64+8+1),UInt64}
FieldPoint193 = FieldPoint{193, (UInt16(1)<<15) + UInt16(1),UInt64} #v1 only
FieldPoint233 = FieldPoint{233, (UInt128(1)<<74) + UInt128(1),UInt64}
FieldPoint239 = FieldPoint{239, (UInt64(1)<<36) + UInt64(1),UInt64}
FieldPoint283 = FieldPoint{283, (UInt16(1)<<12) + UInt16(128+32+1),UInt64}
FieldPoint409 = FieldPoint{409, (UInt128(1)<<87) + UInt128(1),UInt64}
FieldPoint571 = FieldPoint{571, (UInt16(1)<<10) + UInt16(32+4+1),UInt64}
