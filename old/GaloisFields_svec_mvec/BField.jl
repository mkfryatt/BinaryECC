#D is the degree of the reduction polynomial
#R is the reduction polynomial without the x^D term
"""
    BFieldElt{D,R,T}
Represents a point in the binary field which has order ``2^D`` and reduction polynomial

``x^D + x^{r_n} + \\cdots + x^{r_0}``

where ``R = r_n r_{n-1}\\ldots r_1 r_0`` in binary. The polynomial is stored
with an array of T words.

Types for points in the standard fields (taken from SEC 2, table 3)
 are available:
- BFieldElt113
- BFieldElt131
- BFieldElt163
- BFieldElt193
- BFieldElt233
- BFieldElt239
- BFieldElt283
- BFieldElt409
- BFieldElt571
"""
struct BFieldElt{D,R,T<:Unsigned}
    value::SVector{L,T} where L

    BFieldElt{D,R,T}(value::Integer) where {D,R,T} = new{D,R,T}(convert(SVector{ceil(Int,D/bitsize(T)),T}, value))

    BFieldElt{D,R,T}(value::SVector{L,T}) where {D,R,T,L} = new(value)

    BFieldElt{D,R,T}(value::MVector{L,T}) where {D,R,T,L} = new(SVector{L,T}(value))
end

"""
    BFieldElt{D,R,T}(s::String) where {D,R,T}
Using the procedure set out in SEC 1 (version 2) 2.3.6,
this converts a hex string to a field element.
"""
function BFieldElt{D,R,T}(s::String) where {D,R,T}
    s = replace(s, " " => "")
    if length(s)!=ceil(D / 8)*2 throw(ArgumentError("Octet string is of the incorrect length for this field.")) end
    value = SVector{ceil(Int,D/bitsize(T)),T}(s)
    return BFieldElt{D,R,T}(value)
end

function copy(x::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    return x
end

"""
    ==(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}) where {D,R,T}
Returns true if the points ``a`` and ``b`` from the same field are equal,
 and false otherwise.
"""
function ==(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T})::Bool where {D,R,T}
    return a.value==b.value
end

"""
    +(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}) where {D,R,T}
Returns a new element (of the binary field represented by {D,R}) which is the result of ``a+b``.
"""
function +(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    return BFieldElt{D,R,T}(a.value ⊻ b.value)
end

"""
    -(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}) where {D,R,T}
Returns a new element (of the binary field represented by {D,R}) which is the result of ``a-b``.
"""
function -(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    return a+b
end

function -(a::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    return a
end

#note: this is the standard algorithm, but faster specialised versions of it are
#available for each of the standard fields (in Field_fastreduce.jl)
"""
    reduce(a::BFieldElt{D,R,T}) where {D,R,T}
Returns the least element ``b``, such that ``a \\equiv b \\pmod{R}``.
"""
function reduce(a::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    #b will should always be such that a ≡ b (mod R)
    #the loop will modify it until it reaches the smallest value that makes that true
    b = MVector(a.value)
    r = convert(SVector{128÷bitsize(T),T}, R)

    #iterate over the excess bits of a, left to right
    for i in (length(b)*bitsize(T)-1):-1:D
        if getbit(b, i)==1
            flipbit!(b, i)
            shiftedxor!(b, r, i-D)
        end
    end

    #remove excess blocks from b
    b = changelength(b, ceil(Int,D/bitsize(T)))
    return BFieldElt{D,R,T}(b)
end

"""
    *(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}) where {D,R,T}
Returns a new element (of the binary field represented by {D,R}) which is the
 result of ``a \\cdot b``.
"""
function *(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    return mult_shiftandadd_window(a, b, 4)
end

"""
    mult_shiftandadd(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}) where {D,R,T}
Returns ``a \\cdot b`` using the right to left shift and add method.
"""
function mult_shiftandadd(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    #c needs to store a polynomial of degree 2D
    c = zeros(MVector{ceil(Int,2*D/bitsize(T)),T})

    for i in 0:(D-1)
        if getbit(a.value,i)==1
            shiftedxor!(c, b.value, i)
        end
    end

    return reduce(BFieldElt{D,R,T}(c))
end

function mult_shiftandadd_window(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}, window::Int)::BFieldElt{D,R,T} where {D,R,T}
    Bu = [small_mult(b, u) for u::UInt8=UInt8(1):(UInt8(1)<<window -UInt8(1))]
    c = zero(MVector{ceil(Int,2*D/bitsize(T)),T})

    for i in 0:window:D-1
        u = getbits(a.value, i, window)
        if u!=0 shiftedxor!(c, Bu[u], i) end
    end

    return reduce(BFieldElt{D,R,T}(c))
end

"""
    mult_threaded(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}) where {D,R,T}
Returns ``a \\cdot b`` using the right to left shift and add method with multithreading.
"""
function mult_threaded(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    mid = (length(a.value) ÷ 2)*64 -1
    c1 = Threads.@spawn mult_threaded_helper($a, $b, 0, $mid)
    c2 = Threads.@spawn mult_threaded_helper($a, $b, $mid+1, $D-1)
    return reduce(BFieldElt{D,R,T}(fetch(c1) ⊻ fetch(c2)))
end
function mult_threaded_helper(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}, start::Int, stop::Int)::MVector{ceil(Int,2*D/bitsize(T)),T} where {D,R,T}
    c = zeros(MVector{ceil(Int,2*D/bitsize(T)),T})
    for i in start:stop
        if getbit(a.value,i)==1
            shiftedxor!(c, b.value, i)
        end
    end
    return c
end

function mult_threaded_window(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}, w::Int)::BFieldElt{D,R,T} where {D,R,T}
    mid = (length(a.value) ÷ 2)*64 -1
    c1 = Threads.@spawn mult_window_helper($a, $b, $w, 0, $mid)
    c2 = Threads.@spawn mult_window_helper($a, $b, $w, $mid+1, $D-1)
    return reduce(BFieldElt{D,R,T}(fetch(c1) ⊻ fetch(c2)))
end
function mult_window_helper(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}, w::Int,
     start::Int, stop::Int)::MVector{ceil(Int,2*D/bitsize(T)),T} where {D,R,T}

    Bu::Array{MVector{L,T},1} where L = [small_mult(b, u) for u::UInt8=UInt8(1):(UInt8(1)<<w -UInt8(1))]
    c = zeros(MVector{ceil(Int,2*D/bitsize(T)),T})

    extra = (stop+1-start) % w
    for i in start:w:(stop-extra)
        u = getbits(a.value, i, w)
        if u!=0 shiftedxor!(c, Bu[u], i) end
    end

    i = stop - extra +1
    u = getbits(a.value, i, extra)
    if u!=0 shiftedxor!(c, Bu[u], i) end

    return c
end

"""
    mult_ownreduce(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}) where {D,R,T}
Returns ``a \\cdot b`` using the right to left shift and add method,
without needing to call a reduction function.
"""
function mult_ownreduce(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    L = ceil(Int,D/bitsize(T))
    shiftedb = MVector(b.value)
    c = zeros(MVector{L,T})
    r = convert(SVector{128÷bitsize(T),T}, R)

    for i in 0:(D-1)
        if getbit(a.value,i)==1
            xor!(c, shiftedb)
        end
        leftshift!(shiftedb, 1)
        if getbit(shiftedb, D)==1
            flipbit!(shiftedb, D)
            xor!(shiftedb, r)
        end
    end

    return BFieldElt{D,R,T}(c)
end

"""
    mult_comb_rtl(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}) where {D,R,T}
Returns ``a \\cdot b`` using a right to left comb method
(described in Guide to Elliptic Curve Cryptography, algorithm 2.34).
"""
function mult_comb_rtl(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    L = ceil(Int,D/bitsize(T))
    c = zeros(MVector{2*L,T})

    #b needs to store polynomials of degree D+wordsize
    bvalue = changelength(MVector(b.value), L+1)

    for k in 0:(bitsize(T)-1)
        for j in 0:(L-1)
            if getbit(a.value, j*bitsize(T) + k)==1
                shiftedxor!(c, bvalue, j*bitsize(T))
            end
        end
        if k!=(bitsize(T)-1) leftshift!(bvalue,1) end
    end

    return reduce(BFieldElt{D,R,T}(c))
end

"""
    mult_comb_ltr(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}) where {D,R,T}
Returns ``a \\cdot b`` using a left to right comb method
(described in Guide to Elliptic Curve Cryptography, algorithm 2.35).
"""
function mult_comb_ltr(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    L = ceil(Int,D/bitsize(T))
    c = zeros(MVector{2*L,T})

    for k in (bitsize(T)-1):-1:0
        for j in 0:(L-1)
            if getbit(a.value, bitsize(T)*j + k)==1
                shiftedxor!(c, b.value, j*bitsize(T))
            end
        end
        if k!=0 leftshift!(c,1) end
    end

    return reduce(BFieldElt{D,R,T}(c))
end

"""
    mult_comb_window(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}, window::Int) where {D,R,T}
Returns ``a \\cdot b`` using a left to right comb method windowing
(described in Guide to Elliptic Curve Cryptography, algorithm 2.36).

Performs best with a window size of 4.
"""
function mult_comb_window(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}, window::Int)::BFieldElt{D,R,T} where {D,R,T}
    L = ceil(Int,D/bitsize(T))
    Bu = [small_mult(b, u) for u=0:(1<<window -1)]
    c = zeros(MVector{2*L,T})

    for k in ((bitsize(T)÷window)-1):-1:0
        for j in 0:(length(a.value)-1)
            u = getbits(a.value, window*k + bitsize(T)*j, window)
            shiftedxor!(c, Bu[u+1], j*bitsize(T))
        end
        if k!=0
            leftshift!(c, window)
        end
    end

    return reduce(BFieldElt{D,R,T}(c))
end

function small_mult(a::BFieldElt{D,R,T}, b::UInt8)::MVector{((D+8) ÷ bitsize(T)) +1,T} where {D,R,T}
    blen = 8
    newL = ((D+blen) ÷ bitsize(T)) +1
    c = zeros(MVector{newL,T})

    for i in 0:(blen-1)
        if (b>>i)&1==1
            shiftedxor!(c, a.value, i)
        end
    end

    return c
end
function small_mult(a::MVector{L,T}, b::UInt8)::MVector{L+1,T} where {L,T}
    c = zero(MVector{L+1,T})
    for i in 0:(L*bitsize(T))
        if (b>>i)&1==1
            shiftedxor!(c, a, i)
        end
    end
    return c
end

"""
    square(a::BFieldElt{D,R,T}) where {D,R,T}
Returns a new element (of the binary field represented by {D,R}) which is the
result of ``a^2``.
"""
function square(a::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    return square_window(a, 4)
end

#adds a zero between every digit of the original
function square_standard(a::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    b = zeros(MVector{ceil(Int,2*D/bitsize(T)),T})
    for i in 0:(D-1)
        if getbit(a.value,i)==1
            flipbit!(b, i*2)
        end
    end

    return reduce(BFieldElt{D,R,T}(b))
end

"""
    square_window(a::BFieldElt{D,R,T}, window::Int) where {D,R,T}
Returns ``a^2`` by inserting a zero between every bit in the original, using
the specified window size.
"""
function square_window(a::BFieldElt{D,R,T}, window::Int)::BFieldElt{D,R,T} where {D,R,T}
    b = zeros(MVector{ceil(Int,2*D/bitsize(T)),T})
    spread = [SVector{1,T}(spread_bits(i)) for i=1:(1<<window -1)]

    for i in 0:window:D-1
        u = getbits(a.value, i, window)
        if u!=0 shiftedxor!(b, spread[u], i*2) end
    end

    return reduce(BFieldElt{D,R,T}(b))
end

#needed for window_square
function spread_bits(a::Int)::Int
    b = 0
    for i in 0:bits(a)
        if (a>>i)&1==1
            b ⊻= 1<<(2*i)
        end
    end
    return b
end

#uses a version of egcd to invert a
#Algorithm 2.48, Guide to Elliptic Curve Cryptography
"""
    inv(a::BFieldElt{D,R,T}) where {D,R,T}
Returns a new element ``b`` such that ``a b ≡ 1 \\pmod{f_R(x)}``
 (where ``f_R(x)`` is the reduction polynomial for the field).
"""
function inv(a::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    if iszero(a.value) throw(DivideError()) end

    L = ceil(Int,D/bitsize(T))
    u = MVector(a.value)
    v = convert(MVector{L,T}, R)
    flipbit!(v, D)
    g1 = one(MVector{L,T})
    g2 = zero(MVector{L,T})

    while !isone(u)
        j = bits(u) - bits(v)
        if j<0
            u, v = v, u
            g1, g2 = g2, g1
            j = -j
        end
        shiftedxor!(u, v, j)
        shiftedxor!(g1, g2, j)
    end
    return BFieldElt{D,R,T}(g1)
end

"""
    /(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T}) where {D,R,T}
Returns a new element (of the binary field represented by {D,R}) which is the
result of ``\\frac{a}{b}``.
"""
function /(a::BFieldElt{D,R,T}, b::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    return a * inv(b)
end

function div_threaded(b::BFieldElt{D,R,T}, a::BFieldElt{D,R,T}, window::Int)::BFieldElt{D,R,T} where {D,R,T}
    a_inv_task = Threads.@spawn inv($a)
    Bu::Array{MVector{((D+8) ÷ bitsize(T)) +1,T},1} = [small_mult(b, u) for u::UInt8=UInt8(1):(UInt8(1)<<window -UInt8(1))]

    #c needs to store a polynomial of degree 2D
    c = zeros(MVector{ceil(Int,2*D/bitsize(T)),T})
    a = fetch(a_inv_task)

    for i in 0:window:D-1
        u = getbits(a.value, i, window)
        if u!=0 shiftedxor!(c, Bu[u], i) end
    end

    return reduce(BFieldElt{D,R,T}(c))
end

#right to left, square and multiply method
"""
    ^(a::BFieldElt{D,R,T}, b::Integer) where {D,R,T}
Returns a new element (of the binary field represented by {D,R}) which is the
result of ``a^b``.
"""
function ^(a::BFieldElt{D,R,T}, b::Integer)::BFieldElt{D,R,T} where {D,R,T}
    c = one(BFieldElt{D,R,T})
    squaring = a

    while b>0
        if b & 1 == 1
            c *= squaring
        end
        squaring = square(squaring)
        b >>>= 1
    end

    return c
end

"""
    random(::Type{BFieldElt{D,R,T}}) where {D,R,T}
Returns a random element of the specified field.
"""
function random(::Type{BFieldElt{D,R,T}})::BFieldElt{D,R,T} where {D,R,T}
    return BFieldElt{D,R,T}(random(SVector{ceil(Int,D/bitsize(T)),T}, D-1))
end

"""
    iszero(a::BFieldElt)
Returns true if ``a`` is the zero element of the field represented by D and R,
 and false otherwise.
"""
function iszero(a::BFieldElt)::Bool
    return iszero(a.value)
end

"""
    zero(::Type{BFieldElt{D,R,T}}) where {D,R,T}
Returns the zero element of the specified field.
"""
function zero(::Type{BFieldElt{D,R,T}})::BFieldElt{D,R,T} where {D,R,T}
    return BFieldElt{D,R,T}(zero(SVector{ceil(Int,D/bitsize(T)),T}))
end

"""
    isone(a::BFieldElt)
Returns true if ``a`` is equal to one, and false otherwise.
"""
function isone(a::BFieldElt)::Bool
    return isone(a.value)
end

"""
    one(::Type{BFieldElt{D,R,T}}) where {D,R,T}
Returns element 1 of the specified field.
"""
function one(::Type{BFieldElt{D,R,T}})::BFieldElt{D,R,T} where {D,R,T}
    return BFieldElt{D,R,T}(one(SVector{ceil(Int,D/bitsize(T)),T}))
end

"""
    sqrt(a::BFieldElt{D,R,T}) where {D,R,T}
Returns ``b`` such that ``b^2 ≡ a \\pmod(R)``.
"""
function sqrt(a::BFieldElt{D,R,T})::BFieldElt{D,R,T} where {D,R,T}
    #a^{2^{D-1}}
    for i in 1:(D-1)
        a = square(a)
    end
    return a
end

"""
    convert(::Type{BigInt}, a::BFieldElt)
Converts the given field point to a number (of type BigInt), following the procedure
 set out in SEC 1 (version 2) 2.3.9.
"""
function convert(::Type{BigInt}, a::BFieldElt)::BigInt
    return convert(BigInt, a.value)
end

#sec2 v2 (and v1), table 3:
BFieldElt113{T<:Unsigned} = BFieldElt{113, UInt16(512+1),T}#v1 only
BFieldElt131{T<:Unsigned} = BFieldElt{131, UInt16(256+8+4+1),T} #v1 only
BFieldElt163{T<:Unsigned} = BFieldElt{163, UInt16(128+64+8+1),T}
BFieldElt193{T<:Unsigned} = BFieldElt{193, (UInt16(1)<<15) + UInt16(1),T} #v1 only
BFieldElt233{T<:Unsigned} = BFieldElt{233, (UInt128(1)<<74) + UInt128(1),T}
BFieldElt239{T<:Unsigned} = BFieldElt{239, (UInt64(1)<<36) + UInt64(1),T}
BFieldElt283{T<:Unsigned} = BFieldElt{283, (UInt16(1)<<12) + UInt16(128+32+1),T}
BFieldElt409{T<:Unsigned} = BFieldElt{409, (UInt128(1)<<87) + UInt128(1),T}
BFieldElt571{T<:Unsigned} = BFieldElt{571, (UInt16(1)<<10) + UInt16(32+4+1),T}
