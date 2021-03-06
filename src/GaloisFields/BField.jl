#D is the degree of the reduction polynomial
#R is the reduction polynomial without the x^D term
"""
    BFieldElt{D,R,T,L}
Represents an element in the binary field which has order ``2^D`` and reduction polynomial

``x^D + x^{r_n} + \\cdots + x^{r_0}``

where ``R = r_n r_{n-1}\\ldots r_1 r_0`` in binary. The polynomial is stored
with a length `L` array of type `T<:Unsigned` words, using a `StaticUInt{L,T}` object (a wrapper
for the type `MVector{L,T}` from StaticArrays).

Note: binary field arithmetic has been tested with words of type `UInt8`, `UInt16`,
`UInt32`, `UInt64` and `UInt128`. For any other possible word types, it is advisable to perform
additional testing.

Types for points in the standard fields (taken from SEC 2, table 3)
 are available:
- `BFieldElt113{T,L}`
- `BFieldElt131{T,L}`
- `BFieldElt163{T,L}`
- `BFieldElt193{T,L}`
- `BFieldElt233{T,L}`
- `BFieldElt239{T,L}`
- `BFieldElt283{T,L}`
- `BFieldElt409{T,L}`
- `BFieldElt571{T,L}`

Each of these fields can be easily created by the functions `B113(T)`, `B131(T)`, etc.,
which take a word type `T` and they return a binary field element type with the smallest
value of `L`. The function `B(D, R, T)` returns similar types, but for custom fields
defined for `D` and `R`.
"""
struct BFieldElt{D,R,T<:Unsigned,L}
    value::StaticUInt{L,T}
    BFieldElt{D,R,T,L}(value::Integer) where {D,R,T,L} =
        new(StaticUInt{L,T}(value))
    BFieldElt{D,R,T,L}(value::StaticUInt{L,T}) where {D,R,T,L} = new(value)
end

#ceil(Int,D/bitsize(T))
"""
    BFieldElt{D,R,T,L}(s::String) where {D,R,T,L}
Using the procedure set out in SEC 1 (version 2) 2.3.6,
this converts a hex string to a field element.
"""
function BFieldElt{D,R,T,L}(s::String) where {D,R,T,L}
    s = replace(s, " " => "")
    if length(s)!=ceil(D / 8)*2 throw(ArgumentError("Octet string is of the incorrect length for this field.")) end
    value = StaticUInt{ceil(Int,D/bitsize(T)),T}(s)
    return BFieldElt{D,R,T,L}(value)
end

function repr(a::BFieldElt{D,R,T,L}) where {D,R,T,L}
    return repr(a.value)
end

function copy(x::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return BFieldElt{D,R,T,L}(copy(x.value))
end

"""
    ==(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}
Returns true if `a` and `b` represent the same field element,
 and false otherwise.
"""
function ==(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L})::Bool where {D,R,T,L}
    return a.value==b.value
end

"""
    +(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}
Returns a new element which is the result of ``a+b``.
"""
function +(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return BFieldElt{D,R,T,L}(a.value ⊻ b.value)
end

"""
    -(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}
Returns a new element which is the result of ``a-b``.
"""
function -(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return a+b
end

function -(a::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return copy(a)
end

#note: this is the standard algorithm, but faster specialised versions of it are
#available for each of the standard fields (in Field_fastreduce.jl)
"""
    reduce(a::BFieldElt{D,R,T,L}) where {D,R,T,L}
Returns the least element `b`, such that ``a \\equiv b`` in the field represented
by `D` and `R`.
"""
function reduce(a::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,ceil(Int,D/bitsize(T))} where {D,R,T,L}
    #b will should always be such that a ≡ b (mod R)
    #the loop will modify it until it reaches the smallest value that makes that true
    b = copy(a.value)
    r = StaticUInt{128÷bitsize(T),T}(R)

    #iterate over the excess bits of a, left to right
    for i in (2*D):-1:D
        if getbit(b, i)==1
            flipbit!(b, i)
            shiftedxor!(b, r, i-D)
        end
    end

    #remove excess blocks from b
    b = changelength(b, ceil(Int,D/bitsize(T)))
    return BFieldElt{D,R,T,ceil(Int,D/bitsize(T))}(b)
end

"""
    *(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}
Returns a new element which is the result of ``a \\cdot b``. This is the default
multiplication routine for binary field arithmetic, chosen to have high performance.
"""
function *(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return mult_shiftandadd_window(a, b, 4)
end

"""
    mult_shiftandadd(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}
Binary field multiplication using a right-to-left shift-and-add method.
"""
function mult_shiftandadd(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    longL = ceil(Int,2*D/bitsize(T))
    c = zero(StaticUInt{longL,T})

    for i in 0:(D-1)
        if getbit(a.value,i)==1
            shiftedxor!(c, b.value, i)
        end
    end

    return reduce(BFieldElt{D,R,T,longL}(c))
end

"""
    mult_shiftandadd_window(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}, window::Int) where {D,R,T,L}
Binary field multiplication using a right-to-left shift-and-add method with a window
size of `window`. The optimal window size for this routine is 4.
"""
function mult_shiftandadd_window(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}, window::Int)::BFieldElt{D,R,T,L} where {D,R,T,L}
    Bu = [small_mult(b, u) for u::UInt8=UInt8(1):(UInt8(1)<<window -UInt8(1))]
    longL = ceil(Int,2*D/bitsize(T))
    c = zero(StaticUInt{longL,T})

    for i in 0:window:D-1
        u = getbits(a.value, i, window)
        if u!=0 shiftedxor!(c, Bu[u], i) end
    end

    return reduce(BFieldElt{D,R,T,longL}(c))
end

"""
    mult_threaded(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}
Binary field multiplication using a right-to-left shift-and-add method, by spawning
an additional thread.
"""
function mult_threaded(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    mid = (L ÷ 2)*bitsize(T) -1
    c1 = Threads.@spawn mult_threaded_helper($a, $b, 0, $mid)
    c2 = mult_threaded_helper(a, b, mid+1, D-1)
    longL = ceil(Int,2*D/bitsize(T))
    xor!(c2, fetch(c1))
    return reduce(BFieldElt{D,R,T,longL}(c2))
end
function mult_threaded_helper(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L},
    start::Int, stop::Int)::StaticUInt{ceil(Int,2*D/bitsize(T)),T} where {D,R,T,L}
    longL = ceil(Int,2*D/bitsize(T))
    c = zero(StaticUInt{longL,T})

    for i in start:stop
        if getbit(a.value,i)==1
            shiftedxor!(c, b.value, i)
        end
    end

    return c
end

"""
    mult_threaded_window(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}, w::Int) where {D,R,T,L}
Binary field multiplication using a windowed right-to-left shift-and-add method, by spawning
an additional thread.
"""
function mult_threaded_window(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}, w::Int)::BFieldElt{D,R,T,L} where {D,R,T,L}
    mid = (length(a.value) ÷ 2)*64 -1
    c1 = Threads.@spawn mult_window_helper($a, $b, $w, 0, $mid)
    c2 = Threads.@spawn mult_window_helper($a, $b, $w, $mid+1, $D-1)
    longL = ceil(Int,2*D/bitsize(T))
    return reduce(BFieldElt{D,R,T,longL}(fetch(c1) ⊻ fetch(c2)))
end
function mult_window_helper(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}, w::Int,
     start::Int, stop::Int)::StaticUInt{ceil(Int,2*D/bitsize(T)),T} where {D,R,T,L}
    Bu = [small_mult(b, u) for u::UInt8=UInt8(1):(UInt8(1)<<w -UInt8(1))]
    longL = ceil(Int,2*D/bitsize(T))
    c = zero(StaticUInt{longL,T})

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
    mult_ownreduce(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}
Binary field multiplication using a windowed right-to-left shift-and-add method,
in which reduction is performed alongside multiplication.
"""
function mult_ownreduce(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    shiftedb::StaticUInt{L,T} = copy(b.value)
    c = zero(StaticUInt{L,T})
    r = StaticUInt{128÷bitsize(T),T}(R)

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

    return BFieldElt{D,R,T,L}(c)
end

"""
    mult_comb_rtl(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}
Binary field multiplication using a right-to-left comb method.
"""
function mult_comb_rtl(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    c = zero(StaticUInt{2*L,T})
    bvalue = changelength(b.value, L+1)

    for k in 0:(bitsize(T)-1)
        for j in 0:(L-1)
            if getbit(a.value, j*bitsize(T) + k)==1
                shiftedxor!(c, bvalue, j*bitsize(T))
            end
        end
        if k!=(bitsize(T)-1) leftshift!(bvalue,1) end
    end

    return reduce(BFieldElt{D,R,T,2*L}(c))
end

"""
    mult_comb_ltr(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}
Binary field multiplication using a left-to-right comb method.
"""
function mult_comb_ltr(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    c = zero(StaticUInt{2*L,T})

    for k in (bitsize(T)-1):-1:0
        for j in 0:(L-1)
            if getbit(a.value, bitsize(T)*j + k)==1
                shiftedxor!(c, b.value, j*bitsize(T))
            end
        end
        if k!=0 leftshift!(c,1) end
    end

    return reduce(BFieldElt{D,R,T,2*L}(c))
end

"""
    mult_comb_window(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}, window::Int) where {D,R,T,L}
Binary field multiplication using a windowed left-to-right comb method.
Performs best with a window size of 4.
"""
function mult_comb_window(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}, window::Int)::BFieldElt{D,R,T,L} where {D,R,T,L}
    Bu = [small_mult(b, u) for u=UInt8(0):UInt8(1<<window -1)]
    c = zero(StaticUInt{2*L,T})

    for k in ((bitsize(T)÷window)-1):-1:0
        for j in 0:(length(a.value)-1)
            u = getbits(a.value, window*k + bitsize(T)*j, window)
            shiftedxor!(c, Bu[u+1], j*bitsize(T))
        end
        if k!=0
            leftshift!(c, window)
        end
    end

    return reduce(BFieldElt{D,R,T,2*L}(c))
end

function small_mult(a::BFieldElt{D,R,T,L}, b::UInt8)::StaticUInt{((D+8) ÷ bitsize(T)) +1,T} where {D,R,T,L}
    blen = 8
    newL = ((D+blen) ÷ bitsize(T)) +1
    c = zero(StaticUInt{newL,T})

    for i in 0:(blen-1)
        if (b>>i)&1==1
            shiftedxor!(c, a.value, i)
        end
    end

    return c
end

function small_mult(a::StaticUInt{L,T}, b::UInt8)::StaticUInt{L+1,T} where {L,T}
    c = zero(StaticUInt{L+1,T})
    for i in 0:(L*bitsize(T))
        if (b>>i)&1==1
            shiftedxor!(c, a, i)
        end
    end
    return c
end

"""
    square(a::BFieldElt{D,R,T,L}) where {D,R,T,L}
Returns a new element which is the result of ``a^2``, using the default routine for
high performance.
"""
function square(a::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return square_window(a, 4)
end

#adds a zero between every digit of the original
"""
    square_standard(a::BFieldElt{D,R,T,L}) where {D,R,T,L}
Binary field squaring performed by shifting each bit ``b_i`` left by ``i``.
"""
function square_standard(a::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    longL = ceil(Int,2*D/bitsize(T))
    b = zero(StaticUInt{longL,T})
    for i in 0:(D-1)
        if getbit(a.value,i)==1
            flipbit!(b, i*2)
        end
    end

    return reduce(BFieldElt{D,R,T,longL}(b))
end

"""
    square_window(a::BFieldElt{D,R,T,L}, window::Int) where {D,R,T,L}
Binary field squaring performed with a windowed method, in which the square of each
size `window` block is calculated upfront. Performs best with a window size of 4.
"""
function square_window(a::BFieldElt{D,R,T,L}, window::Int)::BFieldElt{D,R,T,L} where {D,R,T,L}
    longL = ceil(Int,2*D/bitsize(T))
    b = zero(StaticUInt{longL,T})
    spread = [StaticUInt{1,T}(spread_bits(i)) for i=0:(1<<window -1)]

    for i in 0:window:D-1
        u = getbits(a.value, i, window)
        shiftedxor!(b, spread[u+1], i*2)
    end

    return reduce(BFieldElt{D,R,T,longL}(b))
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
    inv(a::BFieldElt{D,R,T,L}) where {D,R,T,L}
Returns a new element `b` such that ``a b ≡ 1`` in the field represented by `D`
and `R`.
"""
function inv(a::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    if iszero(a.value) throw(DivideError()) end

    u = copy(a.value)
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
        shiftedxor!(u, v, j)
        shiftedxor!(g1, g2, j)
    end
    return BFieldElt{D,R,T,L}(g1)
end

"""
    /(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}
Returns a new element which is the result of ``\\frac{a}{b}``.
"""
function /(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return a * inv(b)
end

function div_threaded(b::BFieldElt{D,R,T,L}, a::BFieldElt{D,R,T,L}, window::Int)::BFieldElt{D,R,T,L} where {D,R,T,L}
    a_inv_task = Threads.@spawn inv($a)
    Bu = [small_mult(b, u) for u::UInt8=UInt8(1):(UInt8(1)<<window -UInt8(1))]

    #c needs to store a polynomial of degree 2D
    longL = ceil(Int,2*D/bitsize(T))
    c = zero(StaticUInt{longL,T})
    a = fetch(a_inv_task)

    for i in 0:window:D-1
        u = getbits(a.value, i, window)
        if u!=0 shiftedxor!(c, Bu[u], i) end
    end

    return reduce(BFieldElt{D,R,T,longL}(c))
end

#right to left, square and multiply method
"""
    ^(a::BFieldElt{D,R,T,L}, b::Integer) where {D,R,T,L}
Returns a new element which is the result of ``a^b``. If squaring is required, i.e.
`b==2`, it is faster to call `square(a)` directly.
"""
function ^(a::BFieldElt{D,R,T,L}, b::Integer)::BFieldElt{D,R,T,L} where {D,R,T,L}
    c = one(BFieldElt{D,R,T,L})
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
    random(::Type{BFieldElt{D,R,T,L}}) where {D,R,T,L}
Returns a random element of the specified field.
"""
function random(::Type{BFieldElt{D,R,T,L}})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return BFieldElt{D,R,T,L}(random(StaticUInt{L,T}, D-1))
end

"""
    iszero(a::BFieldElt)
Returns true if ``a`` is the zero element of the binary field represented by `D` and `R`,
 and false otherwise.
"""
function iszero(a::BFieldElt)::Bool
    return iszero(a.value)
end

"""
    zero(::Type{BFieldElt{D,R,T,L}}) where {D,R,T,L}
Returns the zero element (additive identity) of the specified field.
"""
function zero(::Type{BFieldElt{D,R,T,L}})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return BFieldElt{D,R,T,L}(zero(StaticUInt{L,T}))
end

function zero(::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return BFieldElt{D,R,T,L}(zero(StaticUInt{L,T}))
end

"""
    isone(a::BFieldElt)
Returns true if ``a`` is equal to one (multiplicative identity), and false otherwise.
"""
function isone(a::BFieldElt)::Bool
    return isone(a.value)
end

"""
    one(::Type{BFieldElt{D,R,T,L}}) where {D,R,T,L}
Returns element one (multiplicative identity) of the specified field.
"""
function one(::Type{BFieldElt{D,R,T,L}})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return BFieldElt{D,R,T,L}(one(StaticUInt{L,T}))
end

function one(::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    return BFieldElt{D,R,T,L}(one(StaticUInt{L,T}))
end

"""
    sqrt(a::BFieldElt{D,R,T,L}) where {D,R,T,L}
Returns `b` such that ``b^2 \\equiv a``.
"""
function sqrt(a::BFieldElt{D,R,T,L})::BFieldElt{D,R,T,L} where {D,R,T,L}
    #a^{2^{D-1}}
    for i in 1:(D-1)
        a *= a
    end
    return a
end

"""
    convert(::Type{BigInt}, a::BFieldElt)
Converts the given field point to a number (of type `BigInt`), following the procedure
 set out in SEC 1 (version 2) 2.3.9.
"""
function convert(::Type{BigInt}, a::BFieldElt)::BigInt
    return convert(BigInt, a.value)
end

#sec2 v2 (and v1), table 3:
BFieldElt113{T<:Unsigned,L} = BFieldElt{113, UInt16(512+1),T,L}#v1 only
BFieldElt131{T<:Unsigned,L} = BFieldElt{131, UInt16(256+8+4+1),T,L} #v1 only
BFieldElt163{T<:Unsigned,L} = BFieldElt{163, UInt16(128+64+8+1),T,L}
BFieldElt193{T<:Unsigned,L} = BFieldElt{193, (UInt16(1)<<15) + UInt16(1),T,L} #v1 only
BFieldElt233{T<:Unsigned,L} = BFieldElt{233, (UInt128(1)<<74) + UInt128(1),T,L}
BFieldElt239{T<:Unsigned,L} = BFieldElt{239, (UInt64(1)<<36) + UInt64(1),T,L}
BFieldElt283{T<:Unsigned,L} = BFieldElt{283, (UInt16(1)<<12) + UInt16(128+32+1),T,L}
BFieldElt409{T<:Unsigned,L} = BFieldElt{409, (UInt128(1)<<87) + UInt128(1),T,L}
BFieldElt571{T<:Unsigned,L} = BFieldElt{571, (UInt16(1)<<10) + UInt16(32+4+1),T,L}

B(D::Integer, R::Integer, T::Type{U}) where U<:Unsigned = BFieldElt{D, R, T, ceil(Int,D/bitsize(T))}
B113(T) = B(113, UInt16(512+1),T)
B131(T) = B(131, UInt16(256+8+4+1),T)
B163(T) = B(163, UInt16(128+64+8+1),T)
B193(T) = B(193, (UInt16(1)<<15) + UInt16(1),T)
B233(T) = B(233, (UInt128(1)<<74) + UInt128(1),T)
B239(T) = B(239, (UInt64(1)<<36) + UInt64(1),T)
B283(T) = B(283, (UInt16(1)<<12) + UInt16(128+32+1),T)
B409(T) = B(409, (UInt128(1)<<87) + UInt128(1),T)
B571(T) = B(571, (UInt16(1)<<10) + UInt16(32+4+1),T)

"""
    @fastreduce(D,R)
A macro to produce a specialised reduction function for the binary field denoted by
``D`` and ``R``. It is strongly recommended this is run for any new user-defined fields,
as it achieves signficantly higher performance than the generic reduction routine.
"""
macro fastreduce(D,R)
    D_val = eval(D)
    R_val = eval(R)

    indices = []
    i = 0
    while R_val>0
        if R_val%2==1
            prepend!(indices, [i])
        end
        i += 1
        R_val >>= 1
    end
    prepend!(indices, [D_val])

    text = """
    function reduce(a::BFieldElt{$D,$R,T,L})::BFieldElt{$D,$R,T,ceil(Int,$D/bitsize(T))} where {T,L}
        b::StaticUInt{L,T} = copy(a.value)
        lastword = ceil(Int,$D/bitsize(T))
        for i in L:-1:(1+lastword)
            t = StaticUInt{1,T}([b.value[i]])
    """

    for i in indices
        text *= "shiftedxor!(b, t, bitsize(T)*(i-1)-$D+$i)\n"
    end

    text *= """end
            extra = $D % bitsize(T)
            t = StaticUInt{1,T}(b.value[lastword]>>>extra)
            """

    for i in indices
        text *= "shiftedxor!(b, t, $i)\n"
    end

    text *= """b.value[lastword] &= (T(1)<<extra)-1
            return BFieldElt{$D,$R,T,lastword}(changelength(b, lastword))
            end"""
    return eval(Meta.parse(text))
end

@fastreduce(113, UInt16(512+1))
@fastreduce(131, UInt16(256+8+4+1))
@fastreduce(163, UInt16(128+64+8+1))
@fastreduce(193, (UInt16(1)<<15) + UInt16(1))
@fastreduce(233, (UInt128(1)<<74) + UInt128(1))
@fastreduce(239, (UInt64(1)<<36) + UInt64(1))
@fastreduce(283, (UInt16(1)<<12) + UInt16(128+32+1))
@fastreduce(409, (UInt128(1)<<87) + UInt128(1))
@fastreduce(571, (UInt16(1)<<10) + UInt16(32+4+1))
