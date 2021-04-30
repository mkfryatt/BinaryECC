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
        for i in length(b):-1:(1+lastword)
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
