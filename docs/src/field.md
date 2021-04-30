# Binary Field Arithmetic

## Types
```@docs
BFieldElt{D,R,T,L}

BFieldElt{D,R,T,L}(s::String) where {D,R,T,L}
```

## Arithmetic

### General Arithmetic
```@docs
==(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}

+(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}

-(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}

*(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}

inv(a::BFieldElt{D,R,T,L}) where {D,R,T,L}

/(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}

^(a::BFieldElt{D,R,T,L}, b::Integer) where {D,R,T,L}

square(a::BFieldElt{D,R,T,L}) where {D,R,T,L}

sqrt(a::BFieldElt{D,R,T,L}) where {D,R,T,L}
```

### Multiplication
```@docs
mult_shiftandadd(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}

mult_shiftandadd_window(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}, w::Int) where {D,R,T,L}

mult_threaded(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}

mult_threaded_window(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}, w::Int) where {D,R,T,L}

mult_ownreduce(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}

mult_comb_rtl(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}

mult_comb_ltr(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}) where {D,R,T,L}

mult_comb_window(a::BFieldElt{D,R,T,L}, b::BFieldElt{D,R,T,L}, window::Int) where {D,R,T,L}

square_standard(a::BFieldElt{D,R,T,L}) where {D,R,T,L}

square_window(a::BFieldElt{D,R,T,L}, window::Int) where {D,R,T,L}
```

### Reduction
```@docs
reduce(a::BFieldElt{D,R,T,L}) where {D,R,T,L}

BinaryECC.@fastreduce
```

### Miscellaneous Arithmetic
```@docs
iszero(a::BFieldElt)

zero(::Type{BFieldElt{D,R,T,L}}) where {D,R,T,L}

isone(a::BFieldElt)

one(::Type{BFieldElt{D,R,T,L}}) where {D,R,T,L}
```

## Miscellaneous Functions
```@docs
random(::Type{BFieldElt{D,R,T,L}}) where {D,R,T,L}

convert(::Type{BigInt}, a::BFieldElt)
```
