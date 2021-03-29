# Binary Field Arithmetic

## Types
```@docs
BFieldPoint{D,R}

BFieldPoint{D,R}(s::String) where {D,R}
```

## Arithmetic

### General Arithmetic
```@docs
==(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

+(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

-(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

*(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

inv(a::BFieldPoint{D,R}) where {D,R}

/(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

^(a::BFieldPoint{D,R}, b::Integer) where {D,R}

sqrt(a::BFieldPoint{D,R}) where {D,R}
```

### Multiplication
```@docs
mult_shiftandadd(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}, w=0) where {D,R}

mult_threaded(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

mult_ownreduce(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

mult_comb_rtl(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

mult_comb_ltr(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

mult_comb_window(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}, window::Int) where {D,R}

square(a::BFieldPoint{D,R}) where {D,R}

square_window(a::BFieldPoint{D,R}, window::Int) where {D,R}
```

### Miscellaneous Arithmetic
```@docs
iszero(a::BFieldPoint)

zero(::Type{BFieldPoint{D,R}}) where {D,R}

isone(a::BFieldPoint)

one(::Type{BFieldPoint{D,R}}) where {D,R}

reduce(a::BFieldPoint{D,R}) where {D,R}
```

## Miscellaneous Functions
```@docs
random(::Type{BFieldPoint{D,R}}) where {D,R}

convert(::Type{BigInt}, a::BFieldPoint)
```
