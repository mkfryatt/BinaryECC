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
right_to_left_mult(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

threads_mult(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

noreduce_mult(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

right_to_left_comb_mult(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

left_to_right_comb_mult(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}) where {D,R}

window_comb_mult(a::BFieldPoint{D,R}, b::BFieldPoint{D,R}, window::Int) where {D,R}

square(a::BFieldPoint{D,R}) where {D,R}

window_square(a::BFieldPoint{D,R}, window::Int) where {D,R}
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
