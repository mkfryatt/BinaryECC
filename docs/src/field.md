# Binary Field Arithmetic

## Types
```@docs
FieldPoint{D,R}

FieldPoint{D,R}(s::String) where {D,R}
```

## Arithmetic

### General Arithmetic
```@docs
==(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

+(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

-(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

*(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

inv(a::FieldPoint{D,R}) where {D,R}

/(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

^(a::FieldPoint{D,R}, b::Integer) where {D,R}

sqrt(a::FieldPoint{D,R}) where {D,R}
```

### Multiplication
```@docs
right_to_left_mult(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

threads_mult(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

noreduce_mult(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

right_to_left_comb_mult(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

left_to_right_comb_mult(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

window_comb_mult(a::FieldPoint{D,R}, b::FieldPoint{D,R}, window::Int) where {D,R}

square(a::FieldPoint{D,R}) where {D,R}

window_square(a::FieldPoint{D,R}, window::Int) where {D,R}
```

### Miscellaneous Arithmetic
```@docs
iszero(a::FieldPoint)

zero(::Type{FieldPoint{D,R}}) where {D,R}

isone(a::FieldPoint)

one(::Type{FieldPoint{D,R}}) where {D,R}

reduce(a::FieldPoint{D,R}) where {D,R}
```

## Miscellaneous Functions
```@docs
random(::Type{FieldPoint{D,R}}) where {D,R}

convert(::Type{BigInt}, a::FieldPoint)
```
