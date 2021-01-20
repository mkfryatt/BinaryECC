# Binary Field Arithmetic

## Types
```@docs
FieldPoint{D,R}

FieldPoint{D,R}(s::String) where {D,R}
```

## Arithmetic
```@docs
==(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

+(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

-(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

-(a::FieldPoint{D,R}) where {D,R}

*(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

inv(a::FieldPoint{D,R}) where {D,R}

/(a::FieldPoint{D,R}, b::FieldPoint{D,R}) where {D,R}

^(a::FieldPoint{D,R}, b::Integer) where {D,R}

iszero(a::FieldPoint)

zero(::Type{FieldPoint{D,R}}) where {D,R}

isone(a::FieldPoint)

one(::Type{FieldPoint{D,R}}) where {D,R}
```

## Miscellaneous Functions
```@docs
random(::Type{FieldPoint{D,R}}) where {D,R}

convert(::Type{BigInt}, a::FieldPoint)
```
