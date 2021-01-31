# Elliptic Curve Arithmetic

## Types
```@docs
EC{D,R}

AbstractECPoint

ECPointAffine{D,R}

ECPointJacobian{D,R}

ECPointLD{D,R}

ECMismatchException
```

## Curve Functions
```@docs
==(ec1::EC{D,R}, ec2::EC{D,R}) where {D,R}

repr(ec::EC)
```

## Point Arithmetic
```@docs
==(p1::ECPointAffine{D,R}, p2::ECPointAffine{D,R}) where {D,R}

+(p1::ECPointAffine{D,R}, p2::ECPointAffine{D,R}) where {D,R}

-(p::ECPointAffine) where {D,R}

-(p1::AbstractECPoint, p2::AbstractECPoint)

*(p::ECPointAffine, n::Integer) where {D,R}
```

### Scalar Point Multiplication
```@docs
mont_pow_ladder

naf_mult(p::ECPointAffine{D,R}, n::Integer) where {D,R}
```

### Additive Identity
```@docs
iszero(p::ECPointAffine)
zero(::Type{ECPointAffine{D,R}}, ec::EC{D,R}) where {D,R}
zero(::Type{ECPointAffine}, ec::EC{D,R}) where {D,R}

iszero(p::ECPointJacobian)
zero(::Type{ECPointJacobian{D,R}}, ec::EC{D,R}) where {D,R}
zero(::Type{ECPointJacobian}, ec::EC{D,R}) where {D,R}

iszero(p::ECPointLD)
zero(::Type{ECPointLD{D,R}}, ec::EC{D,R}) where {D,R}
zero(::Type{ECPointLD}, ec::EC{D,R}) where {D,R}
```

## Mixed Representation Functions
```@docs
+(::Type{ECPointJacobian{D,R}}, p1::ECPointLD{D,R}, p2::ECPointJacobian{D,R}) where {D,R}
+(::Type{ECPointJacobian{D,R}}, p1::ECPointJacobian{D,R}, p2::ECPointLD{D,R}) where {D,R}

+(::Type{ECPointLD{D,R}}, p1::ECPointLD{D,R}, p2::ECPointJacobian{D,R}) where {D,R}
+(::Type{ECPointLD{D,R}}, p1::ECPointJacobian{D,R}, p2::ECPointLD{D,R}) where {D,R}

convert(::Type{ECPointAffine}, p::ECPointJacobian{D,R}) where {D,R}
convert(::Type{ECPointAffine}, p::ECPointLD{D,R}) where {D,R}
convert(::Type{ECPointJacobian}, p::ECPointAffine{D,R}) where {D,R}
convert(::Type{ECPointJacobian}, p::ECPointLD{D,R}) where {D,R}
convert(::Type{ECPointLD}, p::ECPointAffine{D,R}) where {D,R}
convert(::Type{ECPointLD}, p::ECPointJacobian{D,R}) where {D,R}
```

## Miscellaneous Point Functions
```@docs
ECPointAffine(s::String, ec::EC{D,R}) where {D,R}

ECPointAffine(ec::EC{D,R}) where {D,R}

isvalid(p::ECPointAffine)

repr(p::ECPointAffine)

repr(p::ECPointJacobian)

repr(p::ECPointLD)
```
