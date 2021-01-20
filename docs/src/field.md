```@docs
+(x1::FieldPoint{D,R}, x2::FieldPoint{D,R}) where {D,R}
```

# Binary Field Arithmetic

## Types
```@docs
EC{D,R}

AbstractECPoint

ECPointAffine{D,R}

ECPointProjective{D,R}

ECPointJacobian{D,R}

ECPointLD{D,R}

ECMismatchException
```

## Elliptic Curve Functions
```@docs
==(ec1::EC{D,R}, ec2::EC{D,R}) where {D,R}
```

## Elliptic Curve Point Functions
```@docs
ECPointAffine(s::String, ec::EC{D,R}) where {D,R}

ECPointAffine(ec::EC{D,R}) where {D,R}

==(p1::ECPointAffine{D,R}, p2::ECPointAffine{D,R}) where {D,R}

+(::ECPointAffine{D,R},::ECPointAffine{D,R}) where {D,R}

-(p::ECPointAffine) where {D,R}

-(p1::AbstractECPoint, p2::AbstractECPoint)

*(p::ECPointAffine, n::Integer) where {D,R}

montpow(p::AbstractECPoint, n::Integer)

isvalid(p::ECPointAffine)

iszero(p::ECPointAffine)
```
