# Elliptic Curve Arithmetic

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

montmul

iszero(p::ECPointAffine)
```

## Miscellaneous Point Functions
```@docs
ECPointAffine(s::String, ec::EC{D,R}) where {D,R}

ECPointAffine(ec::EC{D,R}) where {D,R}

isvalid(p::ECPointAffine)

repr(p::ECPointAffine)

repr(p::ECPointProjective)

repr(p::ECPointJacobian)

repr(p::ECPointLD)
```
