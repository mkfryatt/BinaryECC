# Elliptic Curve Arithmetic

## Types
```@docs
EC{B}

AbstractECPoint{B}

ECPointAffine{B}

ECPointJacobian{B}

ECPointLD{B}

ECMismatchException
```

## Curve Functions
```@docs
==(ec1::EC{B}, ec2::EC{B}) where B

repr(ec::EC)
```

## Point Arithmetic
```@docs
==(p1::ECPointAffine{B}, p2::ECPointAffine{B}) where B

+(p1::ECPointAffine{B}, p2::ECPointAffine{B}) where B

double(p::ECPointAffine{B}) where B

-(p::ECPointAffine{B}) where B

-(p1::ECPointAffine{B}, p2::ECPointAffine{B}) where B

*(p::AbstractECPoint{B}, n::Integer) where B
```

### Scalar Point Addition
```@docs
double_standard(p::ECPointAffine{B}) where B

double_memo(p::ECPointAffine{B}) where B

double_threaded(p::ECPointAffine{B}) where B
```

### Scalar Point Multiplication
```@docs
mult_standard_rtl(P::AbstractECPoint{B}, k::Integer) where B

mult_standard_ltr(P::AbstractECPoint{B}, k::Integer) where B

mult_window(P::AbstractECPoint{B}, k::Integer, w::Int) where B

mult_bnaf(P::AbstractECPoint{B}, k::Integer) where B

mult_memo(P::AbstractECPoint{B}, k::Integer) where B

mult_bnaf_threaded(P::AbstractECPoint{B}, k::Integer) where B

mult_bnaf_window(P::AbstractECPoint{B}, k::Integer, w::Int) where B

mult_wnaf(P::AbstractECPoint{B}, k::Integer, w::Int) where B

mult_mont_general(p::AbstractECPoint, n::Integer)

mult_mont_affine(p::ECPointAffine{B}, n::Integer) where B
```

### Additive Identity
```@docs
iszero(p::ECPointAffine)
zero(::Type{ECPointAffine{B}}, ec::EC{B}) where B
zero(::Type{ECPointAffine}, ec::EC{B}) where B

iszero(p::ECPointJacobian)
zero(::Type{ECPointJacobian{B}}, ec::EC{B}) where B
zero(::Type{ECPointJacobian}, ec::EC{B}) where B

iszero(p::ECPointLD)
zero(::Type{ECPointLD{B}}, ec::EC{B}) where B
zero(::Type{ECPointLD}, ec::EC{B}) where B
```

## Mixed Representation Functions
```@docs
+(::Type{ECPointJacobian{B}}, p1::ECPointLD{B}, p2::ECPointJacobian{B}) where B
+(::Type{ECPointJacobian{B}}, p1::ECPointJacobian{B}, p2::ECPointLD{B}) where B

+(::Type{ECPointLD{B}}, p1::ECPointLD{B}, p2::ECPointJacobian{B}) where B
+(::Type{ECPointLD{B}}, p1::ECPointJacobian{B}, p2::ECPointLD{B}) where B

convert(::Type{ECPointAffine}, p::ECPointJacobian{B}) where B
convert(::Type{ECPointAffine}, p::ECPointLD{B}) where B
convert(::Type{ECPointJacobian}, p::ECPointAffine{B}) where B
convert(::Type{ECPointJacobian}, p::ECPointLD{B}) where B
convert(::Type{ECPointLD}, p::ECPointAffine{B}) where B
convert(::Type{ECPointLD}, p::ECPointJacobian{B}) where B
```

## Miscellaneous Point Functions
```@docs
ECPointAffine(s::String, ec::EC{BFieldElt{D,R,T,L}}) where {D,R,T,L}

isvalid(p::ECPointAffine)

repr(p::ECPointAffine)

repr(p::ECPointJacobian)

repr(p::ECPointLD)

naf
```
