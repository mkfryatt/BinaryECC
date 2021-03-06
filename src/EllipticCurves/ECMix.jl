"""
    convert(::Type{ECPointAffine}, p::ECPointJacobian{B}) where B
Converts a point from Jacobian coordinates to affine coordinates.
"""
function convert(::Type{ECPointAffine}, p::ECPointJacobian{B}) where B
    if iszero(p) return ECPointAffine{B}(ec) end
    z_inv = inv(p.z)
    z_inv2 = square(z_inv)
    return ECPointAffine{B}(p.x*z_inv2, p.y*z_inv2*z_inv, p.ec)
end

"""
    convert(::Type{ECPointAffine}, p::ECPointLD{B}) where B
Converts a point from Lopez-Dahab coordinates to affine coordinates.
"""
function convert(::Type{ECPointAffine}, p::ECPointLD{B}) where B
    if iszero(p) return ECPointAffine{B}(ec) end
    z_inv = inv(p.z)
    return ECPointAffine{B}(p.x*z_inv, p.y*square(z_inv), p.ec)
end

"""
    convert(::Type{ECPointJacobian}, p::ECPointAffine{B}) where B
Converts a point from affine coordinates to Jacobian coordinates.
"""
function convert(::Type{ECPointJacobian}, p::ECPointAffine{B}) where B
    if iszero(p) return ECPointJacobian{B}(ec) end
    return ECPointJacobian{B}(p.x, p.y, B(1), p.ec)
end

"""
    convert(::Type{ECPointJacobian}, p::ECPointLD{B}) where B
Converts a point from Lopez-Dahab coordinates to Jacobian coordinates.
"""
function convert(::Type{ECPointJacobian}, p::ECPointLD{B}) where B
    if iszero(p) return ECPointJacobian{B}(ec) end
    return ECPointJacobian{B}(p.x*p.z, p.y*p.z, p.z, p.ec)
end

"""
    convert(::Type{ECPointLD}, p::ECPointAffine{B}) where B
Converts a point from affine coordinates to Lopez-Dahab coordinates.
"""
function convert(::Type{ECPointLD}, p::ECPointAffine{B}) where B
    if iszero(p) return ECPointLD{B}(ec) end
    return ECPointLD{B}(p.x, p.y, B(1), p.ec)
end

"""
    convert(::Type{ECPointLD}, p::ECPointJacobian{B}) where B
Converts a point from Jacobian coordinates to Lopez-Dahab coordinates.
"""
function convert(::Type{ECPointLD}, p::ECPointJacobian{B}) where B
    if iszero(p) return ECPointLD{B}(ec) end
    z_inv = inv(p.z)
    return ECPointLD{B}(p.x*z_inv, p.y*z_inv, p.z, p.ec)
end

function ==(p1::ECPointAffine{B}, p2::ECPointLD{B}) where B
    return (iszero(p1) && iszero(p2)) || (p1.x*p2.z == p2.x && p1.y*square(p2.z) == p2.y)
end

function ==(p1::ECPointLD{B}, p2::ECPointAffine{B}) where B
    return p2==p1
end

function ==(p1::ECPointAffine{B}, p2::ECPointJacobian{B}) where B
    z2_2 = square(p2.z)
    return (iszero(p1) && iszero(p2)) || (p1.x*z2_2 == p2.x && p1.y*p2.z*z2_2 == p2.y)
end

function ==(p1::ECPointJacobian{B}, p2::ECPointAffine{B}) where B
    return p2==p1
end

function ==(p1::ECPointLD{B}, p2::ECPointJacobian{B}) where B
    z2_2 = square(p2.z)
    return (iszero(p1) && iszero(p2)) || (p1.x*z2_2 == p2.x*p1.z && p1.y*p2.z*z2_2 == p2.y*square(p1.z))
end

function ==(p1::ECPointJacobian{B}, p2::ECPointLD{B}) where B
    return p2==p1
end

function +(p1::AbstractECPoint{B}, p2::AbstractECPoint{B})::ECPointAffine{B} where B
    return convert(ECPointAffine, p1) + convert(ECPointAffine, p2)
end    

"""
    +(::Type{ECPointJacobian{B}}, p1::ECPointLD{B}, p2::ECPointJacobian{B}) where B
Returns ``p_1 +p_2`` in Jacobian coordinates.
"""
function +(::Type{ECPointJacobian{BF}}, p1::ECPointLD{BF}, p2::ECPointJacobian{BF}) where BF
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return convert(ECPointJacobian{BF}, p1) end
    if p1==-p2 return ECPointJacobian{BF}(p2.ec) end
    if p1==p2 return double(p2) end

    #Adds: 7
    #Mults: 19
    #Sqrs: 5
    #Invs: 0

    z2_2 = square(p2.z)
    z1_2 = square(p1.z)
    A = p1.y*z2_2*p2.z + p2.y*z1_2
    C = p1.x*z2_2 + p2.x*p1.z
    B = p1.z*p2.z*C

    z_new = B*p1.z

    AB = A+B

    x_new = A*AB + p1.ec.a*square(B) + C^3*p1.z
    x_new *= z1_2

    y_new = A*p1.x + C*p1.y*p2.z
    y_new *= square(z_new)
    y_new += x_new*p1.z*AB

    return ECPointJacobian{BF}(x_new, y_new, z_new, p2.ec)
end

"""
    +(::Type{ECPointJacobian{B}}, p1::ECPointJacobian{B}, p2::ECPointLD{B}) where B
Returns ``p_1 +p_2`` in Jacobian coordinates.
"""
function +(::Type{ECPointJacobian{B}}, p1::ECPointJacobian{B}, p2::ECPointLD{B}) where B
    return +(ECPointJacobian{B}, p2, p1)
end

"""
    +(::Type{ECPointLD{B}}, p1::ECPointLD{B}, p2::ECPointJacobian{B}) where B
Returns ``p_1 +p_2`` in Lopez-Dahab coordinates.
"""
function +(::Type{ECPointLD{BF}}, p1::ECPointLD{BF}, p2::ECPointJacobian{BF}) where BF
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return convert(ECPointLD{BF}, p2) end
    if iszero(p2) return p1 end
    if p1==-p2 return ECPointLD{BF}(p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 7
    #Mults: 20
    #Sqrs: 4
    #Invs: 0

    z2_2 = square(p2.z)
    A = p1.y*z2_2*p2.z + p2.y*square(p1.z)
    C = p1.x*z2_2 + p2.x*p1.z
    B = p1.z*p2.z*C

    B2 = square(B)

    z_new = B2*p1.z

    AB = A+B

    x_new = A*AB + p1.ec.a*B2 + C^3*p1.z
    x_new *= p1.z

    y_new = A*p1.x + C*p1.y*p2.z
    y_new *= z_new
    y_new += x_new*p1.z*AB
    y_new *= B

    return ECPointLD{BF}(x_new, y_new, z_new, p1.ec)
end

"""
    +(::Type{ECPointLD{B}}, p1::ECPointJacobian{B}, p2::ECPointLD{B}) where B
Returns ``p_1 +p_2`` in Lopez-Dahab coordinates.
"""
function +(::Type{ECPointLD{B}}, p1::ECPointJacobian{B}, p2::ECPointLD{B}) where B
    return +(ECPointLD{B}, p2, p1)
end
