import Base: convert, +, ==

function convert(::Type{ECPointAffine}, p::ECPointProjective)
    if iszero(p) return ECPointAffine(ec) end
    z_inv = inv(p.z)
    return ECPointAffine(p.x*z_inv, p.y*z_inv, p.ec)
end

function convert(::Type{ECPointAffine}, p::ECPointJacobian)
    if iszero(p) return ECPointAffine(ec) end
    z_inv = inv(p.z)
    return ECPointAffine(p.x*z_inv^2, p.y*z_inv^3, p.ec)
end

function convert(::Type{ECPointAffine}, p::ECPointLD)
    if iszero(p) return ECPointAffine(ec) end
    z_inv = inv(p.z)
    return ECPointAffine(p.x*z_inv, p.y*z_inv^2, p.ec)
end

function convert(::Type{ECPointProjective}, p::ECPointAffine)
    if iszero(p) return ECPointProjective(ec) end
    return ECPointProjective(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointProjective}, p::ECPointJacobian)
    if iszero(p) return ECPointProjective(ec) end
    z_inv = inv(p.z)
    return ECPointProjective(p.x*z_inv, p.y*(z_inv^2), p.z, p.ec)
end

function convert(::Type{ECPointProjective}, p::ECPointLD)
    if iszero(p) return ECPointProjective(ec) end
    return ECPointProjective(p.x, p.y / p.z, p.z, p.ec)
end

function convert(::Type{ECPointJacobian}, p::ECPointAffine)
    if iszero(p) return ECPointJacobian(ec) end
    return ECPointJacobian(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointJacobian}, p::ECPointProjective)
    if iszero(p) return ECPointJacobian(ec) end
    z2 = p.z^2
    return ECPointJacobian(p.x*z2, p.y*z2*p.z, p.z, p.ec)
end

function convert(::Type{ECPointJacobian}, p::ECPointLD)
    if iszero(p) return ECPointJacobian(ec) end
    return ECPointJacobian(p.x*p.z, p.y*p.z, p.z, p.ec)
end

function convert(::Type{ECPointLD}, p::ECPointAffine)
    if iszero(p) return ECPointLD(ec) end
    return ECPointLD(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointLD}, p::ECPointProjective)
    if iszero(p) return ECPointLD(ec) end
    return ECPointLD(p.x, p.y*p.z, p.z, p.ec)
end

function convert(::Type{ECPointLD}, p::ECPointJacobian)
    if iszero(p) return ECPointLD(ec) end
    z_inv = inv(p.z)
    return ECPointLD(p.x*z_inv, p.y*z_inv, p.z, p.ec)
end

function ==(p1::ECPointAffine, p2::ECPointProjective)
    return (iszero(p1) && iszero(p2)) || (p1.x*p2.z == p2.x && p1.y*p2.z == p2.y)
end

function ==(p1::ECPointProjective, p2::ECPointAffine)
    return p2==p1
end

function ==(p1::ECPointAffine, p2::ECPointLD)
    return (iszero(p1) && iszero(p2)) || (p1.x*p2.z == p2.x && p1.y*p2.z^2 == p2.y)
end

function ==(p1::ECPointLD, p2::ECPointAffine)
    return p2==p1
end

function ==(p1::ECPointAffine, p2::ECPointJacobian)
    z2_2 = p2.z^2
    return (iszero(p1) && iszero(p2)) || (p1.x*z2_2 == p2.x && p1.y*p2.z*z2_2 == p2.y)
end

function ==(p1::ECPointJacobian, p2::ECPointAffine)
    return p2==p1
end

function ==(p1::ECPointProjective, p2::ECPointLD)
    return (iszero(p1) && iszero(p2)) || (p1.x*p2.z == p2.x*p1.z && p1.y*p2.z^2 == p2.y*p1.z)
end

function ==(p1::ECPointLD, p2::ECPointProjective)
    return p2==p1
end

function ==(p1::ECPointProjective, p2::ECPointJacobian)
    z2_2 = p2.z^2
    return (iszero(p1) && iszero(p2)) || (p1.x*z2_2 == p2.x*p1.z && p1.y*p2.z*z2_2 == p2.y*p1.z)
end

function ==(p1::ECPointJacobian, p2::ECPointProjective)
    return p2==p1
end

function ==(p1::ECPointLD, p2::ECPointJacobian)
    z2_2 = p2.z^2
    return (iszero(p1) && iszero(p2)) || (p1.x*z2_2 == p2.x*p1.z && p1.y*p2.z*z2_2 == p2.y*p1.z^2)
end

function ==(p1::ECPointJacobian, p2::ECPointLD)
    return p2==p1
end

#add p1 and p2, return a jacobian
function +(::Type{ECPointJacobian}, p1::ECPointLD, p2::ECPointJacobian)
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return convert(ECPointJacobian, p1) end
    if p1==-p2 return ECPointJacobian(p2.ec) end
    if p1==p2 return double(p2) end

    #Adds: 7
    #Mults: 19
    #Sqrs: 5
    #Invs: 0

    z2_2 = p2.z^2
    z1_2 = p1.z^2
    A = p1.y*z2_2*p2.z + p2.y*z1_2
    C = p1.x*z2_2 + p2.x*p1.z
    B = p1.z*p2.z*C

    z_new = B*p1.z
    if iszero(z_new) return ECPointJacobian(p2.ec) end

    AB = A+B

    x_new = A*AB + p1.ec.a*B^2 + C^3*p1.z
    x_new *= z1_2

    y_new = A*p1.x + C*p1.y*p2.z
    y_new *= z_new^2
    y_new += x_new*p1.z*AB

    return ECPointJacobian(x_new, y_new, z_new, p2.ec)
end

function +(::Type{ECPointJacobian}, p1::ECPointJacobian, p2::ECPointLD)
    return +(ECPointJacobian, p2, p1)
end

#add p1 and p2, return an LD
function +(::Type{ECPointLD}, p1::ECPointLD, p2::ECPointJacobian)
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return convert(ECPointLD, p2) end
    if iszero(p2) return p1 end
    if p1==-p2 return ECPointLD(p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 7
    #Mults: 20
    #Sqrs: 4
    #Invs: 0

    z2_2 = p2.z^2
    A = p1.y*z2_2*p2.z + p2.y*p1.z^2
    C = p1.x*z2_2 + p2.x*p1.z
    B = p1.z*p2.z*C

    B2 = B^2

    z_new = B2*p1.z
    if iszero(z_new) return ECPointLD(p1.ec) end

    AB = A+B

    x_new = A*AB + p1.ec.a*B2 + C^3*p1.z
    x_new *= p1.z

    y_new = A*p1.x + C*p1.y*p2.z
    y_new *= z_new
    y_new += x_new*p1.z*AB
    y_new *= B

    if iszero(y_new) return ECPointLD(p1.ec) end

    return ECPointLD(x_new, y_new, z_new, p1.ec)
end

function +(::Type{ECPointLD}, p1::ECPointJacobian, p2::ECPointLD)
    return +(ECPointLD, p2, p1)
end
