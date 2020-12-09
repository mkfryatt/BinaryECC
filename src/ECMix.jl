import Base: convert

function convert(::Type{ECPointAffine}, p::ECPointProjective)
    if p.isId return ECPointAffine(ec) end
    z_inv = inv(p.z)
    return ECPointAffine(p.x*z_inv, p.y*z_inv, p.ec)
end

function convert(::Type{ECPointAffine}, p::ECPointJacobian)
    if p.isId return ECPointAffine(ec) end
    z_inv = inv(p.z)
    return ECPointAffine(p.x*z_inv^2, p.y*z_inv^3, p.ec)
end

function convert(::Type{ECPointAffine}, p::ECPointProjective)
    if p.isId return ECPointAffine(ec) end
    z_inv = inv(p.z)
    return ECPointAffine(p.x*z_inv, p.y*z_inv^2, p.ec)
end

function convert(::Type{ECPointProjective}, p::ECPointAffine)
    if p.isId return ECPointProjective(ec) end
    return ECPointProjective(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointProjective}, p::ECPointJacobian)
    #TODO
    if p.isId return ECPointProjective(ec) end
    return ECPointProjective(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointProjective}, p::ECPointLD)
    #TODO
    if p.isId return ECPointProjective(ec) end
    return ECPointProjective(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointJacobian}, p::ECPointAffine)
    if p.isId return ECPointJacobian(ec) end
    return ECPointJacobian(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointJacobian}, p::ECPointProjective)
    #TODO
    if p.isId return ECPointJacobian(ec) end
    return ECPointJacobian(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointJacobian}, p::ECPointLD)
    #TODO
    if p.isId return ECPointJacobian(ec) end
    return ECPointJacobian(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointLD}, p::ECPointAffine)
    if p.isId return ECPointLD(ec) end
    return ECPointLD(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointLD}, p::ECPointProjective)
    #TODO
    if p.isId return ECPointLD(ec) end
    return ECPointLD(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end

function convert(::Type{ECPointLD}, p::ECPointJacobian)
    #TODO
    if p.isId return ECPointLD(ec) end
    return ECPointLD(p.x, p.y, FieldPoint(1, p.x.field), p.ec)
end
