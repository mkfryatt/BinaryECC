import Base: +, -, *, ==, !=, repr

struct ECAffine
    a::FieldPoint
    b::FieldPoint
    field::Field
    factor::BigInt
    cofactor::BigInt
    ECAffine(a::FieldPoint, b::FieldPoint, field::Field, f::Number, c::Number) =
        new(a, b, field, convert(BigInt, f), convert(BigInt, c))
end

struct ECAffinePoint
    x::FieldPoint
    y::FieldPoint
    isId::Bool
    ec::ECAffine
    ECAffinePoint(x::FieldPoint, y::FieldPoint, isId::Bool, ec::ECAffine) = new(x, y, isId, ec)
    ECAffinePoint(x::FieldPoint, y::FieldPoint, ec::ECAffine) = new(x, y, false, ec)
    ECAffinePoint(isId::Bool, ec::ECAffine) =
        new(FieldPoint(0, ec.field), FieldPoint(0, ec.field), isId, ec)
end

function copy(ec::ECAffine)
    ECAffine(ec.a, ec.b, ec.field, ec.factor, ec.cofactor)
end

function copy(p::ECAffinePoint)
    return ECPoint(p.x, p.y, p.isId, p.ec)
end

function repr(ec::ECAffine)
    return "Equation: y^2 + xy = x^3 + "*repr(ec.a)*"x^2 + "*repr(ec.b)*"\n"*repr(ec.field)
end

function repr(p::ECAffinePoint)
    return "("*repr(p.x)*", "*repr(p.y)*")"
end

function ==(ec1::ECAffine, ec2::ECAffine)
    return ec1.a==ec2.a && ec1.b==ec2.b && ec1.field==ec2.field && ec1.factor==ec2.factor && ec1.cofactor==ec2.cofactor
end

function ==(p1::ECAffinePoint, p2::ECAffinePoint)
    return p1.x==p2.x && p1.y==p2.y && p1.isId==p2.isId && p1.ec==p2.ec
end

function !=(ec1::ECAffine, ec2::ECAffine)
    return ec1.a!=ec2.a || ec1.b!=ec2.b || ec1.field!=ec2.field || ec1.factor!=ec2.factor || ec1.cofactor!=ec2.cofactor
end

function !=(p1::ECAffinePoint, p2::ECAffinePoint)
    return p1.x!=p2.x || p1.y!=p2.y || p1.isId!=p2.isId || p1.ec!=p2.ec
end

function +(p1::ECAffinePoint, p2::ECAffinePoint)
    #TODO
end

function -(p::ECAffinePoint)
    #TODO
end

function -(p1::ECAffinePoint, p2::ECAffinePoint)
    return p1 + (-p2)
end

function double(p::ECAffinePoint)
    #TODO
end

function *(p::ECAffinePoint, n::Number)
    #TODO
end
