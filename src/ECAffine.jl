import Base: +, -, *, ==, !=, repr

struct ECMismatchException <: Exception end

#TODO: add supersingular version?
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
    return "Equation: y^2 + xy = x^3 + "*repr(ec.a.x)*"x^2 + "*repr(ec.b.x)*"\n"*repr(ec.field)
end

function repr(p::ECAffinePoint)
    return "("*repr(p.x)*", "*repr(p.y)*")"
end

function ==(ec1::ECAffine, ec2::ECAffine)
    return ec1.a==ec2.a && ec1.b==ec2.b && ec1.field==ec2.field &&
        ec1.factor==ec2.factor && ec1.cofactor==ec2.cofactor
end

function ==(p1::ECAffinePoint, p2::ECAffinePoint)
    return p1.x==p2.x && p1.y==p2.y && p1.isId==p2.isId && p1.ec==p2.ec
end

function !=(ec1::ECAffine, ec2::ECAffine)
    return ec1.a!=ec2.a || ec1.b!=ec2.b || ec1.field!=ec2.field ||
        ec1.factor!=ec2.factor || ec1.cofactor!=ec2.cofactor
end

function !=(p1::ECAffinePoint, p2::ECAffinePoint)
    return p1.x!=p2.x || p1.y!=p2.y || p1.isId!=p2.isId || p1.ec!=p2.ec
end

function +(p1::ECAffinePoint, p2::ECAffinePoint)
    if p1.ec!=p2.ec throw(ECMismatchException) end
    if p1.isId return p2 end
    if p2.isId return p1 end
    if p1==p2 return double(p1) end

    lambda = (p1.y+p2.y) / (p1.x +p2.x)
    x_new = lambda^2 + lambda + p1.x + p2.x + p1.ec.a
    y_new = lambda*(p1.x+x_new) + x_new + p1.y
    return ECAffinePoint(x_new, y_new, p1.ec)
end

function -(p::ECAffinePoint)
    if p.isId return p end
    return ECAffinePoint(p.x, p.x+p.y, p.ec)
end

function -(p1::ECAffinePoint, p2::ECAffinePoint)
    return p1 + (-p2)
end

function double(p::ECAffinePoint)
    if p.isId return p end

    lambda = p1.x + (p1.y / p1.x)
    x_new = lambda^2 + lambda + p1.ec.a
    y_new = p1.x^2 + lambda*x_new + x_new
    return ECAffinePoint(x_new, y_new, p.ec)
end

function *(p::ECAffinePoint, n::Number)
    if p.isId return p end
    if n==1 return p end

    result = ECAffinePoint(true, p.ec)
    doubling = p
    while n>0
        if n&1==1
            result += doubling
        end
        doubling = double(doubling)
        n >>= 1
    end
    return result
end

function isValid(p::ECAffinePoint)
    return p.isId || (p.y^2 + p.x*p.y == p.x^3 + p.ec.a*p.x^2 + p.ec.b)
end

#f = Field(4, 19)
#a = 8, b = 9
#valid point = (3, 12)
