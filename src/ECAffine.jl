import Base: +, -, *, ==, !=, repr

#TODO: add supersingular version? would need a "c" param
struct ECAffine <: ECAbstract
    a::FieldPoint
    b::FieldPoint
    ECAffine(a::FieldPoint, b::FieldPoint) =
        if a.field!=b.field throw(FieldMismatchException())
        else new(a, b) end
end

struct ECPointAffine <: ECPointAbstract
    x::FieldPoint
    y::FieldPoint
    isId::Bool
    ec::ECAffine

    ECPointAffine(x::FieldPoint, y::FieldPoint, isId::Bool, ec::ECAffine) =
        if x.field!=y.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, isId, ec) end

    ECPointAffine(x::FieldPoint, y::FieldPoint, ec::ECAffine) =
        if x.field!=y.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, false, ec) end

    ECPointAffine(ec::ECAffine) =
        new(FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), true, ec)
end

function repr(ec::ECAffine)
    return "Equation: y^2 + xy = x^3 + "*repr(ec.a.x)*"x^2 + "*repr(ec.b.x)
end

function repr(p::ECPointAffine)
    return "("*repr(p.x)*", "*repr(p.y)*")"
end

function ==(ec1::ECAffine, ec2::ECAffine)
    return ec1.a==ec2.a && ec1.b==ec2.b && ec1.a.field==ec2.a.field
end

function ==(p1::ECPointAffine, p2::ECPointAffine)
    return p1.x==p2.x && p1.y==p2.y && p1.isId==p2.isId && p1.ec==p2.ec
end

function !=(ec1::ECAffine, ec2::ECAffine)
    return ec1.a!=ec2.a || ec1.b!=ec2.b || ec1.a.field!=ec2.a.field
end

function !=(p1::ECPointAffine, p2::ECPointAffine)
    return p1.x!=p2.x || p1.y!=p2.y || p1.isId!=p2.isId || p1.ec!=p2.ec
end

function +(p1::ECPointAffine, p2::ECPointAffine)
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if p1.isId return p2 end
    if p2.isId return p1 end
    if p1==p2 return double(p1) end

    lambda = (p1.y+p2.y) / (p1.x +p2.x)
    x_new = lambda^2 + lambda + p1.x + p2.x + p1.ec.a
    y_new = lambda*(p1.x+x_new) + x_new + p1.y
    return ECPointAffine(x_new, y_new, p1.ec)
end

function -(p::ECPointAffine)
    if p.isId return p end
    return ECPointAffine(p.x, p.x+p.y, p.ec)
end

function -(p1::ECPointAffine, p2::ECPointAffine)
    return p1 + (-p2)
end

function double(p::ECPointAffine)
    if p.isId return p end

    lambda = p1.x + (p1.y / p1.x)
    x_new = lambda^2 + lambda + p1.ec.a
    y_new = p1.x^2 + lambda*x_new + x_new
    return ECPointAffine(x_new, y_new, p.ec)
end

function *(p::ECPointAffine, n::Number)
    if p.isId return p end
    if n==1 return p end

    result = ECPointAffine(p.ec)
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

function is_valid(p::ECPointAffine)
    return p.isId || (p.y^2 + p.x*p.y == p.x^3 + p.ec.a*p.x^2 + p.ec.b)
end

#sec1v2 2.3.4
function from_octet_string(s::String, ec::ECAffine)
    #point is id
    if s=="00" return ECPointAffine(ec) end

    #input string is not of the uncompressed format specified by sec1v2
    if length(s) != 4*floor(Int16, ec.a.field.order / 8)+2
        throw(ArgumentError("Octet string is of the wrong length for this curve."))
    end
    if s[1:2]!="04"
        throw(ArgumentError("Octet string must start with '04'.")) 
    end

    #TODO handle compressed format

    x = from_octet_string(s[3:2+2*floor(Int16, ec.a.field.order / 8)], ec.a.field)
    y = from_octet_string(s[3+2*floor(Int16, ec.a.field.order / 8):2+4*floor(Int16, ec.a.field.order / 8)], ec.a.field)

    return ECPointAffine(x, y, ec)
end

#TODO random point on curve
#TODO find y coord given x
#TODO find x coord given y

#f = Field(4, 19)
#a = 8, b = 9
#valid point = (3, 12)
