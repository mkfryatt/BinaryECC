import Base: +, -, *, ==, !=, repr, ceil, convert

struct ECAffine <: ECAbstract
    a::FieldPoint
    b::FieldPoint
    ECAffine(a::FieldPoint, b::FieldPoint) =
        if a.field!=b.field throw(FieldMismatchException())
        else new(a, b) end
end

function convert(::Type{ECAffine}, ec::ECAbstract)
    return ECAffine(ec.a, ec.b)
end

function repr(ec::ECAffine)
    return "E: y² + xy = x³ + "*repr(ec.a.x)*"x² + "*repr(ec.b.x)
end

struct ECPointAffine <: ECPointAbstract
    x::FieldPoint
    y::FieldPoint
    isId::Bool
    ec::ECAffine

    ECPointAffine(x::FieldPoint, y::FieldPoint, isId::Bool, ec::ECAbstract) =
        if x.field!=y.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, isId, convert(ECAffine, ec)) end

    ECPointAffine(x::FieldPoint, y::FieldPoint, ec::ECAbstract) =
        if x.field!=y.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, false, convert(ECAffine, ec)) end

    ECPointAffine(ec::ECAbstract) =
        new(FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), true, convert(ECAffine, ec))
end

#sec1v2 2.3.4
function ECPointAffine(s::String, ec::ECAbstract)
    s = replace(s, " " => "")
    ec = convert(ECAffine, ec)

    #point is id
    if s=="00" return ECPointAffine(ec) end

    #input string is not of the uncompressed format specified by sec1v2
    if length(s) != 4*ceil(Int16, ec.a.field.order / 8)+2
        throw(ArgumentError("Octet string is of the wrong length for this curve."))
    end
    if s[1:2]!="04"
        throw(ArgumentError("Octet string must start with '04'."))
    end

    #TODO handle compressed format

    x = FieldPoint(s[3:4+2*floor(Int16, ec.a.field.order / 8)], ec.a.field)
    y = FieldPoint(s[5+2*floor(Int16, ec.a.field.order / 8):6+4*floor(Int16, ec.a.field.order / 8)], ec.a.field)

    return ECPointAffine(x, y, ec)
end

function repr(p::ECPointAffine)
    return "("*repr(p.x)*", "*repr(p.y)*")"
end

function ==(p1::ECPointAffine, p2::ECPointAffine)
    return p1.x==p2.x && p1.y==p2.y && p1.isId==p2.isId && p1.ec==p2.ec
end

function !=(p1::ECPointAffine, p2::ECPointAffine)
    return p1.x!=p2.x || p1.y!=p2.y || p1.isId!=p2.isId || p1.ec!=p2.ec
end

function +(p1::ECPointAffine, p2::ECPointAffine)
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if p1.isId return p2 end
    if p2.isId return p1 end
    if p1==-p2 return ECPointAffine(p1.ec) end
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

function double(p::ECPointAffine)
    if p.isId return p end
    if p==-p return ECPointAffine(p.ec) end

    lambda = p.x + (p.y / p.x)
    x_new = lambda^2 + lambda + p.ec.a
    y_new = p.x^2 + lambda*x_new + x_new
    return ECPointAffine(x_new, y_new, p.ec)
end

function *(p::ECPointAffine, n::Integer)
    if p.isId return p end
    if n==0 return ECPointAffine(p.ec) end
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
