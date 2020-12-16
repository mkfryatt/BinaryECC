import Base: +, -, *, ==, repr, isvalid, iszero
using Base: ceil

struct ECPointAffine <: AbstractECPoint
    x::FieldPoint
    y::FieldPoint
    isId::Bool
    ec::EC

    ECPointAffine(x::FieldPoint, y::FieldPoint, isId::Bool, ec::EC) =
        if x.field!=y.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, isId, ec) end

    ECPointAffine(x::FieldPoint, y::FieldPoint, ec::EC) =
        if x.field!=y.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, false, ec) end

    ECPointAffine(ec::EC) =
        new(FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), true, ec)
end

#sec1v2 2.3.4
function ECPointAffine(s::String, ec::EC)
    s = replace(s, " " => "")

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
    return p1.x==p2.x && p1.y==p2.y && iszero(p1)==iszero(p2) && p1.ec==p2.ec
end

function +(p1::ECPointAffine, p2::ECPointAffine)
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
    if p1==-p2 return ECPointAffine(p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 8
    #Mults: 2
    #Sqrs: 1
    #Invs: 1
    x1x2 = p1.x + p2.x
    lambda = (p1.y+p2.y) / x1x2
    x3 = lambda^2 + lambda + x1x2 + p1.ec.a
    y3 = lambda*(p1.x+x_new) + x3 + p1.y
    return ECPointAffine(x3, y3, p1.ec)
end

function -(p::ECPointAffine)
    if iszero(p) return p end
    return ECPointAffine(p.x, p.x+p.y, p.ec)
end

function double(p::ECPointAffine)
    if iszero(p) return p end
    if p==-p return ECPointAffine(p.ec) end

    #Adds: 5
    #Mults: 2
    #Sqrs: 2
    #Invs: 1
    lambda = p.x + (p.y / p.x)
    x_new = lambda^2 + lambda + p.ec.a
    y_new = p.x^2 + lambda*x_new + x_new
    return ECPointAffine(x_new, y_new, p.ec)
end

function *(p::ECPointAffine, n::Integer)
    if iszero(p) return p end
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

function isvalid(p::ECPointAffine)
    return iszero(p) || (p.y^2 + p.x*p.y == p.x^3 + p.ec.a*p.x^2 + p.ec.b)
end

function iszero(p::ECPointAffine)
    return p.isId
end
