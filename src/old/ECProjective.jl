import Base: +, -, *, ==, repr, isvalid, iszero

struct ECPointProjective <: AbstractECPoint
    x::FieldPoint
    y::FieldPoint
    z::FieldPoint
    ec::EC

    ECPointProjective(x::FieldPoint, y::FieldPoint, z::FieldPoint, ec::EC) =
        if x.field!=y.field || x.field!=z.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, z, ec) end

    ECPointProjective(x::FieldPoint, y::FieldPoint, z::FieldPoint, ec::EC) =
        if x.field!=y.field || x.field!=z.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, z, ec) end

    ECPointProjective(ec::EC) =
        new(FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), ec)
end

function repr(p::ECPointProjective)
    return "("*repr(p.x)*", "*repr(p.y)*", "*repr(p.z)*")"
end

function ==(p1::ECPointProjective, p2::ECPointProjective)
    return iszero(p1)==iszero(p2) && p1.ec==p2.ec && p1.x*p2.z==p2.x*p1.z && p1.y*p2.z==p2.y*p1.z
end

function +(p1::ECPointProjective, p2::ECPointProjective)
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
    if p1==-p2 return ECPointProjective(p1.ec) end
    if p1==p2 return double(p1) end
    #TODO
end

function -(p::ECPointProjective)
    if iszero(p) return p end
    return ECPointProjective(p.x, p.x+p.y, p.z, p.ec)
end

function double(p::ECPointProjective)
    if iszero(p) return p end
    if p==-p return ECPointProjective(p.ec) end
    #TODO
end

function *(p::ECPointProjective, n::Integer)
    if iszero(p) return p end
    if n==0 return ECPointProjective(p.ec) end
    if n==1 return p end

    result = ECPointProjective(p.ec)
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

function isvalid(p::ECPointProjective)
    return iszero(p) || (p.y^2*p.z + p.x*p.y*p.z == p.x^3 + p.ec.a*p.x^2*p.z + p.ec.b*p.z^3)
end

function iszero(p::ECPointProjective)
    return iszero(p.x) && iszero(p.z)
end
