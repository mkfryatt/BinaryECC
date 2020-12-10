import Base: +, -, *, ==, !=, repr, convert

struct ECProjective <: ECAbstract
    a::FieldPoint
    b::FieldPoint
    ECAffine(a::FieldPoint, b::FieldPoint) =
        if a.field!=b.field throw(FieldMismatchException())
        else new(a, b) end
end

function convert(::Type{ECProjective}, ec::ECAbstract)
    return ECProjective(ec.a, ec.b)
end

function repr(ec::ECProjective)
    return "E: y²z + xyz = x³ + "*repr(ec.a.x)*"x²z + "*repr(ec.b.x)*"z³"
end

struct ECPointProjective <: ECPointAbstract
    x::FieldPoint
    y::FieldPoint
    z::FieldPoint
    isId::Bool
    ec::ECProjective

    ECPointProjective(x::FieldPoint, y::FieldPoint, z::FieldPoint, isId::Bool, ec::ECAbstract) =
        if x.field!=y.field || x.field!=z.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, z, isId, convert(ECProjective, ec)) end

    ECPointProjective(x::FieldPoint, y::FieldPoint, z::FieldPoint, ec::ECAbstract) =
        if x.field!=y.field || x.field!=z.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, z, false, convert(ECProjective, ec)) end

    ECPointProjective(ec::ECAbstract) =
        new(FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), true, convert(ECProjective, ec))
end

function repr(p::ECPointProjective)
    return "("*repr(p.x)*", "*repr(p.y)*", "*repr(p.z)*")"
end

function ==(p1::ECPointProjective, p2::ECPointProjective)
    return p1.isId==p2.isId && p1.ec==p2.ec && p1.x*p2.z==p2.x*p1.z && p1.y*p2.z==p2.y*p1.z
end

function !=(p1::ECPointProjective, p2::ECPointProjective)
    return p1.isId!=p2.isId || p1.ec!=p2.ec || p1.x*p2.z!=p2.x*p1.z || p1.y*p2.z!=p2.y*p1.z
end

function +(p1::ECPointProjective, p2::ECPointProjective)
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if p1.isId return p2 end
    if p2.isId return p1 end
    if p1==-p2 return ECPointProjective(p1.ec) end
    if p1==p2 return double(p1) end
    #TODO
end

function -(p::ECPointProjective)
    if p.isId return p end
    return ECPointProjective(p.x, p.x+p.y, p.z, p.ec)
end

function double(p::ECPointProjective)
    if p.isId return p end
    if p==-p return ECPointProjective(p.ec) end
    #TODO
end

function *(p::ECPointProjective, n::Integer)
    if p.isId return p end
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

function is_valid(p::ECPointProjective)
    return p.isId || (p.y^2*p.z + p.x*p.y*p.z == p.x^3 + p.ec.a*p.x^2*p.z + p.ec.b*p.z^3)
end
