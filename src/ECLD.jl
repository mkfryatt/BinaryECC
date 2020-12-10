import Base: +, -, *, ==, !=, repr, convert

struct ECLD <: ECAbstract
    a::FieldPoint
    b::FieldPoint
    ECLD(a::FieldPoint, b::FieldPoint) =
        if a.field!=b.field throw(FieldMismatchException())
        else new(a, b) end
end

function convert(::Type{ECLD}, ec::ECAbstract)
    return ECLD(ec.a, ec.b)
end

function repr(ec::ECLD)
    return "E: y² + xyz = x³z + "*repr(ec.a.x)*"x²z² + "*repr(ec.b.x)*"z⁴"
end

struct ECPointLD <: ECPointAbstract
    x::FieldPoint
    y::FieldPoint
    z::FieldPoint
    isId::Bool
    ec::ECLD

    ECPointLD(x::FieldPoint, y::FieldPoint, z::FieldPoint, isId::Bool, ec::ECAbstract) =
        if x.field!=y.field || x.field!=z.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, z, isId, convert(ECLD, ec)) end

    ECPointLD(x::FieldPoint, y::FieldPoint, z::FieldPoint, ec::ECAbstract) =
        if x.field!=y.field || x.field!=z.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, z, false, convert(ECLD, ec)) end

    ECPointLD(ec::ECAbstract) =
        new(FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), true, convert(ECLD, ec))
end

function repr(p::ECPointLD)
    return "("*repr(p.x)*", "*repr(p.y)*", "*repr(p.z)*")"
end

function ==(p1::ECPointLD, p2::ECPointLD)
    return p1.isId==p2.isId && p1.ec==p2.ec && p1.x*p2.z==p2.x*p1.z && p1.y*p2.z^2==p2.y*p1.z^2
end

function !=(p1::ECPointLD, p2::ECPointLD)
    return p1.isId!=p2.isId || p1.ec!=p2.ec || p1.x*p2.z!=p2.x*p1.z || p1.y*p2.z^2!=p2.y*p1.z^2
end

function +(p1::ECPointLD, p2::ECPointLD)
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if p1.isId return p2 end
    if p2.isId return p1 end
    if p1==-p2 return ECPointLD(p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 8
    #Mults: 16
    #Sqrs: 4
    #Invs: 0

    A = p1.y*p2.z^2 + p2.y*p1.z^2
    C = p1.x*p2.z + p2.x*p1.z
    D = p2.z*C
    B = p1.z*D

    z3 = B^2
    if iszero(z3) return ECPointLD(p1.ec) end

    x3 = A*(A+B) + B*(C^2 + p1.ec.a*B)

    t = x3*p1.z

    y3 = A*(p1.x*z3 + t) + (t*p1.z + p1.y*z3)*D
    y3 *= D
    if iszero(y3) return ECPointLD(p1.ec) end

    return ECPointLD(x3, y3, z3, p1.ec)
end

function -(p::ECPointLD)
    if p.isId return p end
    return ECPointLD(p.x, p.x+p.y, p.z, p.ec)
end

function double(p::ECPointLD)
    if p.isId return p end
    if p==-p return ECPointLD(p.ec) end

    #Adds: 4
    #Mults: 6
    #Sqrs: 3
    #Invs: 0

    x_2 = p.x^2
    A = x_2 + p.y
    B = p.x*p.z

    z_new = B^2
    if iszero(z_new) return ECPointLD(p.ec) end

    AB = A+B

    x_new = A*AB + p.ec.a*z_new

    y_new = x_2^2 * z_new + x_new*B*AB
    if iszero(y_new) return ECPointLD(p.ec) end

    return ECPointLD(x_new, y_new, z_new, p.ec)
end

function *(p::ECPointLD, n::Integer)
    if p.isId return p end
    if n==0 return ECPointLD(p.ec) end
    if n==1 return p end

    result = ECPointLD(p.ec)
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

function is_valid(p::ECPointLD)
    return p.isId || (p.y^2 + p.x*p.y*p.z == p.x^3*p.z + p.ec.a*p.x^2*p.z^2 + p.ec.b*p.z^4)
end
