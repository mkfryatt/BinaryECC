import Base: +, -, *, ==, repr, isvalid, iszero

struct ECPointLD <: AbstractECPoint
    x::FieldPoint
    y::FieldPoint
    z::FieldPoint
    ec::EC

    ECPointLD(x::FieldPoint, y::FieldPoint, z::FieldPoint, ec::EC) =
        if x.field!=y.field || x.field!=z.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, z, ec) end

    ECPointLD(x::FieldPoint, y::FieldPoint, z::FieldPoint, ec::EC) =
        if x.field!=y.field || x.field!=z.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, z, ec) end

    ECPointLD(ec::EC) =
        new(FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), ec)
end

function repr(p::ECPointLD)
    return "("*repr(p.x)*", "*repr(p.y)*", "*repr(p.z)*")"
end

function ==(p1::ECPointLD, p2::ECPointLD)
    return iszero(p1)==iszero(p2) && p1.ec==p2.ec && p1.x*p2.z==p2.x*p1.z && p1.y*p2.z^2==p2.y*p1.z^2
end

function +(p1::ECPointLD, p2::ECPointLD)
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
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
    if iszero(p) return p end
    return ECPointLD(p.x, p.x+p.y, p.z, p.ec)
end

function double(p::ECPointLD)
    if iszero(p) return p end
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
    if iszero(p) return p end
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

function isvalid(p::ECPointLD)
    return iszero(p) || (p.y^2 + p.x*p.y*p.z == p.x^3*p.z + p.ec.a*p.x^2*p.z^2 + p.ec.b*p.z^4)
end

function iszero(p::ECPointLD)
    return iszero(p.y) && iszero(p.z)
end
