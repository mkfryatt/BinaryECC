import Base: +, -, *, ==, !=, repr, convert

struct ECJacobian <: ECAbstract
    a::FieldPoint
    b::FieldPoint
    ECJacobian(a::FieldPoint, b::FieldPoint) =
        if a.field!=b.field throw(FieldMismatchException())
        else new(a, b) end
end

function convert(::Type{ECJacobian}, ec::ECAbstract)
    return ECJacobian(ec.a, ec.b)
end

function repr(ec::ECJacobian)
    return "E: y² + xyz = x³ + "*repr(ec.a.x)*"x²z² + "*repr(ec.b.x)*"z⁶"
end

struct ECPointJacobian <: ECPointAbstract
    x::FieldPoint
    y::FieldPoint
    z::FieldPoint
    isId::Bool
    ec::ECAbstract

    ECPointJacobian(x::FieldPoint, y::FieldPoint, z::FieldPoint, isId::Bool, ec::ECAbstract) =
        if x.field!=y.field || x.field!=z.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, z, isId, convert(ECJacobian, ec)) end

    ECPointJacobian(x::FieldPoint, y::FieldPoint, z::FieldPoint, ec::ECAbstract) =
        if x.field!=y.field || x.field!=z.field || x.field!=ec.a.field throw(FieldMismatchException())
        else new(x, y, z, false, convert(ECJacobian, ec)) end

    ECPointJacobian(ec::ECAbstract) =
        new(FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), FieldPoint(0, ec.a.field), true, convert(ECJacobian, ec))
end

function repr(p::ECPointJacobian)
    return "("*repr(p.x)*", "*repr(p.y)*", "*repr(p.z)*")"
end

function ==(p1::ECPointJacobian, p2::ECPointJacobian)
    return p1.isId==p2.isId && p1.ec==p2.ec && p1.x*p2.z^2==p2.x*p1.z^2 && p1.y*p2.z^3==p2.y*p1.z^3
end

function !=(p1::ECPointJacobian, p2::ECPointJacobian)
    return p1.isId!=p2.isId || p1.ec!=p2.ec || p1.x*p2.z^2!=p2.x*p1.z^2 || p1.y*p2.z^3!=p2.y*p1.z^3
end

function +(p1::ECPointJacobian, p2::ECPointJacobian)
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if p1.isId return p2 end
    if p2.isId return p1 end
    if p1==-p2 return ECPointJacobian(p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 8
    #Mults: 20
    #Sqrs: 6
    #Invs: 0

    z1_2 = p1.z^2
    z2_2 = p2.z^2
    A = p1.y*p2.z*z2_2 + p2.y*p1.z*z1_2
    C = p1.x*z2_2 + p2.x*z1_2
    B = p1.z * p2.z * C

    z3 = B * z1_2
    if iszero(z3) return ECPointJacobian(p1.ec) end

    x3 = A*(A+B) + C^3 + p1.ec.a*B^2
    x3 *= z1_2^2

    z3_2 = z3^2
    t = x3*z1_2
    y3 = A*(p1.x*z3_2 + t)
    y3 += C*p2.z*(t*p1.z + p1.y*z3_2)

    return ECPointJacobian(x3, y3, z3, p1.ec)
end

function -(p::ECPointJacobian)
    if p.isId return p end
    return ECPointJacobian(p.x, p.x+p.y, p.z, p.ec)
end

function double(p::ECPointJacobian)
    if p.isId return p end
    if p==-p return ECPointJacobian(p.ec) end

    #Adds: 5
    #Mults: 10
    #Sqrs: 6
    #Invs: 0

    z_2 = p.z^2
    A = p.x^2 + p.y*p.z
    B = p.x*z_2

    z_4 = z_2^2

    z_new = B*z_4
    if iszero(z_new) return ECPointJacobian(p.ec) end

    x_new = A*(A+B)
    x_new += p.ec.a * B^2
    x_new *= z_4^2

    y_new = B * (p.x*z_new)^2
    y_new += (A+B)*x_new*z4

    return ECPointJacobian(x_new, y_new, z_new, p.ec)
end

function *(p::ECPointJacobian, n::Integer)
    if p.isId return p end
    if n==0 return ECPointJacobian(p.ec) end
    if n==1 return p end

    result = ECPointJacobian(p.ec)
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

function is_valid(p::ECPointJacobian)
    return p.isId || (p.y^2 + p.x*p.y*p.z == p.x^3 + p.ec.a*p.x^2*p.z^2 + p.ec.b*p.z^6)
end
