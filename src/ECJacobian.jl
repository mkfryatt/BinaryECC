import Base: +, -, *, ==, !=, repr, ceil

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

    A = p1.y*p2.z^3 + p2.y*p1.z^3
    B = p1.x*p1.z*p2.z^3 + p2.x*p1.z^3*p2.z

    z_new = B * p1.z^3 * p2.z
    if iszero(z_new) return ECPointJacobian(p1.ec) end

    B2 = B^2

    z1z2_2 = (p1.z*p2.z)^2

    x_new = A*z1z2_2*(A+B)
    x_new += B2*(p1.x*p2.z^2 + p2.x*p1.z^2)
    x_new += p1.ec.a*B2*z1z2_2
    x_new *= p1.z^4

    y_new = A * p1.x * p1.z * z_new^2
    y_new += (A+B) * x_new * p1.z^3
    y_new += B * p1.y * z_new^2
    y_new *= p2.z

    return ECPointJacobian(x_new, y_new, z_new, p1.ec)
end

function -(p::ECPointJacobian)
    if p.isId return p end
    return ECPointJacobian(p.x, p.x+p.y, p.z, p.ec)
end

function double(p::ECPointJacobian)
    if p.isId return p end
    if p==-p return ECPointJacobian(p.ec) end

    A = p.x^2 + p.y*p.z
    B = p.x*p.z^2

    z4 = p.z^4

    z_new = B*z4
    if iszero(z_new) return ECPointJacobian(p.ec) end

    x_new = A*(A+B)
    x_new += p.ec.a * B^2
    x_new *= z4^2

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
