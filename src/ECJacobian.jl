import Base: +, -, *, ==, repr, isvalid, iszero

struct ECPointJacobian{D,R} <: AbstractECPoint
    x::FieldPoint{D,R}
    y::FieldPoint{D,R}
    z::FieldPoint{D,R}
    ec::EC{D,R}

    ECPointJacobian{D,R}(x::FieldPoint{D,R}, y::FieldPoint{D,R}, z::FieldPoint{D,R}, ec::EC{D,R}) where D where R =
        new(x, y, z, ec)

    ECPointJacobian{D,R}(ec::EC{D,R}) where D where R =
        new(FieldPoint{D,R}(0), FieldPoint{D,R}(0), FieldPoint{D,R}(0), ec)
end

function repr(p::ECPointJacobian)
    return "("*repr(p.x)*", "*repr(p.y)*", "*repr(p.z)*")"
end

function ==(p1::ECPointJacobian{D,R}, p2::ECPointJacobian{D,R}) where D where R
    z1_2 = p1.z^2
    z2_2 = p2.z^2
    return iszero(p1)==iszero(p2) && p1.ec==p2.ec && p1.x*z2_2==p2.x*z1_2 && p1.y*p2.z*z2_2==p2.y*p1.z*z1_2
end

function +(p1::ECPointJacobian{D,R}, p2::ECPointJacobian{D,R}) where D where R
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
    if p1==-p2 return ECPointJacobian{D,R}(p1.ec) end
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

    x3 = A*(A+B) + C^3 + p1.ec.a*B^2
    x3 *= z1_2^2

    z3_2 = z3^2
    t = x3*z1_2
    y3 = A*(p1.x*z3_2 + t)
    y3 += C*p2.z*(t*p1.z + p1.y*z3_2)

    return ECPointJacobian{D,R}(x3, y3, z3, p1.ec)
end

function -(p::ECPointJacobian{D,R}) where D where R
    if isero(p) return p end
    return ECPointJacobian{D,R}(p.x, p.x+p.y, p.z, p.ec)
end

function double(p::ECPointJacobian{D,R}) where D where R
    if iszero(p) return p end
    if p==-p return ECPointJacobian{D,R}(p.ec) end

    #Adds: 5
    #Mults: 10
    #Sqrs: 6
    #Invs: 0

    z_2 = p.z^2
    A = p.x^2 + p.y*p.z
    B = p.x*z_2

    z_4 = z_2^2

    z_new = B*z_4

    x_new = A*(A+B)
    x_new += p.ec.a * B^2
    x_new *= z_4^2

    y_new = B * (p.x*z_new)^2
    y_new += (A+B)*x_new*z4

    return ECPointJacobian{D,R}(x_new, y_new, z_new, p.ec)
end

function *(p::ECPointJacobian{D,R}, n::Integer) where D where R
    if iszero(p) return p end
    if n==0 return ECPointJacobian{D,R}(p.ec) end
    if n==1 return p end

    result = ECPointJacobian{D,R}(p.ec)
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

function isvalid(p::ECPointJacobian)
    return iszero(p) || (p.y^2 + p.x*p.y*p.z == p.x^3 + p.ec.a*p.x^2*p.z^2 + p.ec.b*p.z^6)
end

function iszero(p::ECPointJacobian)
    return iszero(p.z)
end
