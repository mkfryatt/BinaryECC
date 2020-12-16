import Base: +, -, *, ==, repr, isvalid, iszero
using Base: ceil

struct ECPointAffine{D,R} <: AbstractECPoint
    x::FieldPoint{D,R}
    y::FieldPoint{D,R}
    isId::Bool
    ec::EC{D,R}

    ECPointAffine{D,R}(x::FieldPoint{D,R}, y::FieldPoint{D,R}, isId::Bool, ec::EC{D,R}) where D where R =
        new(x, y, isId, ec)

    ECPointAffine{D,R}(x::FieldPoint{D,R}, y::FieldPoint{D,R}, ec::EC{D,R}) where D where R =
        new(x, y, false, ec)

    ECPointAffine{D,R}(ec::EC{D,R}) where D where R =
        new(FieldPoint{D,R}(0), FieldPoint{D,R}(0), true, ec)
end

#sec1v2 2.3.4
function ECPointAffine{D,R}(s::String, ec::EC{D,R}) where D where R
    s = replace(s, " " => "")

    #point is id
    if s=="00" return ECPointAffine(ec) end

    #input string is not of the uncompressed format specified by sec1v2
    if length(s) != 4*ceil(Int16, D / 8)+2
        throw(ArgumentError("Octet string is of the wrong length for this curve."))
    end
    if s[1:2]!="04"
        throw(ArgumentError("Octet string must start with '04'."))
    end

    #TODO handle compressed format

    x = FieldPoint{D,R}(s[3:4+2*floor(Int16, D / 8)])
    y = FieldPoint{D,R}(s[5+2*floor(Int16, D / 8):6+4*floor(Int16, D / 8)])

    return ECPointAffine{D,R}(x, y, ec)
end

function repr(p::ECPointAffine)
    return "("*repr(p.x)*", "*repr(p.y)*")"
end

function ==(p1::ECPointAffine{D,R}, p2::ECPointAffine{D,R}) where D where R
    return p1.x==p2.x && p1.y==p2.y && iszero(p1)==iszero(p2) && p1.ec==p2.ec
end

function +(p1::ECPointAffine{D,R}, p2::ECPointAffine{D,R}) where D where R
    if p1.ec!=p2.ec throw(ECMismatchException()) end
    if iszero(p1) return p2 end
    if iszero(p2) return p1 end
    if p1==-p2 return ECPointAffine{D,R}(p1.ec) end
    if p1==p2 return double(p1) end

    #Adds: 8
    #Mults: 2
    #Sqrs: 1
    #Invs: 1
    x1x2 = p1.x + p2.x
    lambda = (p1.y+p2.y) / x1x2
    x3 = lambda^2 + lambda + x1x2 + p1.ec.a
    y3 = lambda*(p1.x+x_new) + x3 + p1.y
    return ECPointAffine{D,R}(x3, y3, p1.ec)
end

function -(p::ECPointAffine{D,R}) where D where R
    if iszero(p) return p end
    return ECPointAffine{D,R}(p.x, p.x+p.y, p.ec)
end

function double(p::ECPointAffine{D,R}) where D where R
    if iszero(p) return p end
    if p==-p return ECPointAffine{D,R}(p.ec) end

    #Adds: 5
    #Mults: 2
    #Sqrs: 2
    #Invs: 1
    lambda = p.x + (p.y / p.x)
    x_new = lambda^2 + lambda + p.ec.a
    y_new = p.x^2 + lambda*x_new + x_new
    return ECPointAffine{D,R}(x_new, y_new, p.ec)
end

function *(p::ECPointAffine{D,R}, n::Integer) where D where R
    if iszero(p) return p end
    if n==0 return ECPointAffine{D,R}(p.ec) end
    if n==1 return p end

    result = ECPointAffine{D,R}(p.ec)
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
