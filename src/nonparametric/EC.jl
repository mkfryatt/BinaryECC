import Base: -, *, ==, repr

struct ECMismatchException <: Exception end

abstract type AbstractECPoint end

struct EC
    a::FieldPoint
    b::FieldPoint
    EC(a::FieldPoint, b::FieldPoint) =
        if a.field!=b.field throw(FieldMismatchException())
        else new(a, b) end
end

function repr(ec::EC)
    return "E: y² + xy = x³ + "*repr(ec.a.x)*"x² + "*repr(ec.b.x)
end

function ==(ec1::EC, ec2::EC)
    return ec1.a==ec2.a && ec1.b==ec2.b && ec1.a.field==ec2.a.field
end

function -(p1::AbstractECPoint, p2::AbstractECPoint)
    return p1 + (-p2)
end

function *(n::Integer, p::AbstractECPoint)
    return p*n
end

#scalar multiplication
#performs a fixed sequence of curve and field ops
function montpow(p::AbstractECPoint, n::Integer)
    R0 = p
    R1 = double(p)
    for i in (bits(n)-2):-1:0
        bit = (n>>>i)&1
        if bit==0
            R1, R0 = R0+R1, double(R0)
        else
            R0, R1 = R0+R1, double(R1)
        end
    end
    return R0
end
