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