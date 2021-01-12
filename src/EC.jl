struct ECMismatchException <: Exception end

abstract type AbstractECPoint end

struct EC{D,R}
    a::FieldPoint{D,R}
    b::FieldPoint{D,R}
end

function repr(ec::EC)
    return "E: y² + xy = x³ + "*repr(ec.a.value)*"x² + "*repr(ec.b.value)
end

function ==(ec1::EC, ec2::EC)
    return ec1.a==ec2.a && ec1.b==ec2.b
end

function -(p1::AbstractECPoint, p2::AbstractECPoint)
    return p1 + (-p2)
end

function *(n::Integer, p::AbstractECPoint)
    return p*n
end
