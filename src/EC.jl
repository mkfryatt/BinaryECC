struct ECMismatchException <: Exception end

abstract type ECAbstract end
abstract type ECPointAbstract end

function ==(ec1::ECAbstract, ec2::ECAbstract)
    return ec1.a==ec2.a && ec1.b==ec2.b && ec1.a.field==ec2.a.field
end

function !=(ec1::ECAbstract, ec2::ECAbstract)
    return ec1.a!=ec2.a || ec1.b!=ec2.b || ec1.a.field!=ec2.a.field
end

function -(p1::ECPointAbstract, p2::ECPointAbstract)
    return p1 + (-p2)
end

function *(n::Integer, p::ECPointAbstract)
    return p*n
end
