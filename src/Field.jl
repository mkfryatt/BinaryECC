import Base: +, -, *, /, ^, ==, !=, repr, inv

struct FieldMismatchException <: Exception end

struct Field
    order::Int16 #the order of the field
    reduction::BigInt #the reduction polynomial - x^m
    Field(m::Number, f::Number) = Field(convert(Int16, m), convert(BigInt, f))
end

struct FieldPoint
    x::BigInt
    field::Field
    FieldPoint(x::Number, m::Number, f::Number) = FieldPoint(convert(BigInt, x), Field(m, f))
    FieldPoint(x::BigInt, field:Field) = new(x, field)
end

function copy(f::Field)
    return Field(f.order, f.reduction)
end

function copy(a::FieldPoint)
    return FieldPoint(a.x, a.field)
end

function repr(f::Field)
    return "Order: "*repr(f.order)*"\nReduction: "*repr(reduction)
end

function repr(a::FieldPoint)
    #TODO
    return "TODO"
end

function ==(a::Field, b::Field)
    return a.order==b.order && a.reduction==b.reduction
end

function ==(a::FieldPoint, b::FieldPoint)
    return a.x==b.x && a.field==b.field
end

function !=(a::Field, b::Field)
    return a.order!=b.order || a.reduction!=b.reduction
end

function !=(a::FieldPoint, b::FieldPoint)
    return a.x!=b.x || a.field!=b.field
end

function +(a::FieldPoint, b::FieldPoint)
    if a.field!=b.field throw(FieldMismatchException) end
    #TODO check
    return FieldPoint(a.x ⊻ b.x, a.field)
end

function -(a::FieldPoint, b::FieldPoint)
    return a+b
end

function reduce(a::FieldPoint)
    #TODO
    return copy(a)
end

function *(a::FieldPoint, b::FieldPoint)
    if a.field!=b.field throw(FieldMismatchException) end
    #TODO
    return copy(a)
end

function /(a::FieldPoint, b::FieldPoint)
    if a.field!=b.field throw(FieldMismatchException) end
    #TODO
    return copy(a)
end

function square(a::FieldPoint)
    #TODO
    return copy(a)
end

function ^(a::FieldPoint, b::Number)
    #TODO
end

function random(f::Field)
    #TODO
    return FieldPoint(0, f)
end
