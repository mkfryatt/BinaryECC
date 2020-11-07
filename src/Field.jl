import Base: +, -, *, /, ^, ==, repr, inv

struct FieldMismatchException <: Exception end

struct Field
    x::BigInt #the value of this element
    m::Integer #the order of its field
    f::BigInt #the reduction polynomial - x^m
    Field(x::Number, m::Number, f::Number) = Field(convert(BigInt, x), convert(Integer, m), convert(BigInt, f))
    Field(x::BigInt, m::Integer, f::BigInt) = new(x, m, f)
end

function copy(a::Field)
    return Field(a.x, a.m, a.f)
end

function repr(x::Field)
    return repr(x.val)
end

function ==(a::Field, b::Field)
    return a.m==b.m &&  a.x==b.x && a.f==b.f
end

function +(a::Field, b::Field)
    if a.m != b.m || a.f != b.f
        throw(FieldMismatchException)
    end
    return Field(a.x ⊻ b.x, a.m, a.f)
end

function -(a::Field, b::Field)
    return a+b
end

#TODO
function reduce(a::Field)
    return copy(a)
end

#TODO
function *(a::Field, b::Field)
    if a.m != b.m || a.f != b.f
        throw(FieldMismatchException)
    end
    return copy(a)
end

#TODO
function /(a::Field, b::Field)
    if a.m != b.m || a.f != b.f
        throw(FieldMismatchException)
    end
    return copy(a)
end

#TODO
function square(a::Field)
    return copy(a)
end

function ^(a::Field, b::Number)
    result = copy(a)
    while b > 0
        if b&1
            result *= a
        a = square(a)
        b -= 1
    return result
end
