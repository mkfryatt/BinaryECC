push!(LOAD_PATH,"../src/")
using Documenter, BinaryECC

makedocs(
    sitename = "BinaryECC Documentation",
    modules = [BinaryECC],
    pages = [
        "Home" => "index.md",
        "Binary Fields" => "field.md",
        "Elliptic Curves" => "ec.md",
        "Cryptographic Primitives" => "crypto.md"
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)
