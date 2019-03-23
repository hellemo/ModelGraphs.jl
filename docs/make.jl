using Documenter, AlgebraicGraphs

makedocs(;
    modules=[AlgebraicGraphs],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jalving/AlgebraicGraphs.jl/blob/{commit}{path}#L{line}",
    sitename="AlgebraicGraphs.jl",
    authors="Jordan Jalving",
    assets=[],
)

deploydocs(;
    repo="github.com/jalving/AlgebraicGraphs.jl",
)
