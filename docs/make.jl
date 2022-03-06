using NumericalSolutionOfPartialDifferentialEquations
using Documenter

DocMeta.setdocmeta!(NumericalSolutionOfPartialDifferentialEquations, :DocTestSetup, :(using NumericalSolutionOfPartialDifferentialEquations); recursive=true)

makedocs(;
    modules=[NumericalSolutionOfPartialDifferentialEquations],
    authors="Quejiahao <quejiahao@live.com> and contributors",
    repo="https://github.com/Quejiahao/NumericalSolutionOfPartialDifferentialEquations.jl/blob/{commit}{path}#{line}",
    sitename="NumericalSolutionOfPartialDifferentialEquations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Quejiahao.github.io/NumericalSolutionOfPartialDifferentialEquations.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Quejiahao/NumericalSolutionOfPartialDifferentialEquations.jl",
    devbranch="main",
)
