using Documenter
using ClearSky

DocMeta.setdocmeta!(ClearSky, :DocTestSetup, :(using ClearSky); recursive=true)

makedocs(;
    modules=[ClearSky],
    authors="Mark Baum",
    repo="https://github.com/wordsworthgroup/ClearSky.jl/blob/{commit}{path}#{line}",
    sitename="ClearSky.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wordsworthgroup.github.io/ClearSky.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Reading Absorption Data" => "reading_absorption_data.md",
        "Computing Line Shapes" => "computing_line_shapes.md",
        
    ],
)

deploydocs(;
    repo="github.com/wordsworthgroup/ClearSky.jl",
)
