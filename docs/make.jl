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
        "Absorption Data" => "absorption_data.md",
        "Line Shapes" => "line_shapes.md",
        "Gas Objects" => "gas_objects.md",
        "Atmospheric Profiles" => "atmospheric_profiles.md"
    ],
)

deploydocs(;
    repo="github.com/wordsworthgroup/ClearSky.jl",
    devbranch="main",
    versions=["stable"=>"v^"]
)
