using Documenter
using ClearSky

DocMeta.setdocmeta!(ClearSky, :DocTestSetup, :(using ClearSky); recursive=true)

makedocs(;
    modules=[ClearSky],
    authors="Mark Baum",
    repo="https://github.com/markmbaum/ClearSky.jl/blob/{commit}{path}#{line}",
    sitename="ClearSky.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://markmbaum.github.io/ClearSky.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Absorption Data" => "absorption_data.md",
        "Line Shapes" => "line_shapes.md",
        "Gas Objects" => "gas_objects.md",
        "Atmospheric Profiles" => "atmospheric_profiles.md",
        "Radiation" => "radiation.md",
        "Modeling" => "modeling.md",
        "Orbits & Insolation" => "orbits_insolation.md"
    ],
)

deploydocs(;
    repo="github.com/markmbaum/ClearSky.jl",
    devbranch="main",
    versions=nothing
)
