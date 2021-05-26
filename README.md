# ClearSky

<!--- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wordsworthgroup.github.io/ClearSky.jl/stable) --->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wordsworthgroup.github.io/ClearSky.jl/dev)
[![Build Status](https://github.com/wordsworthgroup/ClearSky.jl/workflows/CI/badge.svg)](https://github.com/wordsworthgroup/ClearSky.jl/actions)
[![Codecov](https://img.shields.io/codecov/c/github/wordsworthgroup/ClearSky.jl?logo=codecov)](https://codecov.io/gh/wordsworthgroup/ClearSky.jl)



under development

to use the code, download or clone the repository, take `.jl` out of the folder name, then start Julia with all your threads. For example, if your compute has 8 threads,
```shell
julia --threads 8
```
Tell Julia where you put the code
```julia
push!(LOAD_PATH, "path/to/repo")
```
replacing `path/to/repo` with the path to the folder above the `ClearSky` folder. Then load the code with
```julia
using ClearSky
```
and you should see something like
```
[ Info: Precompiling ClearSky [5964c129-204c-4c32-bd6e-c8dff7ca179b]
```
