# Generate documentation with this command:
# (cd docs && julia --color=yes make.jl)

push!(LOAD_PATH, "..")

using Documenter
using GridFunctions

makedocs(; sitename="GridFunctions", format=Documenter.HTML(),
         modules=[GridFunctions])

deploydocs(; repo="github.com/eschnett/GridFunctions.jl.git", devbranch="main",
           push_preview=true)
