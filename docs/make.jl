using Documenter, Multitaper

makedocs(
    modules = [Multitaper],
    sitename = "Multitaper",
    repo = "https://github.com/lootie/multitaper.jl/src/{commit}{path}#lines-{line}",
)

deploydocs(
    repo = "github.com/lootie/Multitaper.jl.git",
)
