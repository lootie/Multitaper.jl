using Documenter
using Multitaper

makedocs(
    sitename = "Multitaper",
    format = Documenter.HTML(),
    modules = [Multitaper],
    repo = "https://bitbucket.org/clhaley/multitaper.jl/{commit}{path}#{line}"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

deploydocs(
    repo = "git@bitbucket.org:clhaley/multitaper.jl.git",
)

# deploydocs(
#    repo = "git@bitbucket.org:lootie/clhaley.bitbucket.io.git",
#    branch = "master", 
#    root = "Multitaper.jl"
#)

