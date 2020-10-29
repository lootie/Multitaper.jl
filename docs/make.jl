using Documenter
using Multitaper

makedocs(
    sitename = "Multitaper",
    format = Documenter.HTML(),
    modules = [Multitaper]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#


