
# Lifted very closely and with a great deal of gratitude from: 
#
#   https://bitbucket.org/mortenpi/documenterexample.jl/src/master/docs/make.jl
#

using Documenter, Multitaper

makedocs(
    modules = [Multitaper],
    sitename = "Multitaper",
    repo = "https://bitbucket.org/clhaley/multitaper.jl/src/{commit}{path}#lines-{line}",
)

withenv(
    "TRAVIS_REPO_SLUG" => "bitbucket.org/clhaley/clhaley.bitbucket.io",
    "TRAVIS_PULL_REQUEST" => get(ENV, "BITBUCKET_BRANCH", nothing),
    "TRAVIS_BRANCH" => get(ENV, "BITBUCKET_BRANCH", nothing),
    "TRAVIS_TAG" => get(ENV, "BITBUCKET_TAG", nothing),
    "TRAVIS_PULL_REQUEST" => ("BITBUCKET_PR_ID" in keys(ENV)) ? "true" : "false",
) do
    deploydocs(
        repo = "bitbucket.org/clhaley/clhaley.bitbucket.io.git",
        branch = "gh-pages",
        dirname = "Multitaper.jl",
        devbranch = "master"
    )
end


