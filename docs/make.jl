using Documenter, Multitaper

makedocs(
    modules = [Multitaper],
    sitename = "Multitaper",
    repo = "https://github.com/lootie/multitaper.jl/src/{commit}{path}#lines-{line}",
)

withenv(
    "TRAVIS_REPO_SLUG" => "github.com/lootie/lootie.bitbucket.io",
    "TRAVIS_PULL_REQUEST" => get(ENV, "BITBUCKET_BRANCH", nothing),
    "TRAVIS_BRANCH" => get(ENV, "BITBUCKET_BRANCH", nothing),
    "TRAVIS_TAG" => get(ENV, "BITBUCKET_TAG", nothing),
    "TRAVIS_PULL_REQUEST" => ("BITBUCKET_PR_ID" in keys(ENV)) ? "true" : "false",
) do
    deploydocs(
        repo = "github.com/lootie/lootie.bitbucket.io.git",
        branch = "gh-pages",
        dirname = "Multitaper.jl",
        devbranch = "master"
    )
end


