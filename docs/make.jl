using Documenter
using SingleLinkagePartitions


# # # local.
# makedocs(
#     sitename = "SingleLinkagePartitions",
#     modules = [SingleLinkagePartitions],
#     #format = Documenter.HTML(),
#     pages = [
#         "Overview" => "index.md",
#         "Public API" => "api.md",
#         "Demo: chaining & remedy" => "generated/chaining.md",
#     ],
# )

# github.
makedocs(
    sitename="SingleLinkagePartitions.jl",
    modules=[SingleLinkagePartitions],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
    pages = [
        "Overview" => "index.md",
        "Public API" => "api.md",
        "Demo: image analysis" => "demo.md",
    ],
)
deploydocs(
    repo = "github.com/RoyCCWang/SingleLinkagePartitions.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)