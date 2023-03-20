using Documenter
using SingleLinkagePartitions

makedocs(
    sitename = "SingleLinkagePartitions",
    format = Documenter.HTML(),
    modules = [SingleLinkagePartitions]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
