using Literate

function replace_includes(str)

    included = ["./snippet/plot_partition.jl";]
    path = "./"

    for ex in included
        content = read(path * ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    return str
end

Literate.markdown(
    "chaining.jl",
    "../docs/src/generated/";
    execute = true,
    preprocess = replace_includes,
)
