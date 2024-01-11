
partition, title_string, X, fig_num, marker_size = ARGS

PLT.figure(fig_num)
fig_num += 1

for k in eachindex(partition)
    X_part = X[partition[k]]
    plot_x = collect( X_part[n][begin] for n in eachindex(X_part) )
    plot_y = collect( X_part[n][begin+1] for n in eachindex(X_part) )

    if length(X_part) < 2
        PLT.scatter(plot_x, plot_y, s = marker_size, label = "part $k", marker = "x")
    else
        PLT.scatter(plot_x, plot_y, s = marker_size, label = "part $k")
    end
end

PLT.axis("scaled")
PLT.title(title_string)
PLT.legend()
PLT.gcf()