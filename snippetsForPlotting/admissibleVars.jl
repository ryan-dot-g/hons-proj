using Plots, LaTeXStrings;

a = 1.0; b = 2.0; # start/endpoints of function 
X = a:0.01:b; 

f(x) = cos.(5*x) .+ 0.5*x.^2 .+ 2*x; # sufficiently interesting fn 
df1(x) = cos.(2π*x) .- 1; # sufficiently interesting variation 
df2(x) = 10*(x.-1)*(x.-2)*(x.-1.5); # sufficiently interesting varn 2 

# plt = plot(X, f.(X), label = L"f(x)",
#             linewidth = 3,
#             xlabel = "x", xguidefont = font(16),
#             ylabel = "y", yguidefont = font(16),
#             xticks = ([a,b], [L"a_x", L"b_x"]), 
#             yticks = ([f(a), f(b)], [L"a_y", L"b_y"]),
#             xtickfont = font(16), ytickfont = font(16),
#             legend = :bottomright, legendfont = font(12));
# plot!(X, f.(X) .+ 0.1*df1.(X), label = L"f(x) + \epsilon\cdot\delta f_1(x)",
#         linewidth = 2, linestyle = :dash)
# plot!(X, f.(X) .+ 0.4*df2.(X), label = L"f(x) + \epsilon\cdot\delta f_2(x)",
#         linewidth = 2, linestyle = :dash)


plt = plot(X, f.(X) .+ 0.4*df2.(X), label = L"f(x) + \epsilon\cdot\delta f_1(x)",
            color = :red, linewidth = 2, linestyle = :dash,
            xlabel = "x", xguidefont = font(16),
            ylabel = "y", yguidefont = font(16),
            xticks = ([a,b], [L"a_x", L"b_x"]), 
            yticks = ([f(a), f(b)], [L"a_y", L"b_y"]),
            xtickfont = font(16), ytickfont = font(16),
            legend = :bottomright, legendfont = font(12));
plot!(X, f.(X), label = L"f(x)",
        color = :blue, linewidth = 3)
plot!(X, f.(X) .+ 0.1*df1.(X), label = L"f(x) + \epsilon\cdot\delta f_2(x)",
        color = :green, linewidth = 2, linestyle = :dash)

savefig(plt, "admissibleVars.png");
display(plt);

