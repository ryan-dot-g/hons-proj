using Plots, LaTeXStrings, Subscripts, Measures;

do1 = false

if do1
    C = [1.0, 0.5, 0.3];
    nsol = length(C);
    evec(x, i) = sin(π * x * i);
    soln(x) = sum(C[i] * evec(x, i) for i in 1:nsol);

    X = 0:0.01:1;
    p = plot(X, soln.(X), lw = 2, title = "Eigendecomposition of guitar",
                label = "String vibration", xlabel = "x", ylabel = "Amplitude",
                xticks = [], yticks = [], margin = 4mm,
                titlefontsize = 18, legendfontsize = 14, guidefontsize = 14);
    for i in 1:nsol
        plot!(X, evec.(X, i), ls = :dash, label = "f"*sub(string(i)))
    end
    display(p); savefig(p, "guitarEcomp.pdf")
end

### ----------- EIGENDECOMPOSITION OF HYDRA? -------------- ###

θ = -π/2:0.02:π/2;
r = 1;
Xunpet = r*cos.(θ); Yunpet = r*sin.(θ);
Yfull = [Yunpet, Yunpet];
full(b) = [-b, b];

xfund(ξ) = r * ( 1 - 4/π^2 * ξ^2 + 0.2 * (ξ^2 - π^2/4) * ξ );
yfund(ξ) = r * sin(ξ);

x1(ξ) = r*(1-4/π^2 * ξ^2);
x2(ξ) = r*0.2 * (ξ^2 - π^2/4) * ξ;

p = plot(full(xfund.(θ)), Yfull, lw = 2, title = "Eigendecomposition of hydra?",
            color = 1,
            label = ["Hydra" ""], xlabel = "x", ylabel = "y",
            xticks = [], yticks = [], margin = 4mm,
            aspect_ratio = :equal,
            titlefontsize = 18, legendfontsize = 14, guidefontsize = 14);
plot!(full(Xunpet), Yfull, color = :black, lw = 1, label = ["Steady-state" ""])

x1good = (x1.(θ) .- Xunpet)*4 .+ Xunpet;
plot!(full(x1good), Yfull, lw = 1, ls = :dash, color = 2, label = [L"Z_1" ""]);

x2good = (x2.(θ) .- Xunpet)*0.3 .+ Xunpet;
plot!(full(x2good), Yfull, lw = 1, ls = :dash, color = 3, label = [L"Z_2" ""]);

display(p); savefig(p, "hydraEcomp.pdf")