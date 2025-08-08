using Plots;

Xm(ξ) = ξ; # undeformed surface 
xm(ξ) = ξ .+ 0.1 * sin(2π*ξ); # deformed rod 
Fm(ξ) = 1 .+ (π)/5 * cos(2π*ξ); # deformation gradient 
Em(ξ) = 0.5 * (Fm(ξ).^2 .- 1); # Green strain tensor 

ξarr = 0:0.1:1; # array of values to plot 
ξfine = 0:0.01:1; # finer grid for continuous qtys 

# material points of undeformed 
plt = plot(Xm.(ξarr), [4 for _ in ξarr], label = "",
            color = :blue, linewidth = 2,
            xlabel = "X", xguidefont = font(16),
            xticks = [0, 1], yticks = [0],
            xtickfont = font(16), ytickfont = font(16),
            legend = :topright, legendfont = font(12)
            );
scatter!(Xm.(ξarr), [4 for _ in ξarr], label = "Undeformed material points",
            color = :blue, ms = 6)

# material points of deformed 
plot!(xm.(ξarr), [2 for _ in ξarr], label = "",
            color = :red, linewidth = 2);
scatter!(xm.(ξarr), [2 for _ in ξarr], label = "Deformed material points",
            color = :red, ms = 6)

savefig(plt, "1dstrainHalf.png")

# deformation gradient 
plot!(Xm.(ξfine), Fm.(ξfine), label = "Deformation gradient",
        color = :green, linewidth = 3)

# green strain tensor 
plot!(Xm.(ξfine), Em.(ξfine), label = "Strain tensor",
        color = :purple, linewidth = 3)

# zero 
plot!(Xm.(ξfine), [0 for _ in ξfine], label = "",
        color = :black, linestyle = :dash, linewidth = 2)

savefig(plt, "1Dstrain.png");
display(plt);

