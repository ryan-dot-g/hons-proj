using Plots; 

cmap = cgrad(:viridis);

X = 1:0.001:2; # x-pos, aesthetics only
Y1 = X.*0 .+ 2.0; # y-pos, aesthetics only
Y2 = X.*0 .+ 1.0; 

C1 = 1 .+ 0.25*exp.(-100*(X.-1.5).^2); # sharp peak
C2 = 1 .- (X.-1).*(X.-2); # shallow peak 

plt = plot(X, Y1, line_z = C1, c = cmap,
            linewidth = 50,
            colorbar = false, legend = false,
            grid = false, framestyle = :none,
            xaxis = false, yaxis = false,
            dpi = 500, 
            );
plot!(X, Y2, line_z = C2, c = cmap,
            linewidth = 50);
plot!(X, Y1.+0.8, linewidth = 0) # extra to force the lines closer together 
plot!(X, Y2.-0.8, linewidth = 0)

savefig(plt, "activatorInhibitorSTEP1.png");
display(plt);

