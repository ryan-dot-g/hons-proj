### TODO 
# Fix color mappings
# convert operators to sparse?

using LinearAlgebra, DifferentialEquations, JuMP, Ipopt,
        Plots, LaTeXStrings; 
println("running...");

############ -------------------------------------------------- ############
############ ------------- PHYSICAL PARAMETERS ---------------- ############
############ -------------------------------------------------- ############
R = 130; # radius of force-free reference sphere
r = 150; # TODO radius of initial condition sphere
E0 = 440; # base Young's modulus
h = 20 # bilayer thickness
b = 5; # exponential decrease of stiffness
κ = 1; # TODO morphogen stress upregulation coefficient 
ζ = 0.2; # morphogen decay rate 
D = 0.01; # morphogen diffusion coefficient 

E(phi) = E0 * h * exp(-b*phi); # elasticity function
f(ϵ) = κ * ϵ; # strain-dependent morphogen expression function 

############ -------------------------------------------------- ############
############ ------------- NUMERICAL PARAMETERS --------------- ############
############ -------------------------------------------------- ############
Ndisc = 100; # number of discretisation points on s 
smin = -pi/2; smax = pi/2; # bounds of s values

############ -------------------------------------------------- ############
############ ---------------- NUMERICAL SETUP ----------------- ############
############ -------------------------------------------------- ############
Si = range(smin, smax, Ndisc); # grid of s values 
ds = Si[2] - Si[1]; 

############ -------------------------------------------------- ############
############ ------------- DISCRETISED OPERATORS -------------- ############
############ -------------------------------------------------- ############
# first deriv 
P = Tridiagonal(fill(-1.0, Ndisc-1), fill(0.0, Ndisc), fill(1.0, Ndisc-1)); # main section 
P = Matrix(P); # convert out of sparse to modify elements 
P[1, 1] = -3; P[1, 2] = 4; P[1, 3] = -1; # forward diff top row 
P[Ndisc, Ndisc] = 3; P[Ndisc, Ndisc-1] = -4; P[Ndisc, Ndisc-2] = 1; # backward diff bottom row 
P /= (2*ds); # scale 

# second deriv
Q = Tridiagonal(fill(1.0, Ndisc-1), fill(-2.0, Ndisc), fill(1.0, Ndisc-1)); # main section 
Q = Matrix(Q); # convert out of sparse to modify elements 
Q[1, 1] = 2; Q[1, 2] = -5; Q[1, 3] = 4; Q[1, 4] = -1; # forward diff top row 
Q[Ndisc, Ndisc] = 2; Q[Ndisc, Ndisc-1] = -5; Q[Ndisc, Ndisc-2] = 4;  Q[Ndisc, Ndisc-3] = -1; # backwards diff bottom row 
Q[1,:] /= ds; Q[Ndisc, :] /= ds; # add 1/ds factor to f/b diffs on top/bottom row 
Q /= ds^2; # scale 

############ -------------------------------------------------- ############
############ ------------------ FUNCTIONALS ------------------- ############
############ -------------------------------------------------- ############



############ -------------------------------------------------- ############
############ ------- OTHER AUXILIARY / USEFUL FUNCTIONS ------- ############
############ -------------------------------------------------- ############
ϕ0 = 1/ζ * f(0.5 * ((r^2 - R^2)/R^2)^2); # steady-state morphogen scalar
Phi0 = Si .* 0 .+ ϕ0; # steady-state morphogen over the full grid 

X0 = r * cos.(Si); Y0 = r * sin.(Si); # initial condition shape 
Xundef = R * cos.(Si); Yundef = R * sin.(Si); # undeformed shape

# 2nd-order approximations of r cos(s), r sin(s) to use as initial guesses for 
#   x(s), y(s) at steady-state
x0quad(s) = r * ( 1 - 4/pi^2 * s^2 + 0.2 * (s^2 - pi^2/4) * s );
y0quad(s) = r * (s);

############ -------------------------------------------------- ############
############ ------------- PRESAVING CALCULATIONS ------------- ############
############ -------------------------------------------------- ############
CosS = cos.(Si); Sins = sin.(Si); # presaving trigs 
TanS = 0 * CosS; TanS[2:end-1] = tan(Si[2:end-1]); # presaving Tan, with 0 on bdy to avoid undefineds


############ -------------------------------------------------- ############
############ ------------ VISUALISATION FUNCTIONS ------------- ############
############ -------------------------------------------------- ############
function visualise(X, Y, Phi, titleTxt = false)
    # takes a surface discretised parameterisation and morphogen concentration
    # displays a plot of the deformed hydra, colored by the morphogen concentration 
    # optionally adds a provided title 
    # returns the plot object for displaying and saving

    cmap = cgrad(:viridis); # produce colormap 
    nphi = (Phi .- minimum(Phi)) ./ (maximum(Phi) - minimum(Phi)); # normalise for color 

    # first, plot actual deformed shape 
    plt = plot(X, Y, lw = 2, color = :black, label = "Hydra shape");
    plot!(-X, Y, lw = 2, color = :black, label = "");

    # next, plot initial shape 
    plot!(X0, Y0, lw = 1, color = :grey, ls = :dash, label = "Initial shape"); 
    plot!(-X0, Y0, lw = 1, color = :grey, ls = :dash, label = "");

    # next, color by morphogen concentration 
    for i in 1:Ndisc-1
        plot!(X[i:i+1], Y[i:i+1], lw = 3, color = get(cmap, nphi[i]), label = "");
        plot!(-X[i:i+1], Y[i:i+1], lw = 3, color = get(cmap, nphi[i]), label = "");
    end

    # other plot necessities
    xlabel!("x"); ylabel!("y");
    if titleTxt
        title!(titleTxt)
    end

    return plt 
end

plt = visualise(x0quad.(Si), y0quad.(Si), Phi0 .* Si);
display(plt)






