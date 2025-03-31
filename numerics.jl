### TODO 
# convert operators to sparse?
# fix color mappings so they r fixed despite changing phimin, phimax, and add colorbar 

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
Ndisc = 1000; # number of discretisation points on s 
smin = -pi/2; smax = pi/2; # bounds of s values
dt = 0.1; # time discretisation
tmax = 1; # max time 


############ -------------------------------------------------- ############
############ ----------- VISUALISATION PARAMETERS ------------- ############
############ -------------------------------------------------- ############
plotRes = 5; # how many timesteps per plot
cmap = cgrad(:viridis); # colormap for morphogen concentration 

############ -------------------------------------------------- ############
############ ---------------- NUMERICAL SETUP ----------------- ############
############ -------------------------------------------------- ############
Si = range(smin, smax, Ndisc); # grid of s values 
ds = Si[2] - Si[1]; 
Tn = 0:dt:tmax; # grid of t values 
Ntimes = length(Tn); 

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
############ ------------- FUNCTIONALS AND QTYs --------------- ############
############ -------------------------------------------------- ############
function trStrSq(X, Y, Xdash, Ydash)
    # takes vectors of the surface parameterisation x(s), y(s) and derivatives w.r.t. s
    # returns a vector containing the trace of the square of the strain tensor at each s
    t1 = (P*Xdash).^2 .+ (P*Ydash).^2 .- R^2;
    t2 = t1*0; 
    t2[2:end-1] = X[2:end-1].^2 ./ (CosS[2:end-1].^2); # interior is just x/cos^2 
    t2[1] = r^2; t2[end] = r^2; # boundary is when x(s)~r cos(s)
    t2 = t2 .- R^2; 

    en = 1/(4*R^4) * (t1.^2 .+ t2.^2);
    return en; 
end


############ -------------------------------------------------- ############
############ ---------- AUXILIARY / USEFUL FUNCTIONS ---------- ############
############ -------------------------------------------------- ############
ϕ0sc = 1/ζ * f(0.5 * ((r^2 - R^2)/R^2)^2); # steady-state morphogen scalar
ϕ0 = Si .* 0 .+ ϕ0sc; # steady-state morphogen over the full grid 
trStrSq0 = 0.5 * ((r^2 - R^2)/R^2)^2; # initial (constant) trace of strain tensor squared

X0 = r * cos.(Si); Y0 = r * sin.(Si); # initial condition shape 
X0dash = -Y0; Y0dash = X0; # initial derivatives 
Xundef = R * cos.(Si); Yundef = R * sin.(Si); # undeformed shape

# 2nd-order approximations of r cos(s), r sin(s) to use as initial guesses for 
#   x(s), y(s) at steady-state
x0quad(s) = r * ( 1 - 4/pi^2 * s^2 + 0.2 * (s^2 - pi^2/4) * s );
y0quad(s) = r * (s);

############ -------------------------------------------------- ############
############ ------------- PRESAVING CALCULATIONS ------------- ############
############ -------------------------------------------------- ############
CosS = cos.(Si); SinS = sin.(Si); # presaving trigs 
TanS = 0 * CosS; TanS[2:end-1] = tan.(Si[2:end-1]); # presaving Tan, with 0 on bdy to avoid undefineds

############ -------------------------------------------------- ############
############ ------------------ UPDATE STEPS ------------------ ############
############ -------------------------------------------------- ############
function UpdateShape!(ϕ, X, Y, Xdash, Ydash)
    # takes the current state data (shape and morphogen concentration), and updates X, Y 
    # and derivatives IN PLACE, by minimising the energy functional 

    # finally, update derivatives 
    Xdash = P*X; Ydash = P*Y;
end;

function UpdateMorphogen!(ϕ, X, Y, ϵ)
    # takes the current state data, and updates the morphogen concentration 
    # IN PLACE, by performing an implicit update step 
    ϕ = ϕ * 0 .+ ϕ;
end;

function ApplyBCs!(ϕ, X, Y)
    # takes the current state, and enforces boundary conditions 
    # Dirichelet BCs in X, and Neumann in ϕ and Y 
    X = X*0 .+ X; 
end;

############ -------------------------------------------------- ############
############ ------------ VISUALISATION FUNCTIONS ------------- ############
############ -------------------------------------------------- ############
function visualise(ϕ, X, Y, titleTxt = false)
    # takes a surface discretised parameterisation and morphogen concentration
    # displays a plot of the deformed hydra, colored by the morphogen concentration 
    # optionally adds a provided title 
    # returns the plot object for displaying and saving

    nphi = (ϕ .- minimum(ϕ)) ./ (maximum(ϕ) - minimum(ϕ)); # normalise for color 

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
    lm = 2 * r; xlims!(-lm, lm); ylims!(-lm, lm);
    if titleTxt isa String
        title!(titleTxt)
    end

    return plt;
end

plt = visualise(ϕ0 .* Si, x0quad.(Si), y0quad.(Si));
display(plt)


############ -------------------------------------------------- ############
############ --------------- CHECKING FORMULAS ---------------- ############
############ -------------------------------------------------- ############
checkFormulas = false
if checkFormulas 
    X = r * CosS; Y = r * SinS; 
    Xdash = P*X; Ydash = P*Y; 

    # check trace of strain tensor calc (CHECKED)
    en = trStrSq(X, Y, Xdash, Ydash); 
    # println(maximum( abs.(en .- trStrSq0) )); 

    # check first-derivative matrix (CHECKED)
    # error on bdy is about double error on interior, but still scales with ds as expected
    XdashTrue = - r * SinS;
    # println(XdashTrue .- Xdash);
    # println(maximum( abs.(XdashTrue .- Xdash) ));

    # check second-derivative matrix (CHECKED)
    # error on bdy is much (10^6) higher than error on interior but does scale well 
    Xddash = Q*X;
    XddashTrue = -X;
    println(XddashTrue .- Xddash);
    println(maximum( abs.(XddashTrue .- Xddash) ));
end



############ -------------------------------------------------- ############
############ ---------------- TIME EVOLUTION ------------------ ############
############ -------------------------------------------------- ############
X = X0; Y = Y0; 
Xdash = X0dash; Ydash = Y0dash;
ϕ = ϕ0 .* Si; 
runsim = false 
if runsim
    for n = 1:Ntimes 
        t = Tn[n];

        # Step 1: update shape of X(s), Y(s) and derivatives 
        UpdateShape!(ϕ, X, Y, Xdash, Ydash); 

        # Step 2: compute new strain tensor squared 
        ϵ = trStrSq(X, Y, Xdash, Ydash); 

        # Step 3: evolve morphogen
        UpdateMorphogen!(ϕ, X, Y, ϵ); 

        # Step 4: apply BCs 
        ApplyBCs!(ϕ, X, Y); 

        # Step 5: visualise if necessary 
        if mod(n, plotRes) == 0
            tstr = round(t); 
            pltAnim = visualise(ϕ, X, Y, "t = $tstr");
            display(pltAnim);
        end
    end
end

