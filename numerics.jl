### TODO 
# more sophisticated numerical integral?

using LinearAlgebra, DifferentialEquations, JuMP, Ipopt, HiGHS
        Plots, LaTeXStrings; 
println("running...");

############ -------------------------------------------------- ############
############ ------------- PHYSICAL PARAMETERS ---------------- ############
############ -------------------------------------------------- ############
R = 130; # radius of force-free reference sphere (130)
r = 150; # TODO radius of initial condition sphere 
E0 = 440; # base Young's modulus (440)
h = 20 # bilayer thickness (20)
b = 5; # exponential decrease of stiffness (5)
κ = 1; # TODO morphogen stress upregulation coefficient 
ζ = 0.2; # morphogen decay rate (0.2)
D = 0.0; # 0.01; # morphogen diffusion coefficient (0.01)

E(ϕ) = E0 * h * exp(-b*ϕ); # elasticity function
f(ϵ) = κ * ϵ; # strain-dependent morphogen expression function 

############ -------------------------------------------------- ############
############ ------------- NUMERICAL PARAMETERS --------------- ############
############ -------------------------------------------------- ############
Ndisc = 100; # number of discretisation points on s 
smin = -π/2; smax = π/2; # bounds of s values
dt = 0.1; # time discretisation
tmax = 50.0; # max time 
Ω = 1e8; # punishing potential for volume deviation 

############ -------------------------------------------------- ############
############ ----------- VISUALISATION PARAMETERS ------------- ############
############ -------------------------------------------------- ############
plotRes = Inf; # how many timesteps per plot
cmap = cgrad(:viridis); # colormap for morphogen concentration 

############ -------------------------------------------------- ############
############ ---------------- NUMERICAL SETUP ----------------- ############
############ -------------------------------------------------- ############
Si = range(smin, smax, Ndisc); # grid of s values 
ds = Si[2] - Si[1]; 
Tn = 0:dt:tmax; # grid of t values 
Ntimes = length(Tn); 

# preparing quantities to be saved. Each qty is saved BEFORE and AFTER sim
ϕtot = zeros(Ntimes); 

############ -------------------------------------------------- ############
############ ------------- DISCRETISED OPERATORS -------------- ############
############ -------------------------------------------------- ############
# identity 
In = I(Ndisc) * 1.0;

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
############ ---------- AUXILIARY / USEFUL FUNCTIONS ---------- ############
############ -------------------------------------------------- ############
trStrSq0 = 0.5 * ((r^2 - R^2)/R^2)^2; # initial (constant) trace of strain tensor squared
ϕ0sc = 1/ζ * f(trStrSq0); # steady-state morphogen scalar
ϕ0 = zeros(Ndisc) .+ ϕ0sc; # steady-state morphogen over the full grid 
vol0 = 4/3 * π * r^3;

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
V = ((1/dt) + ζ)*In .+ D/R^2*(diagm(TanS)*P .- Q); # morphogen update matrix 
Vinv = inv(V); # inverse save speeds up calc 

############ -------------------------------------------------- ############
############ ------------ FUNCTIONALS AND QTYs --------------- ############
############ -------------------------------------------------- ############
function trStrSq(X, Y, Xdash, Ydash)
    # takes vectors of the surface parameterisation x(s), y(s) and derivatives w.r.t. s
    # returns a vector containing the trace of the square of the strain tensor at each s
    t1 = (P*Xdash).^2 .+ (P*Ydash).^2 .- R^2;
    t2 = t1 * 0; 
    t2[2:end-1] = X[2:end-1].^2 ./ (CosS[2:end-1].^2); # interior is just x/cos^2 
    t2[1] = r^2; t2[end] = r^2; # boundary is when x(s)~r cos(s)
    t2 = t2 .- R^2; 

    en = 1/(4*R^4) * (t1.^2 .+ t2.^2);
    return en; 
end

function ElEn(ϕ, X, Y, Xdash, Ydash)
    # Elastic energy functional 
    fn = E.(ϕ)/2 .* trStrSq(X, Y, Xdash, Ydash); # integrand
    return sum(fn) * ds; 
end

function VolInt(X, Y, Xdash, Ydash)
    # volume integral 
    fn = X.^2 .* Ydash; # integrand 
    return sum(fn) * ds;
end

function EnergyFunctional(ϕ, X, Y)
    # total energy functional, including elastic energy and a punishing potential 
    Xdash = P*X; Ydash = P*Y;
    return ElEn(ϕ, X, Y, Xdash, Ydash) + Ω * (VolInt(X, Y, Xdash, Ydash) - vol0)^2;
end

############ -------------------------------------------------- ############
############ ------------------ UPDATE STEPS ------------------ ############
############ -------------------------------------------------- ############
function UpdateShape!(ϕ, X, Y, Xdash, Ydash)
    # takes the current state data (shape and morphogen concentration), and updates X, Y 
    # and derivatives IN PLACE, by minimising the energy functional 

    # finally, update derivatives 
    Xdash .= P*X; Ydash .= P*Y;
end;

function UpdateMorphogen!(ϕ, X, Y, ϵ)
    # takes the current state data, and updates the morphogen concentration 
    # IN PLACE, by performing an implicit update step 
    ϕ .= Vinv * ( (1/dt)*ϕ .+ f.(ϵ));
end;

function ApplyBCs!(ϕ, X, Y)
    # takes the current state, and enforces boundary conditions 
    # Dirichelet BCs in X, and Neumann in ϕ and Y 
    X[1] = 0; X[end] = 0; # Dirichelet BCs in X 

    Y[1] = 1/3 * (4*Y[2] - Y[3]);
    Y[end] = 1/3 * (4*Y[end-1] - Y[end-2]);

    ϕ[1] = 1/3 * (4*ϕ[2] - ϕ[3]);
    ϕ[end] = 1/3 * (4*ϕ[end-1] - ϕ[end-2]);
end;

function SaveData!(n, ϕ, X, Y, Xdash, Ydash, ϵ,
                    ϕtot)
    # Saves data relevant to the current state IN PLACE. Takes state data,
    #   and relevant data vecs to save into
    # ϕtot: total morphogen present in the organism, appropriately normalised 

    ϕtot[n+1] = sum(ϕ) * ds; 
end

############ -------------------------------------------------- ############
############ ------------ VISUALISATION FUNCTIONS ------------- ############
############ -------------------------------------------------- ############
function visualise(ϕ, X, Y, titleTxt = false)
    # takes a surface discretised parameterisation and morphogen concentration
    # displays a plot of the deformed hydra, colored by the morphogen concentration 
    # optionally adds a provided title 
    # returns the plot object for displaying and saving

    # Plot actual deformed shape, colored by morphogen concentration 
    plt = plot(X, Y, line_z = ϕ, lw = 3, c = cmap, label = "Hydra shape");
    plot!(-X, Y, line_z = ϕ, lw = 3, c = cmap, label = "")

    # plot initial shape 
    plot!(X0, Y0, lw = 1, color = :grey, ls = :dash, label = "Initial shape"); 
    plot!(-X0, Y0, lw = 1, color = :grey, ls = :dash, label = "");

    # other plot necessities
    xlabel!("x"); ylabel!("y");
    lm = 2 * r; xlims!(-lm, lm); ylims!(-lm, lm);
    if titleTxt isa String
        title!(titleTxt)
    end

    return plt;
end

# plt = visualise(ϕ0 .* Si, x0quad.(Si), y0quad.(Si));
# display(plt)


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
ϵ = trStrSq(X, Y, Xdash, Ydash); 
ϕ = ϕ0; ϕ = zeros(Ndisc) .+ Si .+ smax;
runsim = false 
if runsim
    global ϵ
    SaveData!(0, ϕ, X, Y, Xdash, Ydash, ϵ,
                    ϕtot)

    for n = 1:Ntimes-1
        t = Tn[n];

        # Step A1: update shape of X(s), Y(s) and derivatives 
        UpdateShape!(ϕ, X, Y, Xdash, Ydash); 

        # Step A2: compute new strain tensor squared 
        ϵ = trStrSq(X, Y, Xdash, Ydash); 

        # Step A3: evolve morphogen
        UpdateMorphogen!(ϕ, X, Y, ϵ); 

        # Step B1: apply BCs 
        ApplyBCs!(ϕ, X, Y); 

        # Step C1: save necessary data 
        SaveData!(n, ϕ, X, Y, Xdash, Ydash, ϵ,
                    ϕtot)

        # Step C2: visualise if necessary, always visualising first and last 
        if mod(n-1, plotRes) == 0 || n == Ntimes-1
            tstr = round(t, digits = 1); 
            pltAnim = visualise(ϕ, X, Y, "t = $tstr");
            display(pltAnim);
        end
    end

    plt = plot(Tn, ϕtot);
    display(plt); 
end

ff(x) = sum(x.^2);

shapeModel = Model(Ipopt.Optimizer)
@variable(shapeModel, xModel[1:Ndisc]); @variable(shapeModel, yModel[1:Ndisc]);
@constraint(shapeModel, xModel[1] == 0); @constraint(shapeModel, xModel[end] == 0);
@constraint(shapeModel, (P*yModel)[1] == 0); @constraint(shapeModel, (P*yModel)[end] == 0);
@objective(shapeModel, Min, EnergyFunctional(ϕ, xModel, yModel));
optimize!(shapeModel);
x = JuMP.value.(xModel); y = JuMP.value.(yModel);

println("Done");
