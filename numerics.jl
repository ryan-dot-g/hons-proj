### TODO 

using LinearAlgebra, DifferentialEquations, JuMP, Ipopt, HiGHS, Plots, LaTeXStrings; 
println("running...");

############ -------------------------------------------------- ############
############ ---------------- FILE HANDLING ------------------- ############
############ -------------------------------------------------- ############
gitFP = "C:/Users/rgray/OneDrive/ryan/Uni/HONOURS/hons-proj";
suppressFP = gitFP * "/suppress.txt";

############ -------------------------------------------------- ############
############ ------------- PHYSICAL PARAMETERS ---------------- ############
############ -------------------------------------------------- ############
R = 130; # radius of force-free reference sphere (130)
r = 150; # TODO radius of initial condition sphere (150 for now) 
E0 = 440; # base Young's modulus (440)
h = 20 # bilayer thickness (20)
b = 5; # exponential decrease of stiffness (5)
κ = 1.0; # TODO morphogen stress upregulation coefficient (1 for now)
ζ = 0.2; # morphogen decay rate (0.2)
D = 5.0; # morphogen diffusion coefficient (0.01)

E(ϕ) = E0 * h * exp(-b*ϕ); # elasticity function
f(ϵ) = κ * ϵ; # strain-dependent morphogen expression function 

############ -------------------------------------------------- ############
############ ------------- NUMERICAL PARAMETERS --------------- ############
############ -------------------------------------------------- ############
Ndisc = 50; # number of discretisation points on s (50)
smin = -π/2; smax = π/2; # bounds of s values (-π/2, π/2)
dt = 0.3; # time discretisation (0.1)
tmax = 10*dt; # max time (1e3 * dt)
Ω = 1e4; # large number: punishing potential for volume deviation (1e8)
ω = 0.0; # small number: surface friction 
dsint = 0.0; # small number: distance inside the s grid to start at to avoid div0

############ -------------------------------------------------- ############
############ ----------- VISUALISATION PARAMETERS ------------- ############
############ -------------------------------------------------- ############
plotRes = Inf; # how many timesteps per plot
cmap = cgrad(:viridis); # colormap for morphogen concentration 

############ -------------------------------------------------- ############
############ ---------------- NUMERICAL SETUP ----------------- ############
############ -------------------------------------------------- ############
Si = range(smin+dsint, smax-dsint, Ndisc); # grid of s values 
ds = Si[2] - Si[1]; 
Tn = 0:dt:tmax; # grid of t values 
Ntimes = length(Tn); 
interior = 2:Ndisc-1; # index set for the interior of the array

# preparing quantities to be saved. Each qty is saved BEFORE and AFTER sim
ϕtot = zeros(Ntimes); 
ϵtot = zeros(Ntimes);

############ -------------------------------------------------- ############
############ ------------- DISCRETISED OPERATORS -------------- ############
############ -------------------------------------------------- ############

# first deriv 
P = Tridiagonal(fill(-1.0, Ndisc-1), fill(0.0, Ndisc), fill(1.0, Ndisc-1)); # main section 
P = Matrix(P); # convert out of sparse to modify elements 
P[1, 1] = -3; P[1, 2] = 4; P[1, 3] = -1; # forward diff top row 
P[Ndisc, Ndisc] = 3; P[Ndisc, Ndisc-1] = -4; P[Ndisc, Ndisc-2] = 1; # backward diff bottom row 
P /= (2*ds); # scale 

# integral over total s domain, using trapezoid rule
integrateSdom(Qty) = ds * ( sum(Qty) - (Qty[1] + Qty[end])/2 );

############ -------------------------------------------------- ############
############ -------- STEADY-STATE / INITIAL FUNCTIONS -------- ############
############ -------------------------------------------------- ############
ϵ0sc = 0.5 * ((r^2 - R^2)/R^2)^2; # scalar initial (constant) trace of strain tensor squared
ϵ0 = zeros(Ndisc) .+ ϵ0sc; # over s-grid
V0 = 4/3 * r^3; # initial volume, normalised by pi

ϕ0sc = 1/ζ * f.(ϵ0sc); # steady-state morphogen scalar
ϕ0 = zeros(Ndisc) .+ ϕ0sc; # steady-state morphogen over the full grid 
ϕ0test = zeros(Ndisc) .+ SinS * 3 .+ 3; # non-uniform but obeys BCs

X0 = r * cos.(Si); Y0 = r * sin.(Si); # initial condition shape 
X0dash = -Y0[:]; Y0dash = X0[:]; # initial derivatives 
Xundef = R * cos.(Si); Yundef = R * sin.(Si); # undeformed shape

# test approximations of r cos(s), r sin(s) to use as initial guesses for x(s), y(s) to test
# shape evolution. 
# To ensure BCs, using quadratic approximation of x(s) with roots at +-pi/2 plus an interesting perturbation 
# Just using the actual y(s) function 
x0temp(s) = r * ( 1 - 4/π^2 * s^2 + 0.2 * (s^2 - π^2/4) * s );
y0test(s) = r * sin(s);
x0tempDash(s) = r * ( -8/π^2 * s + 0.2 * (3*s^2 - π^2/4) );
y0testDash(s) = r * cos(s); 
# normalise so the volume is satisfied
vtemp = VolInt(x0temp.(Si), y0test.(Si), x0tempDash.(Si), y0testDash.(Si));
# x^2 term in integral so normalise by square root of desired volume
# and add slight noise so punishing potential is similar order for this and X0 
normK = sqrt(V0/vtemp) * 1.000001; 
x0test(s) = normK * x0temp(s); x0testDash(s) = normK * x0tempDash(s);
# finally create them as arrays 
X0test = x0test.(Si); Y0test = y0test.(Si); 
X0testDash = x0testDash.(Si); Y0testDash = y0testDash.(Si);

############ -------------------------------------------------- ############
############ ------------- PRESAVING CALCULATIONS ------------- ############
############ -------------------------------------------------- ############
CosS = cos.(Si); SinS = sin.(Si); # presaving trigs 
volEl = R^2 * CosS; # volume element

############ -------------------------------------------------- ############
############ ------------ FUNCTIONALS AND QTYs --------------- ############
############ -------------------------------------------------- ############
function trStrSq(X, Y, Xdash, Ydash)
    # takes vectors of the surface parameterisation x(s), y(s) and derivatives w.r.t. s
    # returns a vector containing the trace of the square of the strain tensor at each s
    t1 = (Xdash).^2 .+ (Ydash).^2 .- R^2;
    t2 = t1 * 0; 
    # handle boundary different if s coords go all the way to the boundary 
    if Si[1] == -π/2
        t2[interior] = X[interior].^2 ./ (CosS[interior].^2); # interior is just x/cos^2 
        t2[[1,end]] .= r^2; # boundary is when x(s)~r cos(s)
    else
        t2 = X.^2 ./ CosS.^2;
    end
    t2 = t2 .- R^2; 

    return 1/(4*R^4) * (t1.^2 .+ t2.^2); 
end

function ElEn(ϕ, X, Y, Xdash, Ydash)
    # Elastic energy functional 
    fn = E.(ϕ)/2 .* trStrSq(X, Y, Xdash, Ydash) .* volEl; # integrand over S0 including volel
    return integrateSdom(fn); 
end

function VolInt(X, Y, Xdash, Ydash)
    # volume integral 
    fn = X.^2 .* Ydash; # integrand normalised by pi 
    return integrateSdom(fn);
end

function SurfaceFriction(X, Y)
    # number that models friction, as the difference between the COM and zero
    # center of mass in x direction is naturally 0 by symm
    ycom = integrateSdom(Y);
    return (ycom)^2
end

function EnergyFunctional(ϕ, X, Y, Xdash, Ydash, Ω)
    # total energy functional, including elastic energy and a punishing potential 
    # takes punishing potential as input so it can be run with/without volume constraint at various strengths
    # additionally takes surface friction 
    return ElEn(ϕ, X, Y, Xdash, Ydash) + Ω * (VolInt(X, Y, Xdash, Ydash) - V0)^2 +
                ω * SurfaceFriction(X, Y);
end

function MorphogenFunctional(ϕ, ϕn, X, Y, ϵ)
    # total energy for the morphogen equation 
    integrand1 = (ϕ .- ϕn).^2 /(2*dt) .* volEl; # time-deriv bit
    integrand2 = -f.(ϵ) .* ϕ .* volEl; # morphogen production bit
    integrand3 = ζ * (ϕ).^2 /2 .* volEl; # disassociation bit
    integrand4 = D * (P*ϕ).^2 /2 .* CosS; # diffusion bit, volEl already incorporated
    
    integrandTot = integrand1 .+ integrand2 .+ integrand3 .+ integrand4;
    return integrateSdom(integrandTot);
end

############ -------------------------------------------------- ############
############ ------------------ UPDATE STEPS ------------------ ############
############ -------------------------------------------------- ############
function Renormalise!(X, Y, Xdash, Ydash, V0)
    # takes a shape parameterisation, and modifies in-place to ensure that the volume 
    # remains equal to some V0. 
    # uses the fact that integral propto x^2 y', so cube root normalisation constant 
    # to apply equally to x, y, x', y' and ensure volume works
    vtemp = VolInt(X, Y, Xdash, Ydash);
    normk = cbrt(V0/vtemp);
    X .= X * normk; Y .= Y * normk; 
    Xdash .= Xdash * normk; Ydash .= Ydash * normk;
end 

function UpdateShape!(ϕ, X, Y, Xdash, Ydash)
    # takes the current state data (shape and morphogen concentration), and updates X, Y 
    # and derivatives IN PLACE, by minimising the energy functional 
    # shapeModel = Model(Ipopt.Optimizer)
    # @variable(shapeModel, xModel[1:Ndisc]); @variable(shapeModel, yModel[1:Ndisc]);
    # set_start_value.(xModel, X); set_start_value.(yModel, Y);

    # @constraint(shapeModel, xModel[1] == 0); @constraint(shapeModel, xModel[end] == 0);
    # @constraint(shapeModel, (P*yModel)[1] == 0); @constraint(shapeModel, (P*yModel)[end] == 0);

    # @objective(shapeModel, Min, EnergyFunctional(ϕ, xModel, yModel));
    # redirect_stdout((()->optimize!(shapeModel)),open(suppressFP, "w"));
    # X .= JuMP.value.(xModel); Y .= JuMP.value.(yModel);

    # finally, update derivatives 
    Xdash .= P*X; Ydash .= P*Y;
end;

function UpdateShapeTest!(ϕ, X, Y, Xdash, Ydash, Ω)
    # takes the current state data (shape and morphogen concentration), and updates X, Y 
    # and derivatives IN PLACE, by minimising the energy functional 
    shapeModel = Model(Ipopt.Optimizer)
    @variable(shapeModel, xModel[1:Ndisc]); @variable(shapeModel, yModel[1:Ndisc]);
    set_start_value.(xModel, X); set_start_value.(yModel, Y);

    @constraint(shapeModel, xModel[1] == 0); @constraint(shapeModel, xModel[end] == 0);
    @constraint(shapeModel, (P*yModel)[1] == 0); @constraint(shapeModel, (P*yModel)[end] == 0);
    # @constraint(shapeModel, VolInt(xModel, yModel, P*xModel, P*yModel) == V0);

    @objective(shapeModel, Min, EnergyFunctional(ϕ, xModel, yModel, P*xModel, P*yModel, Ω));
    redirect_stdout((()->optimize!(shapeModel)),open(suppressFP, "w"));
    X .= JuMP.value.(xModel); Y .= JuMP.value.(yModel);

    # finally, update derivatives 
    Xdash .= P*X; Ydash .= P*Y;
end;

function UpdateMorphogen!(ϕ, X, Y, ϵ)
    # takes the current state data, and updates the morphogen concentration 
    # IN PLACE, by performing an implicit update step 
    ϕn = ϕ;
    morphogenModel = Model(Ipopt.Optimizer);
    @variable(morphogenModel, ϕModel[1:Ndisc]);
    set_start_value.(ϕModel, ϕn);

    @constraint(morphogenModel, (P*ϕModel)[1] == 0);
    @constraint(morphogenModel, (P*ϕModel)[end] == 0);

    @objective(morphogenModel, Min, MorphogenFunctional(ϕModel, ϕn, X, Y, ϵ));
    redirect_stdout((()->optimize!(morphogenModel)),open(suppressFP, "w"));

    ϕ .= JuMP.value.(ϕModel);
end;

function SaveData!(n, ϕ, X, Y, Xdash, Ydash, ϵ,
                    ϕtot, ϵtot)
    # Saves data relevant to the current state IN PLACE. Takes state data,
    #   and relevant data vecs to save into
    # ϕtot: total morphogen present in the organism
    # ϵtot: total strain
    ϕtot[n+1] = integrateSdom(ϕ); 
    ϵtot[n+1] = integrateSdom(ϵ);
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

    # dummy bit to get a good color scale 
    if ϕ0[1] != Inf
        plot!([0.0, 0.01], [0.0, 0.01], lw = 0, 
                line_z = [0.0, 2*ϕ0[1]], c = cmap, label = "");
    end

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
# initialise X, Y, ϕ arrays from the various possibilities 
# shape arrays first
X = zeros(Ndisc); Y = zeros(Ndisc); Xdash = zeros(Ndisc); Ydash = zeros(Ndisc);
X .= X0; Y .= Y0; 
Xdash .= X0dash; Ydash .= Y0dash;
ϵ = trStrSq(X, Y, Xdash, Ydash)

# now morphogen array
ϕ = zeros(Ndisc);
ϕ .= ϕ0test;

runsim = true 
if runsim
    global ϵ; 
    SaveData!(0, ϕ, X, Y, Xdash, Ydash, ϵ,
                    ϕtot, ϵtot)

    for n = 1:Ntimes-1
        t = Tn[n];

        # Step A1: update shape of X(s), Y(s) and derivatives 
        UpdateShape!(ϕ, X, Y, Xdash, Ydash); 

        # Step A2: compute new strain tensor squared 
        ϵ = trStrSq(X, Y, Xdash, Ydash); 

        # Step A3: evolve morphogen
        UpdateMorphogen!(ϕ, X, Y, ϵ);

        # Step C1: save necessary data 
        SaveData!(n, ϕ, X, Y, Xdash, Ydash, ϵ,
                    ϕtot, ϵtot)

        # Step C2: visualise if necessary, always visualising first and last 
        if mod(n-1, plotRes) == 0 || n == Ntimes-1
            tstr = round(t, digits = 1); 
            pltAnim = visualise(ϕ, X, Y, "t = $tstr");
            display(pltAnim);
        end
    end

    plt = plot(Tn, ϕtot, title = "Morphogen vs time", xlabel = "t", ylabel = "Total morphogen",
                label = L"\varphi");
    display(plt); 
end


plt = visualise(ϕ0, X0test, Y0test, "before"); display(plt);
UpdateShapeTest!(ϕ0, X0test, Y0test, X0testDash, Y0testDash, 0)
plt = visualise(ϕ0, X0test, Y0test, "middle"); display(plt);
UpdateShapeTest!(ϕ0, X0test, Y0test, X0testDash, Y0testDash, Ω)
plt = visualise(ϕ0, X0test, Y0test, "after"); display(plt);

# try again
X0test = x0test.(Si); Y0test = y0test.(Si); X0testDash = x0testDash.(Si); Y0testDash = y0testDash.(Si);
plt = visualise(ϕ0, X0test, Y0test, "before 2"); display(plt);
UpdateShapeTest!(ϕ0, X0test, Y0test, X0testDash, Y0testDash, 0)
plt = visualise(ϕ0, X0test, Y0test, "middle 2"); display(plt);

vtemp = VolInt(X0test, Y0test, X0testDash, Y0testDash);
println(vtemp)
normk = cbrt(V0/vtemp);
X0test *= normk; Y0test *= normk;
vtemp2 = VolInt(X0test, Y0test, X0testDash, Y0testDash);
println(vtemp2)
plt = visualise(ϕ0, X0test, Y0test, "after 2"); display(plt);

# 3rd try
X0test .= X0; Y0test .= Y0; X0testDash .= X0dash; Y0testDash .= Y0dash;
plt = visualise(ϕ0test, X0test, Y0test, "before 3"); display(plt);
UpdateShapeTest!(ϕ0test, X0test, Y0test, X0testDash, Y0testDash, 0)
plt = visualise(ϕ0test, X0test, Y0test, "middle 3"); display(plt);

vtemp = VolInt(X0test, Y0test, X0testDash, Y0testDash);
println(vtemp)
normk = cbrt(V0/vtemp);
X0test *= normk; Y0test *= normk; X0testDash *= normk; Y0testDash *= normk;
vtemp2 = VolInt(X0test, Y0test, X0testDash, Y0testDash);
println(vtemp2)
plt = visualise(ϕ0test, X0test, Y0test, "after 3"); display(plt);

