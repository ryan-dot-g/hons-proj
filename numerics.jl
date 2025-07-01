### TODO 
# SHORT TERM:
# - number of steps until failure, not total number for time thing
# LONG TERM:
# - Get surface friction working better: could pin the base (current fix), or penalise movement
# - Relax convergence requirements to speed shape up (tried loosely but could try harder)

using LinearAlgebra, DifferentialEquations, Plots, LaTeXStrings;
using Optim, ForwardDiff;
println("running...");

############ -------------------------------------------------- ############
############ ---------------- FILE HANDLING ------------------- ############
############ -------------------------------------------------- ############
gitFP = "C:/Users/rgray/OneDrive/ryan/Uni/HONOURS/hons-proj";
simFP = "C:/Users/rgray/OneDrive/ryan/Uni/HONOURS/hons-proj/vids"

############ -------------------------------------------------- ############
############ ------------- PHYSICAL PARAMETERS ---------------- ############
############ -------------------------------------------------- ############
# Length units in μm, time units in hr 
R0 = 130; # radius of force-free reference sphere (130)
r = 160; # radius of initial condition sphere (160 for now) 
E0 = 440.0; # base Young's modulus (440). NOTE: this interacts with Ω
h = 20.0 # bilayer thickness (20). NOTE: this interacts with Ω
b = 5.0; # exponential decrease of stiffness (5). NOTE: this interacts with Ω
κ = 1.0; # TODO morphogen stress upregulation coefficient (1 for now)
ζ = 0.2; # morphogen decay rate (0.2)
D = 0.01; # morphogen diffusion coefficient (0.01)

# PARAM SETS OF NOTE: 
# - (A) 

E(ϕ) = E0 * h * exp(-b*ϕ); # elasticity function
f(ϵ) = κ * ϵ; # strain-dependent morphogen expression function 

############ -------------------------------------------------- ############
############ ------------- NUMERICAL PARAMETERS --------------- ############
############ -------------------------------------------------- ############
Ndisc = 40; # number of discretisation points on s (40)
smin = -π/2; smax = π/2; # bounds of s values (-π/2, π/2)
dt = 0.03; # time discretisation (0.04)
tmax = 1e2*dt; # max time (5e2 * dt for evec). Sims dont really get here tho
Ω = 1e2; # not too large number: punishing potential for volume deviation (1e2)
ω = 1e1; # not too large number: surface friction (1e1)
dsint = 0.01; # small number: distance inside the s grid to start at to avoid div0
# dsint = (π/Ndisc)/2;

cfl = 2 * D * dt / (π/Ndisc)^2; # cfl condition for diffusion, checking sensible. ds is approximate 
# println("CFL number (should be <<1): $(round(cfl, digits = 5))")

############ -------------------------------------------------- ############
############ ----------- VISUALISATION PARAMETERS ------------- ############
############ -------------------------------------------------- ############
plotRes = 1; # how many timesteps per plot (Inf for just 1st and last plots)
cmap = cgrad(:viridis); # colormap for morphogen concentration 

############ -------------------------------------------------- ############
############ ---------------- NUMERICAL SETUP ----------------- ############
############ -------------------------------------------------- ############
Si = range(smin+dsint, smax-dsint, Ndisc); # grid of s values 
ds = Si[2] - Si[1]; 
Tn = 0:dt:tmax; # grid of t values 
Ntimes = length(Tn); 
Ncutoff = Ntimes; # where the sim ends because it breaks
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

P[1, 1] = -3; P[1, 2] = 4; P[1, 3] = -1; # forward diff top row 3rd order
P[Ndisc, Ndisc] = 3; P[Ndisc, Ndisc-1] = -4; P[Ndisc, Ndisc-2] = 1; # backward diff bottom row 

# P[1,1:2] = 2*[-1, 1] # testing first-order difference 
# P[Ndisc, Ndisc-1:Ndisc] = 2*[-1, 1]

P /= (2*ds); # scale 

# integral over total s domain, using trapezoid rule
integrateSdom(Qty) = ds * ( sum(Qty) - (Qty[1] + Qty[end])/2 );

############ -------------------------------------------------- ############
############ ------------- PRESAVING CALCULATIONS ------------- ############
############ -------------------------------------------------- ############
CosS = cos.(Si); SinS = sin.(Si); # presaving trigs 
volEl = R0^2 * CosS; # volume element

############ -------------------------------------------------- ############
############ -------- STEADY-STATE / INITIAL FUNCTIONS -------- ############
############ -------------------------------------------------- ############
ϵ0sc = 0.5 * ((r^2 - R0^2)/R0^2)^2; # scalar initial (constant) trace of strain tensor squared
ϵ0 = zeros(Ndisc) .+ ϵ0sc; # over s-grid
V0 = 4/3 * r^3; # initial volume, normalised by pi

ϕ0sc = 1/ζ * f.(ϵ0sc); # steady-state morphogen scalar
ϕ0 = zeros(Ndisc) .+ ϕ0sc; # steady-state morphogen over the full grid 

X0 = r * CosS; Y0 = r * SinS; # initial condition shape 
X0dash = -Y0[:]; Y0dash = X0[:]; # initial derivatives 
Xundef = R0 * CosS; Yundef = R0 * SinS; # undeformed shape
rXu =  X0 ./ r; rYu =  Y0 ./ r; # radial unit vectors

ZZ0 = hcat(ϕ0, X0, Y0); # full state vector, as Ndisc * 3 matrix

############ -------------------------------------------------- ############
############ ------------ FUNCTIONALS AND QTYs --------------- ############
############ -------------------------------------------------- ############
function trStrSq(X, Y, Xdash, Ydash)
    # takes vectors of the surface parameterisation x(s), y(s) and derivatives w.r.t. s
    # returns a vector containing the trace of the square of the strain tensor at each s
    t1 = (Xdash).^2 .+ (Ydash).^2 .- R0^2;
    t2 = t1 * 0; 
    # handle boundary different if s coords go all the way to the boundary 
    if Si[1] == -π/2
        t2[interior] = X[interior].^2 ./ (CosS[interior].^2); # interior is just x/cos^2 
        t2[[1,end]] .= r^2; # boundary is when x(s)~r cos(s)
    else
        t2 = X.^2 ./ CosS.^2;
    end
    t2 = t2 .- R0^2; 

    return 1/(4*R0^4) * (t1.^2 .+ t2.^2); 
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

# function VolInt(X, Y, Xdash, Ydash)
#     # volume itnegral take 2
#     # normalised by π
#     fn = X .* (X.*Ydash .- Xdash.*Y)
#     return 2/3 * integrateSdom(fn);
# end

function SurfaceFriction(X, Y)
    # number that models friction, as the difference between the COM and zero
    # center of mass in x direction is naturally 0 by symm
    # return sum(Y)^2;

    # alternatively, pin the base 
    return (Y[1] - Y0[1])^2
end

κ_bend = 1;
function Bending(X)
    ΔX = X .- X0;
    return sum( (P*ΔX).^2 );
end

function EnergyFunctional(ϕ, X, Y, Xdash, Ydash)
    # total energy functional, including elastic energy and a punishing potential and surface friction
    return ElEn(ϕ, X, Y, Xdash, Ydash) + Ω * (VolInt(X, Y, Xdash, Ydash) - V0)^2 +
                ω * SurfaceFriction(X, Y) +
                κ_bend * (Bending(X) + Bending(Y));
end

function EFwrapped(fullData, ϕ)
    # wrapper for the energy functional. Unpacks the data, computes derivatives, then
    # returns the energy functional to optimise 
    X = fullData[1:Ndisc]; Y = fullData[Ndisc+1:end];
    Xdash = P*X; Ydash = P*Y;
    return EnergyFunctional(ϕ, X, Y, Xdash, Ydash);
end

function MorphogenFunctional(ϕ, ϕn, X, Y, ϵsq)
    # total energy for the morphogen equation 
    # ϕn is the previous state 
    integrand1 = (ϕ .- ϕn).^2 /(2*dt) .* volEl; # time-deriv bit
    integrand2 = -f.(ϵsq) .* ϕ .* volEl; # morphogen production bit
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
    # uses the fact that integral propto 3 factors of x, y, x', y', so a cube 
    # root normalisation constant can apply equally to x,y,x',y' and ensure volume works
    vtemp = VolInt(X, Y, Xdash, Ydash);
    normk = cbrt(V0/vtemp);
    X .= X * normk; Y .= Y * normk; 
    Xdash .= Xdash * normk; Ydash .= Ydash * normk;
end 

function UpdateShape!(ϕ, X, Y, Xdash, Ydash)
    # takes the current state data (shape and morphogen concentration), and updates X, Y 
    # and derivatives IN PLACE, by minimising the energy functional 
    # returns results of the optim
    lbounds = [dsint * 0.9; Vector(1:Ndisc-2)*-Inf; dsint*0.9; Vector(1:Ndisc)*-Inf ];
    ubounds = [dsint * 1.1; Vector(1:Ndisc-2)*Inf; dsint*1.1; Vector(1:Ndisc)*Inf];
    result = optimize(x -> EFwrapped(x, ϕ), 
                    [X; Y], 
                    # lbounds, ubounds, Fminbox(GradientDescent()),
                    method = LBFGS(), 
                    autodiff = :forward, iterations = 15000);  # 100 000
    fullX = result.minimizer;
    X .= fullX[1:Ndisc]; Y .= fullX[Ndisc+1:end];
    # finally, update derivatives 
    Xdash .= P*X; Ydash .= P*Y;
    return result;
end; 

function UpdateMorphogen!(ϕ, X, Y, ϵsq)
    # takes the current state data, and updates the morphogen concentration 
    # IN PLACE, by performing an implicit update step 
    ϕn = ϕ; # previous value of morphogen 
    result = optimize(ϕ -> MorphogenFunctional(ϕ, ϕn, X, Y, ϵsq), 
                    ϕn, method = LBFGS(), 
                    autodiff = :forward, iterations = 1000); 
    ϕ .= result.minimizer;
    return result; 
end

function SaveData!(n, ϕ, X, Y, Xdash, Ydash, ϵsq,
                    ϕtot, ϵtot)
    # Saves data relevant to the current state IN PLACE. Takes state data,
    #   and relevant data vecs to save into
    # ϕtot: total morphogen present in the organism
    # ϵtot: total strain
    ϕtot[n+1] = integrateSdom(ϕ); 
    ϵtot[n+1] = integrateSdom(ϵsq);
end

############ -------------------------------------------------- ############
############ --------------- EIGENVECTOR STEPS ---------------- ############
############ -------------------------------------------------- ############
function ProjectEvec!(ϕ, X, Y)
    # renormalises state to allow it to slowly project onto the leading eigenvector 
    # modifies in place, replacing the existing state with one that has the perturbation normalised
    ZZ = hcat(ϕ, X, Y) # full state 
    dZZ = ZZ .- ZZ0; # change in state 
    dZZmag = dZZ[:,1].^2 .+ (dZZ[:,2]./r).^2 .+ (dZZ[:,3]./r).^2  # length Ndisc; # Δϕ^2 + (ΔX/R)^2 for each s
    globalNorm = sqrt(integrateSdom(dZZmag .* volEl)); # norm over entire domain
    # println(globalNorm)
    dZZ = dZZ / globalNorm; # normalise the perturbation 
    ZZ = ZZ0 .+ dZZ; # get the full state again 
    ϕ .= ZZ[:,1]; X .= ZZ[:,2]; Y .= ZZ[:,3]; # save new state 
end;

############ -------------------------------------------------- ############
############ ------------------ BETTER ICs -------------------- ############
############ -------------------------------------------------- ############

# test approximations of r cos(s), r sin(s) to use as initial guesses for x(s), y(s) to test
# shape evolution. 
# To ensure BCs, using quadratic approximation of x(s) with roots at +-pi/2 plus an interesting perturbation 
# Just using the actual y(s) function 
x0temp(s) = r * ( 1 - 4/π^2 * s^2 + 0.2 * (s^2 - π^2/4) * s );
y0temp(s) = r * sin(s);
x0tempDash(s) = r * ( -8/π^2 * s + 0.2 * (3*s^2 - π^2/4) );
y0tempDash(s) = r * cos(s); 

# so that the guess isn't too far off the actual soln, construct the test as a weighted sum of 
# the correct IC and the 'interesting' IC, then adjust for correct volume
α_shape = 0.2; # weight in the 'interesting' initial condition (0.2)
X0test = α_shape*x0temp.(Si) + (1-α_shape)*X0; Y0test = α_shape*y0temp.(Si) + (1-α_shape)*Y0;
X0testDash = α_shape*x0tempDash.(Si) + (1-α_shape)*X0dash; Y0testDash = α_shape*y0tempDash.(Si) + (1-α_shape)*Y0dash;
Renormalise!(X0test, Y0test, X0testDash, Y0testDash, V0);

### --- ϕ ICs --- ###
# nonuniform but obey BCs
α_morph = 0.125; # weight of the nonuniform 'interesting' component (0.125)

# bump on the top of the domain
ϕ0tempTop = zeros(Ndisc) .+ SinS * ϕ0sc .+ ϕ0sc; 
ϕ0topBump = α_morph*ϕ0tempTop .+ (1-α_morph)*ϕ0; 

# bump on the side of the domain 
ϕ0tempSide = zeros(Ndisc) .+ 2*exp.(-2 * Si.^2) * ϕ0sc;
ϕ0sideBump = α_morph*ϕ0tempSide .+ (1-α_morph)*ϕ0; 

############ -------------------------------------------------- ############
############ ------------ VISUALISATION FUNCTIONS ------------- ############
############ -------------------------------------------------- ############
function visShape(ϕ, X, Y, titleTxt = false)
    # takes a surface discretised parameterisation and morphogen concentration
    # displays a plot of the deformed hydra, colored by the morphogen concentration 
    # optionally adds a provided title 
    # returns the plot object for displaying and saving

    # Plot actual deformed shape, colored by morphogen concentration 
    plt = plot(X, Y, line_z = ϕ, lw = 4, alpha = 0.7,
                c = cmap, label = "Hydra shape", 
                aspect_ratio = :equal, legend = :bottomright);
    plot!(-X, Y, line_z = ϕ, lw = 4, alpha = 0.7,
                c = cmap, label = "")

    # Plot material points 
    scatter!(X, Y, ms = 2.0, color = :black, label = "")
    scatter!(-X, Y, ms = 2.0, color = :black, label = "")

    # plot initial shape 
    plot!(X0, Y0, lw = 1, color = :grey, ls = :dash, label = "Initial shape"); 
    plot!(-X0, Y0, lw = 1, color = :grey, ls = :dash, label = "");

    # dummy bit to get a good color scale 
    # if ϕ0[1] != Inf
    #     plot!([0.0, 0.01], [0.0, 0.01], lw = 0, 
    #             line_z = [0.0, 2*ϕ0[1]], c = cmap, label = "");
    # end

    # other plot necessities
    xlabel!("x"); ylabel!("y");
    lm = 1.33 * r; xlims!(-lm, lm); ylims!(-lm, lm);
    if titleTxt isa String
        title!(titleTxt)
    end

    return plt;
end

function visDeform(ϕ, X, Y)
    # takes a surface discretised parameterisation and morphogen concentration
    # displays a plot of just the deformation, i.e. the difference between main shape and IC

    ΔX = X .- X0; ΔY = Y .- Y0; # compute deformation

    # Plot deformation, colored by morphogen concentration 
    plt = plot(ΔX, ΔY, line_z = ϕ, lw = 4, alpha = 0.7,
                c = cmap, label = "Hydra deformation", 
                legend = :bottomright);
    plot!(-ΔX, ΔY, line_z = ϕ, lw = 4, alpha = 0.7,
                c = cmap, label = "")

    # Plot material points 
    scatter!(ΔX, ΔY, ms = 2.0, color = :black, label = "")
    scatter!(-ΔX, ΔY, ms = 2.0, color = :black, label = "")

    # other plot necessities
    xlabel!("x"); ylabel!("y");
    # lm = 0.1 * r; xlims!(-lm, lm); ylims!(-lm, lm);

    return plt;
end

function visRdef(ϕ, X, Y)
    # plots deformation in radial direction against s, returning plot object

    # first, calculate Δr 
    ΔX = X .- X0; ΔY = Y .- Y0; # compute deformation
    Δr = ΔX .* rXu .+ ΔY .* rYu; # dot-product to get scalar radial deflection

    # then, plot dr against S 
    plt = plot(Si, Δr, label = L"\Delta r", lw = 2, line_z = ϕ, c = cmap,
                xlabel = "s", ylabel = "r", legend = :topright, colorbar = false,
                legendfontsize = 14);
    return plt;
end

function visQtys(ϕ, ϵsq)
    # plots morphogen and strain against S, returning plot object 
    plt = plot(Si, ϕ .- ϕ0, label = L"Δ\varphi", lw = 1.5, line_z = ϕ, c = cmap, ls = :dot, 
                xlabel = "s", ylabel = "Morphogen", legend = :topright, colorbar = false,
                legendfontsize = 14);
    plot!(twinx(), Si, ϵsq, label = L"\epsilon^2", lw = 2, line_z = ϕ, c = cmap,
            ylabel = "Strain", legend = :right, colorbar = false,
            legendfontsize = 14)
    return plt;
end

function visualise(ϕ, X, Y, ϵsq, titleTxt = false)
    # wrapper visualisation function to display all relevant plots at once 
    pltA = visShape(ϕ, X, Y); 
    pltB = visRdef(ϕ, X, Y);
    pltC = visDeform(ϕ, X, Y);
    pltD = visQtys(ϕ, ϵsq);
    combined = plot(pltA, pltB, pltC, pltD, layout = (2,2), size=(800,600)); 
    if titleTxt isa String
        title!(titleTxt)
    end
    return combined;
end

function postPlot(Tn, ϕtot, ϵtot, Ncutoff)
    # plots on the same plot, the total morphogen and total strain 
    # returns the plot object 
    plt = plot(Tn[1:Ncutoff], ϕtot[1:Ncutoff], label = L"\varphi", 
                color = :blue, lw = 2,
                title = "Hydra dynamics", xlabel = "t", ylabel = "Total morphogen",
                legend = :topleft);
    plot!(twinx(), Tn[1:Ncutoff], ϵtot[1:Ncutoff], label = L"\epsilon^2", 
            color = :red, lw = 2,
            ylabel = "Total strain", legend = :left)
    return plt;
end

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
X .= X0; Y .= Y0; Xdash .= X0dash; Ydash .= Y0dash;
# X .= X0test; Y .= Y0test; Xdash .= X0testDash; Ydash .= Y0testDash

# now morphogen array
ϕ = zeros(Ndisc);
ϕ .= ϕ0topBump;
# ϕ .= ϕ0sideBump;
# ϕ .= ϕ0; 

# save initial strain squared
ϵsq = trStrSq(X, Y, Xdash, Ydash)

# prepare the animation
fullAnim = Animation();
shapeAnim = Animation();
simName = "temp"; # "temp"
runAnim = true;

runsim = true;
tstart = time();
if runsim
    global ϵsq; 
    SaveData!(0, ϕ, X, Y, Xdash, Ydash, ϵsq,
                    ϕtot, ϵtot)

    for n = 1:Ntimes-1
        t = Tn[n];

        # Step A1: update shape of X(s), Y(s) and derivatives 
        shapeRes = UpdateShape!(ϕ, X, Y, Xdash, Ydash); 

        # Step A2: compute new strain tensor squared 
        ϵsq = trStrSq(X, Y, Xdash, Ydash); 

        # Step A3: evolve morphogen
        morphRes = UpdateMorphogen!(ϕ, X, Y, ϵsq);

        # Step B1: project onto eigenvector 
        ProjectEvec!(ϕ, X, Y)

        # Step C1: save necessary data 
        SaveData!(n, ϕ, X, Y, Xdash, Ydash, ϵsq,
                    ϕtot, ϵtot)

        # Step C2: visualise 
        tstr = round(t, digits = 2); 
        frameFullAnim = visualise(ϕ, X, Y, ϵsq,"t = $tstr");
        frameShapeAnim = visShape(ϕ, X, Y, "t = $tstr");
        if runAnim
            frame(fullAnim, frameFullAnim)
            frame(shapeAnim, frameShapeAnim) 
        end
        # if mod(n-1, plotRes) == 0 || n == Ntimes-1 # always visualise first and last 
            # display(shapeAnim);
        # end

        # Step D: raise a warning if either of the updates failed.
        if !Optim.converged(shapeRes)
            println("$shapeRes \n Shape update failed to converge: Step $n")
            global Ncutoff = n;
            break
        end
        if !Optim.converged(morphRes)
            println("$morphRes \n Morphogen update failed to converge: Step $n")
            global Ncutoff = n;
            break
        end

        println("$n, $(round(time() - tstart)), $(shapeRes.iterations)")
    end

    # output finish data 
    println("Total animation time: $(round(time() - tstart)) seconds for $Ncutoff timesteps.")

    # save and plot relevant objects
    mp4(fullAnim, simFP * "/sim_$simName-full.mp4", fps = 4);
    mp4(shapeAnim, simFP * "/sim_$simName-shape.mp4", fps = 4);
    plt = postPlot(Tn, ϕtot, ϵtot, Ncutoff); display(plt); 
end


############ -------------------------------------------------- ############
############ ---------------- TESTING SHAPE UPDATES ------------------ ############
############ -------------------------------------------------- ############
# X .= X0test; Y .= Y0test; Xdash .= X0testDash; Ydash .= Y0testDash
# X .= X0; Y .= Y0; Xdash .= X0dash; Ydash .= Y0dash;
# ϕ .= ϕ0sideBump;
# ϵsqtot = ElEn(ϕ, X, Y, Xdash, Ydash)
# plt = visShape(ϕ, X, Y,"before"); display(plt)
# # println("Elastic energy: ", ϵsqtot)
# # println("Total volume energy: ", Ω * (VolInt(X, Y, Xdash, Ydash) - V0)^2)

# res = UpdateShape!(ϕ, X, Y, Xdash, Ydash)
# plt = visShape(ϕ, X, Y, "after"); display(plt)

# ProjectEvec!(ϕ, X, Y)
# plt = visShape(ϕ, X, Y, "after proj"); display(plt)
