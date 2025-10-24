

using LinearAlgebra, DifferentialEquations, BandedMatrices; # calculations
using Plots, LaTeXStrings, Measures, Subscripts; # plots
import GLMakie, Makie, FileIO; # 3D plots 
using Optim, ForwardDiff; # optimisation
using JLD2; # loading and saving data 
using Polynomials; # testing 
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
F0 = 440.0; # base Young's modulus (440). NOTE: this interacts with Ω
h = 20.0 # bilayer thickness (20). NOTE: this interacts with Ω
b = 5.0; # exponential decrease of stiffness (5). NOTE: this interacts with Ω
κ = 1.0; # TODO morphogen stress upregulation coefficient (1 for now)  
ζ = 0.2; # morphogen decay rate (0.2)  
D = 1000.0; # morphogen diffusion coefficient (0.01) 

# PARAM SETS OF NOTE: 
# - (A) 

F(ϕ) = F0 * h * exp(-b*ϕ); # elasticity function
f(ϵsq) = κ * ϵsq; # strain-dependent morphogen expression function 

############ -------------------------------------------------- ############
############ ------------- NUMERICAL PARAMETERS --------------- ############
############ -------------------------------------------------- ############
Ndisc = 40; # number of discretisation points on ξ (40)
ξmin = -π/2; ξmax = π/2; # bounds of ξ values (-π/2, π/2)
dt = 0.03; # time discretisation (0.03) 
tmax = 80*dt; # max time (400 * dt for evec). 
Ω = 1e1; # not too large number: punishing potential for volume deviation (1e2)
ω = 1e-1; # not too large number: surface friction (1e1 for base pin, 1e-1 for X-X0)
dξint = 0.01; # small number: distance inside the ξ grid to start at to avoid div0
# dξint = (π/Ndisc)/2;
maxIter = 15000; # iteratinos for the shape update (v expensive), 15000 usually 
lenEvalEst = 10; # how many final iterations to average eigenvalue over 

cfl = 2 * D * dt / (π/Ndisc)^2; # cfl condition for diffusion, checking sensible. ds is approximate 
# println("CFL number (should be <<1): $(round(cfl, digits = 5))")

############ -------------------------------------------------- ############
############ ----------- VISUALISATION PARAMETERS ------------- ############
############ -------------------------------------------------- ############
plotRes = 1; # how many timesteps per plot (Inf for just 1st and last plots)
cmap = cgrad(:viridis); # colormap for morphogen concentration 
αboost = 2; # factor to plot the eigenvector (1 to plot it exactly as seen) (2)
spacer() = plot(legend=false, framestyle = :none, grid=false, xlims=(0,.1), ylims=(0,.1))

############ -------------------------------------------------- ############
############ ---------------- NUMERICAL SETUP ----------------- ############
############ -------------------------------------------------- ############
avg(vec)=[(vec[i]+vec[i+1])/2 for i=1:length(vec)-1]

ξi = range(ξmin+dξint, ξmax-dξint, Ndisc); # grid of ξ values 
ξiCt = avg(ξi); # svalscentre
dξ = ξi[2] - ξi[1]; 
Tn = 0:dt:tmax; # grid of t values 
Ntimes = length(Tn); 
Ncutoff = Ntimes; # where the sim ends because it breaks
interior = 2:Ndisc-1; # index set for the interior of the array

# preparing quantities to be saved. Each qty is saved BEFORE and AFTER sim
ϕtot = zeros(Ntimes); 
ϵtot = zeros(Ntimes);
evalEstFull = zeros(lenEvalEst);

# quantities necessary for 3D 
Ndisc3D = 100; # discretisation level in 3D
θ = range(0, 2*π, length = Ndisc3D);


############ -------------------------------------------------- ############
############ ------------- DISCRETISED OPERATORS -------------- ############
############ -------------------------------------------------- ############

# first deriv - no boundary conditions, just a higher order forward/backward difference on the ends
# used for x deriv 
Px = Matrix( Tridiagonal(fill(-1.0, Ndisc-1), fill(0.0, Ndisc), fill(1.0, Ndisc-1)) ) # main section 
Px[1, 1] = -3; Px[1, 2] = 4; Px[1, 3] = -1; # forward diff top row 3rd order
Px[Ndisc, Ndisc] = 3; Px[Ndisc, Ndisc-1] = -4; Px[Ndisc, Ndisc-2] = 1; # backward diff bottom row 
# Px[1,1:2] = 2*[-1, 1] # trying first-order difference 
# Px[Ndisc, Ndisc-1:Ndisc] = 2*[-1, 1]
Px /= (2*dξ); # scale 

# first deriv - 0 bcs, used for y and morphogen. Top and bottom row = 0
# Py = Matrix( Tridiagonal(fill(-1.0, Ndisc-1), fill(0.0, Ndisc), fill(1.0, Ndisc-1)) );
# Py[1, :] .= 0.0; Py[Ndisc, :] .= 0.0
# Py /= (2*dξ);
Py = Px;

# first deriv periodic bcs 
# Py = Matrix( Tridiagonal(fill(-1.0, Ndisc-1), fill(0.0, Ndisc), fill(1.0, Ndisc-1)) );
# Py[1, :] .= 0.0; Py[Ndisc, :] .= 0.0;
# Py[1, 2] = 1.0; Py[1, end] = -1.0; 
# Py[Ndisc, 1] = 1.0; Py[Ndisc, Ndisc-1] = -1.0;
# Py /= (2*dξ);

# convert to functions 
FDX(J) = Px*J;
FDY(J) = Py*J;
FDCT(J) = diff(J)/dξ; # centered first derivative


# integral over total ξ domain, using trapezoid rule
integrateξdom(Qty) = dξ * ( sum(Qty) - (Qty[1] + Qty[end])/2 );
integrateξdomCt(QtyCt) = Ndisc/(Ndisc-1) * integrateξdom(QtyCt); # centered needs correction

# inner product for eigenvector calculation 
inp(ZZ1, ZZ2) = integrateξdom((ZZ1[:,1].*ZZ2[:,1] .+ ZZ1[:,2].*ZZ2[:,2]/r^2 
                    .+ ZZ1[:,3].*ZZ2[:,3]/r^2 ) .* volEl)

############ -------------------------------------------------- ############
############ ------------- PRESAVING CALCULATIONS ------------- ############
############ -------------------------------------------------- ############
Cosξ = cos.(ξi); Sinξ = sin.(ξi); # presaving trigs 
volEl = R0^2 * Cosξ; # volume element

# centered versions
CosξCt = cos.(ξiCt);
volElCt = R0^2 * CosξCt;

# 3D requirement 
Cosθ = cos.(θ); Sinθ = sin.(θ);

############ -------------------------------------------------- ############
############ -------- STEADY-STATE / INITIAL FUNCTIONS -------- ############
############ -------------------------------------------------- ############
ϵ0sc = 0.5 * ((r^2 - R0^2)/R0^2)^2; # scalar initial (constant) trace of strain tensor squared
ϵ0 = zeros(Ndisc) .+ ϵ0sc; # over ξ-grid
V0 = 4/3 * r^3; # initial volume, normalised by pi

ϕ0sc = 1/ζ * f.(ϵ0sc); # steady-state morphogen scalar
ϕ0sc = (ζ==0) ? 1.0 : ϕ0sc;
ϕ0 = zeros(Ndisc) .+ ϕ0sc; # steady-state morphogen over the full grid 

X0 = r * Cosξ; Y0 = r * Sinξ; # initial condition shape 
X0dash = -Y0[:]; Y0dash = X0[:]; # initial derivatives 

Xundef = R0 * Cosξ; Yundef = R0 * Sinξ; # undeformed shape
rXu =  X0 ./ r; rYu =  Y0 ./ r; # radial unit vectors

trE0sc = (r^2-R0^2)/R0^2; # scalar initial (constant) trace of strain tensor
trE0 = zeros(Ndisc) .+ trE0sc;

trσ0sc = trE0sc * F(ϕ0sc); # scalar initial (constant) trace of stress tensor 
trσ0 = zeros(Ndisc) .+ trσ0sc; # over ξ-grid

trb0sc = -2/r; # scalar initial (constant) trace of curvature tensor 
trb0 = zeros(Ndisc) .+ trb0sc;

ZZ0 = hcat(ϕ0, X0, Y0); # full state vector, as Ndisc * 3 matrix 

# -- 3D ICs --- #
X03D = [X0[i] * Cosθ[j] for i in 1:Ndisc, j in 1:Ndisc3D];
Y03D = [Y0[i] * 1       for i in 1:Ndisc, j in 1:Ndisc3D];
Z03D = [X0[i] * Sinθ[j] for i in 1:Ndisc, j in 1:Ndisc3D];
ϕ03D = [ϕ0[i] * 1       for i in 1:Ndisc, j in 1:Ndisc3D];

############ -------------------------------------------------- ############
############ ------------ FUNCTIONALS AND QTYs --------------- ############
############ -------------------------------------------------- ############
function trStrSq(X, Y, Xdash, Ydash)
    # takes vectors of the surface parameterisation x(ξ), y(ξ) and derivatives w.r.t. s
    # returns a vector containing the trace of the square of the strain tensor at each s
    t1 = (Xdash).^2 .+ (Ydash).^2 .- R0^2;
    # handle boundary different if ξ coords go all the way to the boundary 
    # if Si[1] == -π/2
    #     t2[interior] = X[interior].^2 ./ (Cosξ[interior].^2); # interior is just x/cos^2 
    #     t2[[1,end]] .= r^2; # boundary is when x(ξ)~r cos(ξ)
    # else
    #     t2 = X.^2 ./ Cosξ.^2;
    # end
    # t2 = t2 .- R0^2; 
    t2 = (X.^2 ./ Cosξ.^2) .- R0^2;

    return 1/(4*R0^4) * (t1.^2 .+ t2.^2); 
end

function trStrain(X, Y, Xdash, Ydash)
    # trace of the strain tensor, ϵ_α^α
    return 1/(2*R0^2) .* (Xdash.^2 .+ Ydash.^2 .+ (X.^2)./(Cosξ.^2) .- 2*R0^2);
end

function trStress(Fϕ, X, Y, Xdash, Ydash)
    # trace of the stress tensor, σ_α^α = F(ϕ)ϵ_α^α 
    return Fϕ .* trStrain(X, Y, Xdash, Ydash);
end

function trCurv(X, Y, Xdash, Ydash)
    # trace of the curvature tensor 
    rsq = Xdash.^2 .+ Ydash.^2;
    d1 = Ydash.*FDX(Xdash) .- Xdash.*FDY(Ydash); # denominator of first term 
    d2 = Ydash;
    return d1./(rsq.^(3/2)) .- d2./(X.*sqrt.(rsq));
end

function ElEn(Fϕ, X, Y, Xdash, Ydash)
    # Elastic energy functional 
    return integrateξdom(Fϕ/2 .* trStrSq(X, Y, Xdash, Ydash) .* volEl); # integrand over S0 including volel
end

function VolInt(X, Ydash)
    # volume integral 
    fn = X.^2 .* Ydash; # integrand normalised by pi 
    return integrateξdom(fn);
end

function SurfaceFriction(X, Y)
    # number that models friction to stop entire hydra running away

    # Keep it close to the IC 
    return sum( (Y .- Y0).^2 + (X .- X0).^2 );
end

κ_bend = 0.0;
function BendingPenalty(X)
    # try to stop bending by minimising derivative 
    ΔX = X .- X0;
    return sum( (FDX(ΔX)).^2 );
end

κ_fric = 0.0;
function Friction(X, Y)
    ΔX = X .- dZZ[:,2];
    ΔY = Y .- dZZ[:,3];
    return sum(ΔX.^2 + ΔY.^2);
end

function EnergyFunctional(Fϕ, X, Y, Xdash, Ydash)
    # total energy functional, including elastic energy and a punishing potential and surface friction
    return ElEn(Fϕ, X, Y, Xdash, Ydash) + Ω * (VolInt(X, Ydash) - V0)^2 +
                ω * SurfaceFriction(X, Y);
                # κ_fric * Friction(X, Y) + 
                # κ_bend * (BendingPenalty(X) + BendingPenalty(Y));
end

function EFwrapped(fullData, Fϕ)
    # wrapper for the energy functional. Unpacks the data, computes derivatives, then
    # returns the energy functional to optimise 
    # trying to minimise extra space allocations
    return EnergyFunctional(Fϕ, fullData[1:Ndisc], fullData[Ndisc+1:end], # X, Y 
                                FDX(fullData[1:Ndisc]), FDY(fullData[Ndisc+1:end]));
end

function MorphogenFunctional(ϕ, ϕn, X, Y, ϵsq)
    # total energy for the morphogen equation 
    # ϕn is the previous state 
    integrand1 = (ϕ .- ϕn).^2 /(2*dt) .* volEl; # time-deriv bit
    integrand2 = -f.(ϵsq) .* ϕ .* volEl; # morphogen production bit
    integrand3 = ζ * (ϕ).^2 /2 .* volEl; # disassociation bit
    integrand4 = D * FDCT(ϕ).^2 /2 .* CosξCt; # diffusion bit, volEl already incorporated
    
    integrandTot = integrand1 .+ integrand2 .+ integrand3;
    return integrateξdom(integrandTot) + integrateξdom(integrand4);
end

############ -------------------------------------------------- ############
############ ------------------ UPDATE STEPS ------------------ ############
############ -------------------------------------------------- ############
function Renormalise!(X, Y, Xdash, Ydash, V0)
    # takes a shape parameterisation, and modifies in-place to ensure that the volume 
    # remains equal to some V0. 
    # uses the fact that integral propto 3 factors of x, y, x', y', so a cube 
    # root normalisation constant can apply equally to x,y,x',y' and ensure volume works
    vtemp = VolInt(X, Ydash);
    normk = cbrt(V0/vtemp);
    X .= X * normk; Y .= Y * normk; 
    Xdash .= Xdash * normk; Ydash .= Ydash * normk;
end 

function UpdateShape!(Fϕ, X, Y, Xdash, Ydash)
    # takes the current state data (shape and morphogen concentration), and updates X, Y 
    # and derivatives IN PLACE, by minimising the energy functional 
    # returns results of the optim

    # currently no BCs in the optim, but these could be useful if implemented 
    # lbounds = [dsint * 0.9; Vector(1:Ndisc-2)*-Inf; dsint*0.9; Vector(1:Ndisc)*-Inf ];
    # ubounds = [dsint * 1.1; Vector(1:Ndisc-2)*Inf; dsint*1.1; Vector(1:Ndisc)*Inf];
    
    result = optimize(x -> EFwrapped(x, Fϕ), 
                    [X; Y], 
                    # lbounds, ubounds, Fminbox(GradientDescent()),
                    method = LBFGS(), 
                    autodiff = :forward, iterations = maxIter);  
    X .= result.minimizer[1:Ndisc]; Y .= result.minimizer[Ndisc+1:end];

    # finally, update derivatives 
    Xdash .= FDX(X); Ydash .= FDY(Y);
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
    ϕtot[n+1] = integrateξdom(ϕ); 
    ϵtot[n+1] = integrateξdom(ϵsq);
end

############ -------------------------------------------------- ############
############ --------------- EIGENVECTOR STEPS ---------------- ############
############ -------------------------------------------------- ############
function ProjectEvec!(ϕ, X, Y; dEVs = [], dEVnorms = [])
    # renormalises state to allow it to slowly project onto the leading eigenvector 
    # modifies in place, replacing the existing state with one that has the perturbation normalised
    # returns the normalised perturbation 
    # also returns a flag if it 'converged' correctly (original norm close to target norm)
    # finally, returns estimated eigenvalue of time evolution
    # optional argument to first subtract projections onto each eigenvector in EV
    # takes the deviation from steady-state, not the full state/
    ZZ = hcat(ϕ, X, Y) # full state 
    dZZ = ZZ .- ZZ0; # change in state 
    for (dEv, dEvnorm) in zip(dEVs, dEVnorms) # subtract projection onto each existing eigenvector 
        coeff = inp(dZZ, dEv) / dEvnorm; # coefficient of projection 
        dZZ .-= coeff .* dEv; # subtract projection
    end
    
    globalNorm = sqrt(inp(dZZ, dZZ)); # norm over entire domain 
    # println("Global norm: $(round(globalNorm, digits = 4))") 
    targetNorm = 6; # target norm 
    projRes = ( (globalNorm / targetNorm) <= 10 ) 
    # projRes = true;

    evalEst = log((globalNorm/targetNorm))/dt;

    dZZ = dZZ / globalNorm * targetNorm; # normalise the perturbation 
    ZZ = ZZ0 .+ dZZ; # get the full state again 
    ϕ .= ZZ[:,1]; X .= ZZ[:,2]; Y .= ZZ[:,3]; # save new state 
    return (dZZ, projRes, evalEst);
end;

############ -------------------------------------------------- ############
############ ------------------ BETTER ICs -------------------- ############
############ -------------------------------------------------- ############

# test approximations of r cos(ξ), r sin(ξ) to use as initial guesses for x(ξ), y(ξ) to test
# shape evolution. 
# To ensure BCs, using quadratic approximation of x(ξ) with roots at +-pi/2 plus an interesting perturbation 
# Just using the actual y(ξ) function 
x0temp(ξ) = r * ( 1 - 4/π^2 * ξ^2 + 0.2 * (ξ^2 - π^2/4) * ξ );
y0temp(ξ) = r * sin(ξ);
x0tempDash(ξ) = r * ( -8/π^2 * ξ + 0.2 * (3*ξ^2 - π^2/4) );
y0tempDash(ξ) = r * cos(ξ); 

# so that the guess isn't too far off the actual soln, construct the test as a weighted sum of 
# the correct IC and the 'interesting' IC, then adjust for correct volume
α_shape = 0.2; # weight in the 'interesting' initial condition (0.2)
X0test = α_shape*x0temp.(ξi) + (1-α_shape)*X0; Y0test = α_shape*y0temp.(ξi) + (1-α_shape)*Y0;
X0testDash = α_shape*x0tempDash.(ξi) + (1-α_shape)*X0dash; Y0testDash = α_shape*y0tempDash.(ξi) + (1-α_shape)*Y0dash;
Renormalise!(X0test, Y0test, X0testDash, Y0testDash, V0);

### --- ϕ ICs --- ###
# nonuniform but obey BCs
α_morph = 0.125; # weight of the nonuniform 'interesting' component (0.125)
# bump on the top of the domain
ϕ0tempTop = zeros(Ndisc) .+ Sinξ * ϕ0sc .+ ϕ0sc; 
ϕ0topBump = α_morph*ϕ0tempTop .+ (1-α_morph)*ϕ0; 

# bump on the side of the domain 
ϕ0tempSide = zeros(Ndisc) .+ 2*exp.(-2 * ξi.^2) * ϕ0sc;
ϕ0sideBump = α_morph*ϕ0tempSide .+ (1-α_morph)*ϕ0; 

# bump on the bottom of the domain 
ϕ0tempBottom = zeros(Ndisc) .+ Sinξ * ϕ0sc .+ ϕ0sc; 
ϕ0bottomBump = -α_morph*ϕ0tempBottom .+ (1-α_morph)*ϕ0; 

# double bump 
ϕ0tempDouble = zeros(Ndisc) .+ sin.(2*ξi).^2 * ϕ0sc .+ ϕ0sc;
ϕ0doubleBump = α_morph*ϕ0tempDouble .+ (1-α_morph)*ϕ0;

# very random not odd or even 
ϕ0tempBumpy = zeros(Ndisc) .+ (sin.(4*ξi).^2 .+ 5*cos.(ξi .+ 0.5) .+ 3)/8 * ϕ0sc;
ϕ0bumpy = α_morph * ϕ0tempBumpy .+ (1-α_morph)*ϕ0;

############ -------------------------------------------------- ############
############ ------------ VISUALISATION FUNCTIONS ------------- ############
############ -------------------------------------------------- ############
function visShape(ϕ, X, Y, titleTxt = false; αboost = 1)
    # takes a surface discretised parameterisation and morphogen concentration
    # displays a plot of the deformed hydra, colored by the morphogen concentration 
    # optionally adds a provided title 
    # returns the plot object for displaying and saving
    # αboost: enhances perturbation, set to 1 for no perturbation
    X = X0 .+ αboost*(X .- X0); Y = Y0 .+ αboost*(Y .- Y0);

    ticks = collect(range(minimum(ϕ), maximum(ϕ), 5))

    # Plot actual deformed shape, colored by morphogen concentration 
    plt = plot(X, Y, line_z = ϕ, lw = 5, alpha = 0.7,
                c = cmap, colorbar_title = "φ",
                label = "", 
                aspect_ratio = :equal, legend = :topright);
    plot!(-X, Y, line_z = ϕ, lw = 3, alpha = 0.7, c = cmap, label = "")

    # Plot material points 
    scatter!(X, Y, ms = 1.5, color = :black, label = "")
    scatter!(-X, Y, ms = 1.5, color = :black, label = "")

    # dummy bit to get hydra shape in legend 
    scatter!([1000, 2000], [1000, 2000], ms = 0.3, color = :black, label = "Hydra shape")

    # plot initial shape 
    plot!(X0, Y0, lw = 1, color = :grey, ls = :dash, label = "Steady-state"); 
    plot!(-X0, Y0, lw = 1, color = :grey, ls = :dash, label = "");

    # dummy bit to get a good color scale 
    # if ϕ0[1] != Inf
    #     plot!([0.0, 0.01], [0.0, 0.01], lw = 0, 
    #             line_z = [0.0, 2*ϕ0[1]], c = cmap, label = "");
    # end

    # other plot necessities
    xlabel!("x (μm)"); ylabel!("y (μm)");
    lm = 1.33 * r; xlims!(-lm, lm); ylims!(-lm, lm);
    if titleTxt isa String
        title!(titleTxt)
    end

    return plt;
end

function visShapeSimple(ϕ, X, Y, titleTxt = false; αboost = 1, ψ = 0, undef = false)
    # same as above but simpler
    # allows to rotate by psi (in degrees)
    # allows to plot undeformed shape as well 
    X = X0 .+ αboost*(X .- X0); Y = Y0 .+ αboost*(Y .- Y0);
    X = [X;-reverse(X)]; Y = [Y;reverse(Y)]; ϕ = [ϕ;reverse(ϕ)] # create full shape 
    Y = Y.+7;
    
    ψ = ψ * π/180 .+ π/2; # convert to rad and polar angle
    cψ = cos(ψ); sψ = sin(ψ)

    # rotate 
    rpts = sqrt.(X.^2 .+ Y.^2);
    θpts = atan.(X, Y);
    X = rpts .* cos.(θpts .+ ψ); Y = rpts .* sin.(θpts .+ ψ);

    plt = plot(X, Y, line_z = ϕ, lw = 6, alpha = 0.7,
                c = cmap, colorbar_title = "Morphogen concentration (φ)",
                colorbar_ticks = false, colorbar_ticklabels = false,   
                label = "", 
                aspect_ratio = :equal, legend = :bottomright,
                legendfontsize = 11, titlefontsize = 20, tickfontsize = 14, 
                guidefontsize = 14, colorbar_titlefontsize = 14,
                xticks = [], yticks = [],
                size = (460, 400),
    );

    # Plot material points 
    scatter!(X, Y, ms = 2, color = :black, label = "")

    # dummy bit to get hydra shape in legend 
    scatter!([1000, 2000], [1000, 2000], ms = 0.3, color = :black, label = "Hydra shape")

    # plot initial shape 
    if undef 
        plot!(R0*Cosξ, R0*Sinξ, lw = 1, color = :grey, ls = :dash, label = "Reference surface"); 
        plot!(-R0*Cosξ, R0*Sinξ, lw = 1, color = :grey, ls = :dash, label = "");
    else
        plot!(X0, Y0, lw = 1, color = :grey, ls = :dash, label = "Steady-state"); 
        plot!(-X0, Y0, lw = 1, color = :grey, ls = :dash, label = "");
    end

    xlabel!("x"); ylabel!("y");
    lm = 1.33 * r; xlims!(-lm, lm); ylims!(-lm, lm);
    
    title!(titleTxt)

    return plt;
end

function plotHSS()
    # simplifies to just plot homogeneous steady-state
    X = X0; Y = Y0; ϕ = ϕ0; # set to steady-state
    X = [X;-reverse(X)]; Y = [Y;reverse(Y)]; ϕ = [ϕ;reverse(ϕ)] # create full shape 

    plt = plot(X, Y, line_z = ϕ, lw = 6, alpha = 0.7,
                c = cmap, colorbar_title = "Morphogen concentration (φ)",
                colorbar_ticks = false, colorbar_ticklabels = false,   
                label = "", 
                aspect_ratio = :equal, legend = :topright,
                legendfontsize = 11, titlefontsize = 20, tickfontsize = 14, 
                guidefontsize = 14, colorbar_titlefontsize = 14,
                xticks = [], yticks = [],
                size = (460, 400),
                yaxis = false
    );

    # Plot material points 
    scatter!(X, Y, ms = 2, color = :black, label = "")

    # dummy bit to get hydra shape in legend 
    scatter!([1000, 2000], [1000, 2000], ms = 0.3, color = :black, label = "Hydra shape")
 
    plot!(R0*Cosξ, R0*Sinξ, lw = 1, color = :grey, ls = :dash, label = "Reference surface"); 
    plot!(-R0*Cosξ, R0*Sinξ, lw = 1, color = :grey, ls = :dash, label = "");

    # dummy ahh y axis 
    vline!([0.0], lw = 1.3, ls = :solid, label = "", color = :black)
    annotate!(20, 0, ("y", 14, :black, :center))

    xlabel!("x"); # ylabel!("y");
    lm = 1.33 * r; xlims!(-lm, lm); ylims!(-lm, lm);

    return plt;

end

function visDeform(ϕ, X, Y)
    # takes a surface discretised parameterisation and morphogen concentration
    # displays a plot of just the deformation, i.e. the difference between main shape and IC

    ΔX = X .- X0; ΔY = Y .- Y0; # compute deformation

    # Plot deformation, colored by morphogen concentration 
    plt = plot(ΔX, ΔY, line_z = ϕ, lw = 4, alpha = 0.7,
                c = cmap, label = "Deformation " * L"\mathbf{x}-\mathbf{x}_0", 
                legend = :topleft);
    plot!(-ΔX, ΔY, line_z = ϕ, lw = 4, alpha = 0.7,
                c = cmap, label = "")

    # Plot material points 
    scatter!(ΔX, ΔY, ms = 2.0, color = :black, label = "")
    scatter!(-ΔX, ΔY, ms = 2.0, color = :black, label = "")

    # other plot necessities
    xlabel!("x (μm)"); ylabel!("y (μm)");
    # lm = 0.1 * r; xlims!(-lm, lm); ylims!(-lm, lm);

    return plt;
end

function visVecDeform(ϕ, X, Y)
    # visualises defomration as a bunch of vectors plotted from the origin,
    # colored by ξ. 
    ΔX = X .- X0; ΔY = Y .- Y0; # compute deformation
    og = 0 * ΔX; # origin

    plt = quiver(og, og, quiver = (ΔX, ΔY),
                    line_z = repeat(ξi, inner=4), c = :magma,
                    linewidth = 2,
                    colorbar_title = "ξ",
                    );
    xlabel!("Δx (μm)"); ylabel!("Δy (μm)");

    return plt;
end

function visRdef(ϕ, X, Y)
    # plots deformation in radial direction against s, returning plot object

    # first, calculate Δr 
    ΔX = X .- X0; ΔY = Y .- Y0; # compute deformation
    Δr = ΔX .* rXu .+ ΔY .* rYu; # dot-product to get scalar radial deflection

    # then, plot dr against ξ 
    plt = plot(ξi, Δr, label = L"|\Delta \mathbf{r}|", lw = 2, line_z = ϕ, c = cmap,
                xlabel = "ξ", ylabel = "r (μm)", legend = :topleft, colorbar = false,
                legendfontsize = 14);
    return plt;
end

function visQtys(ϕ, ϵsq, X, Y,
                    titleTxt = false;
                    plotE = true, plotTrig = false, plotStress = false)
    # plots morphogen and strain against S, returning plot object 
    # additional functionality to plot a trigonometric function along with it to compare 
    # spherical harmonics 
    # finally, additinoal functionality to calculate and plot stress 
    intr = 2; # how many edge points to exclude
    ξcut = ξi[intr+1:end-intr];

    plt = plot(ξcut, ϕ[intr+1:end-intr], label = L"\varphi", lw = 1, line_z = ϕ, c = cmap, 
                xlabel = "ξ", ylabel = "Morphogen", legend = :topleft, colorbar = false,
                legendfontsize = 12);
    if plotE
        plot!(twinx(), ξcut, ϵsq[intr+1:end-intr], label = L"E^2", lw = 2, color = :black, ls = :dash,
                ylabel = "Strain", legend = :left,
                legendfontsize = 12)
    end

    if plotStress
        trσ = trStress(F.(ϕ), X, Y, FDX(Xdash), FDY(Ydash));
        plot!(twinx(), ξi[2:end-1], trσ[2:end-1], label = L"\sigma_{\alpha}^{\alpha}", lw = 2, color = :black,
                ylabel = "Stress", legend = :left, colorbar = false,
                legendfontsize = 12)
    end

    if plotTrig
        # sinScaled = -Cosξ * (maximum(ϕ)-minimum(ϕ)) .+ (maximum(ϕ)+minimum(ϕ))/2;
        sinScaled = -cos.(2*ξi) * (maximum(ϕ)-minimum(ϕ))/2;
        sinScaled .-= sinScaled[end÷2] - ϕ[end÷2];
        plot!(ξi, sinScaled, label = L"-\cos(2ξ)", lw = 1, color = :red)
    end

    if titleTxt isa String
        title!(titleTxt)
    end
    return plt;
end

function visDqtys(ϕ, X, Y, titleTxt = false;
                    plotTrE = false, plotStress = false, plotCurv = false)
    # same as above but only plots perturbation in qtys, and stacks plots 

    # first calculate necessary quantities 
    Xdash = FDX(X); Ydash = FDY(Y);

    dϕ = ϕ .- ϕ0;
    trE = trStrain(X, Y, Xdash, Ydash);
    dtrE = trE .- trE0;
    trσ = trStress(F.(ϕ), X, Y, Xdash, Ydash);
    dtrσ = trσ .- trσ0; 
    trb = trCurv(X, Y, Xdash, Ydash);
    dtrb = trb .- trb0;

    plts = Vector{Any}();
    nplots = 1;

    intr = 2; # how many edge points to exclude
    ξcut = ξi[intr+1:end-intr];
    
    plt = plot(ξcut, dϕ[intr+1:end-intr]/ϕ0sc*100, label = L"\delta\varphi", lw = 2, color = nplots;
                xlabel = "", ylabel = "% φ", legend = :left,
                legendfontsize = 10);
    plot!(ξi, [0 for i in ξi], label = "", color = :grey, linestyle = :dash)
    # if titleTxt isa String
    #     title!(titleTxt)
    # end
    push!(plts, plt); nplots += 1;

    if plotTrE
        plt = plot(ξcut, dtrE[intr+1:end-intr]/trE0sc*100, label = L"\delta E_{\alpha}^{\alpha}", lw = 2, color = nplots;
                xlabel = "", ylabel = "% E", legend = :left,
                legendfontsize = 10);
        plot!(ξi, [0 for i in ξi], label = "", color = :grey, linestyle = :dash)
        push!(plts, plt); nplots += 1;
    end

    if plotStress
        plt = plot(ξcut, dtrσ[intr+1:end-intr]/trσ0sc*100, label = L"\delta\sigma_{\alpha}^{\alpha}", lw = 2, color = nplots;
                xlabel = "", ylabel = "% σ", legend = :left,
                legendfontsize = 10);
        plot!(ξi, [0 for i in ξi], label = "", color = :grey, linestyle = :dash)
        push!(plts, plt); nplots += 1;
    end

    if plotCurv
        plt = plot(ξcut, dtrb[intr+1:end-intr]/trb0sc*100, label = L"\delta b_{\alpha}^{\alpha}", lw = 2, color = nplots;
                xlabel = "ξ", ylabel = "% b", legend = :left,
                legendfontsize = 10);
        plot!(ξi, [0 for i in ξi], label = "", color = :grey, linestyle = :dash)
        push!(plts, plt); nplots += 1;
    end
    
    wd = 660;
    ht = 125 * (nplots-1);
    pltbig = plot(plts...; layout = (nplots-1, 1), size = (wd, ht))

    # experimenting 
    # pexp = scatter(dϕ[intr+1:end-intr], dtrE[intr+1:end-intr], 
    #                 xlabel = L"\delta\varphi", ylabel = L"\delta E_{\alpha}^{\alpha}",
    #                 legend = false)
    # scatter!([0],[0], color = 2);
    # if titleTxt isa String
    #     title!(titleTxt)
    # end
    # linfit = fit(dϕ[intr+1:end-intr], dtrE[intr+1:end-intr], 1)   # degree-1 polynomial
    # print(coeffs(linfit));
    # display(pexp)


    return pltbig;
end



# function visualise(ϕ, X, Y, ϵsq, titleTxt = false; αboost = 1)
#     # wrapper visualisation function to display all relevant plots at once 
#     # αboost: increases perturbation. Set to 1 for no boost
#     X = X0 .+ αboost*(X .- X0); Y = Y0 .+ αboost*(Y .- Y0);
#     pltA = visShape(ϕ, X, Y); 
#     pltB = visRdef(ϕ, X, Y);
#     # pltC = visDeform(ϕ, X, Y);
#     pltC = visVecDeform(ϕ, X, Y);
#     pltD = visQtys(ϕ, ϵsq, X, Y, plotE = false, plotStress = true);
#     combined = plot(pltA, pltB, pltC, pltD, layout = (2,2), size=(800,600)); 
#     if titleTxt isa String
#         title!(titleTxt)
#     end
#     return combined;
# end

function visualise(ϕ, X, Y, ϵsq, titleTxt = false; αboost = 1)
    # wrapper visualisation function to display all relevant plots at once 
    # αboost: increases perturbation. Set to 1 for no boost
    Xb = X0 .+ αboost*(X .- X0); Yb = Y0 .+ αboost*(Y .- Y0);
    pltA = visShape(ϕ, Xb, Yb, titleTxt); 
    pltB = visDqtys(ϕ, X, Y; plotTrE = true, plotStress = false, plotCurv = false);
    combined = plot(pltA, pltB, layout = (1,2)); 
    return combined;
end

function postPlot(Tn, ϕtot, ϵtot, Ncutoff)
    # plots on the same plot, the total morphogen and total strain 
    # returns the plot object 
    plt = plot(Tn[1:Ncutoff], ϕtot[1:Ncutoff], label = L"\varphi", 
                color = :blue, lw = 2,
                title = "Hydra dynamics", xlabel = "t (hr)" * L"^{-1}", ylabel = "Total morphogen",
                legend = :topleft);
    plot!(twinx(), Tn[1:Ncutoff], ϵtot[1:Ncutoff], label = L"E^2", 
            color = :red, lw = 2,
            ylabel = "Total strain", legend = :left)
    return plt;
end

function vis3D(ϕ, X, Y, titleTxt = false; αboost3D = 1, deform = false)
    # creates and returns a 3D plot created by rotating the (X, Y) parameterised 
    # surface through 3D space 
    # αboost3D enhances the perturbation to make it more visible. Set to 1 for no boost
    # deform: only plots deformation as its own surface 

    # first, prepare the shape
    X3D = [X[i] * Cosθ[j] for i in 1:Ndisc, j in 1:Ndisc3D];
    Y3D = [Y[i] * 1       for i in 1:Ndisc, j in 1:Ndisc3D];
    Z3D = [X[i] * Sinθ[j] for i in 1:Ndisc, j in 1:Ndisc3D];
    ϕ3D = [ϕ[i] * 1       for i in 1:Ndisc, j in 1:Ndisc3D];

    # next, boost it 
    (dX3D, dY3D, dZ3D) = X3D .- X03D, Y3D .- Y03D, Z3D .- Z03D;
    (X3D, Y3D, Z3D) = X03D .+ αboost3D*dX3D, Y03D .+ αboost3D*dY3D, Z03D .+ αboost3D*dZ3D;

    # try sharper boost 
    # bt(x) = sign.(x) .* abs(x).^2;
    # (X3D, Y3D, Z3D) = X03D .+ αboost3D*bt.(dX3D), Y03D .+ αboost3D*bt.(dY3D), Z03D .+ αboost3D*bt.(dZ3D);

    # convert to just deformed if required
    X3D = deform ? dX3D : X3D;
    Y3D = deform ? dY3D : Y3D; 

    plt = GLMakie.surface(X3D, Y3D, Z3D,
                axis = (type = GLMakie.Axis3, azimuth = pi/4,
                           aspect = :equal),
                color = ϕ3D, colormap = cmap);

    if titleTxt isa String
        title!(titleTxt)
    end
    return plt;
end

############ -------------------------------------------------- ############
############ ------------- LOADING EIGENVECTORS --------------- ############
############ -------------------------------------------------- ############

@load "EVs//EV1.jld2" EV1;
dEV1 = EV1 .- ZZ0; # the difference to steady state 
dEV1Norm = inp(dEV1, dEV1);

@load "EVs//EV1r.jld2" EV1r;
dEV1r = EV1r .- ZZ0; 
dEV1rNorm = inp(dEV1r, dEV1r);

@load "EVs//EV2side.jld2" EV2side;
dEV2side = EV2side .- ZZ0; 
dEV2sideNorm = inp(dEV2side, dEV2side);

@load "EVs//EV2top.jld2" EV2top;
dEV2top = EV2top .- ZZ0; 
dEV2topNorm = inp(dEV2top, dEV2top);

@load "EVs//EV3r.jld2" EV3r;
dEV3r = EV3r .- ZZ0; 
dEV3rNorm = inp(dEV3r, dEV3r);


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
# ϕ .= ϕ0topBump;  
# ϕ .= ϕ0sideBump;
# ϕ .= ϕ0bottomBump;
# ϕ .= ϕ0doubleBump;
ϕ .= ϕ0bumpy;
# ϕ .= ϕ0; 

# save initial strain squared
ϵsq = trStrSq(X, Y, Xdash, Ydash)

# prepare the animation
fullAnim = Animation();
shapeAnim = Animation();
simName = "temp"; # "temp"
runAnim = true;

dZZ = zeros(Ndisc, 3); dϕ = zeros(Ndisc); dX = zeros(Ndisc); dY = zeros(Ndisc);
EV = zeros(Ndisc, 3);

runsim = false;
doProj = true; # whether to project vs just do time evolution 

tstart = time();
if runsim
    global ϵsq; 
    global dZZ, dϕ, dX, dY;
    SaveData!(0, ϕ, X, Y, Xdash, Ydash, ϵsq,
                    ϕtot, ϵtot)

    for n = 1:Ntimes-1
        t = Tn[n];

        # Step A1: update shape of x(ξ), y(ξ) and derivatives 
        shapeRes = UpdateShape!(F.(ϕ), X, Y, Xdash, Ydash); 

        # Step A2: compute new strain tensor squared 
        ϵsq = trStrSq(X, Y, Xdash, Ydash); 

        # Step A3: evolve morphogen
        morphRes = UpdateMorphogen!(ϕ, X, Y, ϵsq);

        # Step B1: project onto eigenvector 
        if doProj || n==1
            (dZZ, projRes, evalEst) = 
                ProjectEvec!(ϕ, X, Y; 
                    dEVs =      [], 
                    dEVnorms =  []) 
            dϕ = dZZ[:,1]; dX = dZZ[:,2]; dY = dZZ[:,3]; # unpack 
            if Ntimes - n <= lenEvalEst
                evalEstFull[Ntimes-n] = evalEst;
            end
        end

        # Step C1: save necessary data 
        SaveData!(n, ϕ, X, Y, Xdash, Ydash, ϵsq,
                    ϕtot, ϵtot)

        # Step C2: visualise 
        tstr = round(t, digits = 2); 

        if runAnim
            frameFullAnim = visualise(ϕ, X, Y, ϵsq, "t = $tstr", αboost = αboost); 
            frameShapeAnim = visShape(ϕ, X, Y, "t = $tstr", αboost = αboost); 
            frame(fullAnim, frameFullAnim)
            frame(shapeAnim, frameShapeAnim) 
        end

        # special code for getting midyear review plots. Show some frames, then all info, then blowup frame
        if n in [1, 37, 70, Ntimes-1]
            savefig(frameFullAnim, "$n.png")
        end

        # Step D: raise a warning if any of the updates failed.
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
        if doProj
            if !projRes 
                println("Global norm exceeded \n Eigenvector projection failed to converge: Step $n")
                global Ncutoff = n;
                break
            end
        end

        println("$n, $(round(time() - tstart)), $(shapeRes.iterations)") 
        # println("$n, $(round(time() - tstart))")
    end

    # output finish data 
    println("Total sim time: $(round(time() - tstart)) seconds for $Ncutoff timesteps.")
    if doProj
        println("Estimated eigenvalue: $(sum(evalEstFull)/lenEvalEst)")
    end

    # save and plot relevant objects
    if runAnim
        mp4(fullAnim, simFP * "/sim_$simName-full.mp4", fps = 4);
        mp4(shapeAnim, simFP * "/sim_$simName-shape.mp4", fps = 4);
    end

    # plt = postPlot(Tn, ϕtot, ϵtot, Ncutoff); display(plt); 

    if doProj
        pltLE = visualise(ϕ, X, Y, trStrSq(X, Y, Xdash, Ydash), αboost = αboost);
        display(pltLE); savefig(pltLE, simFP * "/sim_$simName-evec.png")
    end

    # EV2top = hcat(ϕ, X, Y);
    # @save "EVs//EV2top.jld2" EV2top;
end
# EV .= EV1;
# pltEVEC = visualise(EV[:,1], EV[:,2], EV[:,3], 
#                     trStrSq(EV[:,2], EV[:,3], FDX(EV[:,2]), FDX(EV[:,3])), 
#                     "Eigenvector 1", αboost = 2)
# display(pltEVEC); savefig(pltEVEC, "TEMP.png") # savefig(pltEVEC, "EVs//EV1.png")

# EV .= EV3r;
# pltEVEC3D = vis3D(EV[:,1], EV[:,2], EV[:,3], "Eigenvector", αboost3D = 10, deform = false);
# scr = display(pltEVEC3D);
# sleep(0.1); FileIO.save("snapshot.png", Makie.colorbuffer(scr)); # to save particular camera view etc.

############ -------------------------------------------------- ############
############ ---------------- TROUBLESHOOTING FULL STEPS ------------------ ############
############ -------------------------------------------------- ############
# X .= X0test; Y .= Y0test; Xdash .= X0testDash; Ydash .= Y0testDash
# X .= X0; Y .= Y0; Xdash .= X0dash; Ydash .= Y0dash;
# ϕ .= ϕ0topBump;

# shapeRes_test = UpdateShape!(E.(ϕ), X, Y, Xdash, Ydash); 
# plt = visShape(ϕ, X, Y, "shape"); display(plt)
# plt = visDeform(ϕ, X, Y); display(plt);
# # plt = plot(FDX(X.-X0), FDY(Y.-Y0)); display(plt);

# # Step A2: compute new strain tensor squared 
# ϵsq = trStrSq(X, Y, Xdash, Ydash); 

# # Step A3: evolve morphogen
# morphRes_test = UpdateMorphogen!(ϕ, X, Y, ϵsq);
# # plt = plot(ξi, ϕ.-ϕ0); display(plt);

# # Step B1: project onto eigenvector 
# (dZZ_test, projRes_test, eval_test) = ProjectEvec!(ϕ, X, Y);



############ -------------------------------------------------- ############
############ ---------------- OTHER TESTING CODE ------------------ ############
############ -------------------------------------------------- ############
# pt = visDqtys(ϕ, X, Y, "Time evolution, t = 3 steps, different IC";
#         plotTrE = true, plotStress = false);
# display(pt);

# EV .= EV1r;
# pt = visDqtys(EV[:,1], EV[:,2], EV[:,3], "EV1";
#         plotTrE = true, plotStress = true, plotCurv = true);
# display(pt);

# EV .= EV3r;
# pt = visualise(EV[:,1], EV[:,2], EV[:,3], 3, "EV3", αboost = 3)
# display(pt); savefig(pt, "EV3.png"); 

# EV .= EV1r;
# pt = visShapeSimple(EV[:,1], EV[:,2], EV[:,3], L"Z_{1}", αboost = 5, ψ = 180); # ψ = -130
# display(pt); savefig(pt, "Z1.png")

# pt = plotHSS()
# display(pt); savefig(pt, "HSS.pdf")

EV .= EV3r;
pt = vis3D(EV[:,1], EV[:,2], EV[:,3]; αboost3D = 8, deform = false);
display(pt);
