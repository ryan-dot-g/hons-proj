using LinearAlgebra, Plots, Optim;
gitFP = "C:/Users/rgray/OneDrive/ryan/Uni/HONOURS/hons-proj";
suppressFP = gitFP * "/suppress.txt";

Ndisc = 50; # number of discretisation points on s (50)
smin = -π/2; smax = π/2; # bounds of s values (-π/2, π/2)
R = 80; # radius of force-free reference sphere (80)
r = 120; # TODO radius of initial condition sphere (120 for now) 
Ω = 1e2; # large number: punishing potential for volume deviation (1e4)
ω = 0.0; # small number: surface friction 
dsint = 0.01; # small number: distance inside the s grid to start at to avoid div0
Si = range(smin+dsint, smax-dsint, Ndisc); # grid of s values 
ds = Si[2] - Si[1]; 
interior = 2:Ndisc-1; # index set for the interior of the array

# integral over total s domain, using trapezoid rule
volEl = R^2 * cos.(Si);
integrateSdom(Qty) = ds * ( sum(Qty) - (Qty[1] + Qty[end])/2 );

ϵ0sc = 0.5 * ((r^2 - R^2)/R^2)^2; # scalar initial (constant) trace of strain tensor squared
ϵ0 = zeros(Ndisc) .+ ϵ0sc; # over s-grid
V0 = 4/3 * r^3; # initial volume, normalised by pi
ϕ = zeros(Ndisc) .+ 2;

X0 = r * cos.(Si); Y0 = r * sin.(Si); # initial condition shape 
X0dash = -Y0[:]; Y0dash = X0[:]; # initial derivatives 
Xundef = R * cos.(Si); Yundef = R * sin.(Si); # undeformed shape

function trStrSq(X, Y, Xdash, Ydash)
    # takes vectors of the surface parameterisation x(s), y(s) and derivatives w.r.t. s
    # returns a vector containing the trace of the square of the strain tensor at each s
    t1 = (Xdash).^2 .+ (Ydash).^2 .- R^2;
    
    t2 = X.^2 ./ cos.(Si).^2 .- R^2;

    return 1/(4*R^4) * (t1.^2 .+ t2.^2); 
end

function ElEn(ϕ, X, Y, Xdash, Ydash)
    # Elastic energy functional 
    fn = trStrSq(X, Y, Xdash, Ydash) .* volEl; # integrand over S0 including volel
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
    # additionally uses surface friction 
    return ElEn(ϕ, X, Y, Xdash, Ydash) + Ω * (VolInt(X, Y, Xdash, Ydash) - V0)^2 +
                ω * SurfaceFriction(X, Y);
end

function eng(fullX)
    # wrapper to optimise over just 1 x
    X = fullX[1:Ndisc]; Y = fullX[Ndisc+1:end];
    Xdash .= P*X; Ydash .= P*Y;
    return EnergyFunctional(ϕ, X, Y, Xdash, Ydash, Ω);
end

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

α = 0.2;
X0test2 = α*X0test + (1-α)*X0; Y0test2 = α*Y0test + (1-α) * Y0;

result = optimize(eng, [X0test2; Y0test2], LBFGS());
fullX = result.minimizer;
X = fullX[1:Ndisc]; Y = fullX[Ndisc+1:end];

ϕ0test = zeros(Ndisc) .+ sin.(Si) * 3 .+ 3;


plt = visualise(ϕ, X, Y); display(plt)

#### TO DO TMRW
# - move the X0test2 test case to the other code, its a better test case
# - port the tetsing code over, remove Popt, upload to git BEFORE AND AFTER
# - test with non-uniform phi.
# - make sure volume constraint stays small!!

