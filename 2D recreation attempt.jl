using LinearAlgebra, DifferentialEquations, Plots, LaTeXStrings, BandedMatrices;

### ----------- HELPER FUNCTIONS ------------- ###
function randbetween(a, b, size)
    # generates an array of specified size, of random unif numbers
    # between a and b 
    r = rand(size);
    r = r * (b-a); # adjust width
    r = r .+ a; # add bottom amount 
    return r 
end

### ----------- PHYSICAL PARAMETERS ------------- ###
zeta = 1.0; # dissipation constant
D = 0.1; # diffusion constant
f(e) = 1.3 * e; # morphogen response to strain
E(phi) = exp.(-1.2 * phi); # strain response to morphogen 
Xi = 1e-5; # small amount of friction
L = 1; l = 1.2; # arclengths

### ----------- SIMULATION PARAMETERS ------------- ###
Npoints = 100; # number of points to discretise terrain into
dt = 0.1; # timestep
tmax = 100; # maximum time

### ----------- SIMULATION SETUP ------------- ###
I_n = Matrix(I, Npoints, Npoints); # explicitly construct identity 
dx = 1/Npoints; # grid spacing
# construct diffusion matrix with periodic BCs
M = Matrix(BandedMatrix(Ones(Npoints,Npoints),(1,1)) - 3I); 
M[1, Npoints] = 1; M[Npoints, 1] = 1;
M = M/dx^2;
# construct first-derivative matrix with periodic BCs
Q = Matrix(Tridiagonal(zeros(Npoints-1) .-1, zeros(Npoints), zeros(Npoints-1) .+ 1));
Q[1, Npoints] = -1; Q[Npoints, 1] = 1;
Q = Q/(2*dx);
# construct helper matrix
K = ((1/dt) + zeta)*I_n - D*M;
Kinv = inv(K); # presave inverse for performance 
S0 = (0:dx:L-dx) * 2*pi; # initial arclengths 

### ----------- PLOTTING PARAMETERS ------------- ###
plotres = 1e100; # how many timesteps per plot 
pointstot = 10; # how many total points to show
pointres = Npoints ÷ pointstot; # ratio to index list with 
sleepytime = 0.1; # how much to pause in between each plot 
xlim = (-1.2, 1.2); ylim = (-1.2, 1.2); # constant plot limits

function main()
    ### ----------- ICs ------------- ###
    Eps = zeros(Npoints) .+ (l/L -1);
    Phi = 1/zeta .* f(Eps); 
    Phi = Phi .* randbetween(0.9, 1.1, Npoints); # perturb by 10%
    ### ----------- SIMULATION RUN ------------- ###
    for t in 0:dt:tmax
        P = Xi*I_n .- dt*M*diagm(E(Phi)); # construct helper matrix 
        Eps = inv(P) * (Xi*I_n) * Eps; # evolve Eps in time 
        # Eps = inv(I_n .- dt/Xi*Q*diagm(E(Phi))) * Eps # first derivative option
        Phi = Kinv * (1/dt * Phi .+ f(Eps)); # evolve Phi in time 

        # plot some iterations and always last iter
        if mod(t÷dt, plotres) == 0 || t > tmax-dt || t == 0;
            U = cumsum(Eps) * dx; # displacements, sum replaced by integral
            S = S0 .+ U; # positions
            Xplt = cos.(S); Yplt = sin.(S); # x, y parameterised by arclength (= angle)
            plt = scatter(Xplt[1:pointres:end], Yplt[1:pointres:end], 
                    xlim = xlim, ylim = ylim, 
                    show = true, legend = false,
                    aspect_ratio = :equal,
                    xlabel = "x", ylabel = "y",
                    title = "t = "*string(round(t,digits=3)));
            display(plt)
            sleep(sleepytime);
        end
    end

    # show final result 
    println(Phi);
    plt = plot(S0, Phi, label = L"\varphi", 
                xlabel = L"\xi");
    plot!(S0, Eps, label = L"\varepsilon");
    display(plt)
end;


