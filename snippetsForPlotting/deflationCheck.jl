using LinearAlgebra;

# setup of problem
M = [0 1; 1 0];
Id = Matrix{Float64}(I, 2, 2)
X = [1;2];
dt = 0.1;
Nsteps = 100

λest = zeros(10); lenλest = length(λest);

for i in 1:Nsteps
    global X, λ;
    X = M*X; # evolve in time 
    nx = norm(X);
    X = X/nx;
    println(nx)

    if Nsteps - i < lenλest
        λest[Nsteps-i+1] = log(nx)/dt;
    end
end
λ = sum(λest)/lenλest;
η = exp(λ*dt)

println(X, "\t", η)