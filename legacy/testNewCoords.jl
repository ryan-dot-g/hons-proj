using QuadGK, Plots

R = 10;
Phi(y) = -asin(y/R);
f(y) = tan(Phi(y)); # integrand

Xint(y) = quadgk(f, 0, y)[1]; # x coord as integral 

Y = -R:0.1:R;
X = Xint.(Y) .+ R;

plt = plot(X, Y, xlabel="X", ylabel="Y")
display(plt);

