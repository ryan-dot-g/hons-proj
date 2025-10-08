############ -------------------------------------------------- ############
############ ---------------- PREV CODE STUFF ----------------- ############
############ -------------------------------------------------- ############

# function VolInt(X, Y, Xdash, Ydash)
#     # volume itnegral take 2
#     # normalised by π
#     fn = X .* (X.*Ydash .- Xdash.*Y)
#     return 2/3 * integrateξdom(fn);
# end

# within SurfaceFriction
    # Option 1: The difference between the COM and zero
    # center of mass in x direction is naturally 0 by symm
    # return sum(Y)^2;

    # Option 2: Pin the base 
    # return (Y[1] - Y0[1])^2

############ -------------------------------------------------- ############
############ --------------- CHECKING FORMULAS ---------------- ############
############ -------------------------------------------------- ############
checkFormulas = false
if checkFormulas 
    X = r * Cosξ; Y = r * Sinξ; 
    Xdash = FDX(X); Ydash = FDY(Y); 

    # check trace of strain tensor calc (CHECKED)
    en = trStrSq(X, Y, Xdash, Ydash); 
    # println(maximum( abs.(en .- trStrSq0) )); 

    # check first-derivative matrix (CHECKED)
    # error on bdy is about double error on interior, but still scales with ds as expected
    XdashTrue = - r * Sinξ;
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
############ ---------------- TROUBLESHOOTING FNS ------------------ ############
############ -------------------------------------------------- ############
function TroubleshootMorphUpdates()

    # wrapper to test morphogen updates on a nonuniform shape 
    X .= X0test; Y .= Y0test; Xdash .= X0testDash; Ydash .= Y0testDash
    ϕ .= ϕ0;
    pltA = visShape(ϕ, X, Y,"Before morphogen updates"); # display(pltA); savefig(pltA,"testMorph-shapeB4.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    pltB = visQtys(ϕ, ϵsq, X, Y, false); # display(pltA); savefig(pltB, "testMorph-morphB4.png")

    for i = 1:5;
        res = UpdateMorphogen!(ϕ, X, Y, ϵsq);
    end

    pltC = visShape(ϕ, X, Y, "After morphogen updates"); # display(pltC); savefig(pltC,"testMorph-shapeAF.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    pltD = visQtys(ϕ, ϵsq, X, Y, false); # display(pltD); savefig(pltD,"testMorph-morphAF.png")

    combined = plot(pltA, pltC, 
                    pltB, pltD, 
                    layout = (2,2), size = (800, 500),  margin = 6mm); 
    display(combined); savefig(combined, "testMorph.png")

end

function TroubleshootShapeUpdates()

    # code wrapper for when shape updates need to be tested/plotted 
    X .= X0; Y .= Y0; Xdash .= X0dash; Ydash .= Y0dash;
    ϕ .= ϕ0doubleBump;
    pltA = visShape(ϕ, X, Y,"Before shape update"); # display(plt); savefig(plt,"testShape-shapeB4.png")
    ϵsq = ϵ0;
    pltB = visQtys(ϕ, ϵsq, X, Y, "Before shape update"); # display(plt); savefig(plt, "testShape-morphB4.png")

    res = UpdateShape!(F.(ϕ), X, Y, Xdash, Ydash);

    pltC = visShape(ϕ, X, Y, "After shape update"); # display(plt); savefig(plt,"testShape-shapeAF.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    pltD = visQtys(ϕ, ϵsq, X, Y, "After shape update"); # display(plt); savefig(plt,"testShape-morphAF.png")

    combined = plot(pltA, pltC, 
                    pltB, pltD, 
                    layout = (2,2), size = (800, 500),  margin = 6mm); 
    display(combined); savefig(combined, "testShape.png")

end

function TroubleshootSphericalHarmonics()
    # code wrapper to get spherical harmonics from diffusion 
    X .= X0; Y .= Y0; Xdash .= X0dash; Ydash .= Y0dash;
    ϕ .= ϕ0doubleBump;

    pltA = visShape(ϕ, X, Y,"Before diffusion eigenfinding"); # display(plt); savefig(plt,"testSphHm-shapeB4.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    pltB = visQtys(ϕ, ϵsq, X, Y, "Before diffusion eigenfinding", plotE = false); # display(plt); savefig(plt, "testSphHm-morphB4.png")

    for i = 1:1000;
        res = UpdateMorphogen!(ϕ, X, Y, ϵsq);
        (dZZ, projRes, evalEst) = ProjectEvec!(ϕ, X, Y);
        dϕ = dZZ[:,1]; dX = dZZ[:,2]; dY = dZZ[:,3]; 
    end

    ϕ = ϕ * 100000; ϕ .-= minimum(ϕ); # trying to scale it sensibly
    pltC = visShape(ϕ, X, Y, "After diffusion eigenfinding"); # display(plt); savefig(plt,"testSphHm-shapeAF.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    pltD = visQtys(ϕ, ϵsq, X, Y, "After diffusion eigenfinding", plotE = false, plotTrig = true); # display(plt); savefig(plt,"testSphHm-morphAF.png")

    combined = plot(pltB, pltD, 
                    layout = (1,2), size = (800, 300),  margin = 2mm); 
    display(combined); savefig(combined, "testSHM.png")

end

function TroubleshootSPHM2()
    # better, finds all evals and evecs 
    X .= X0; Y .= Y0; Xdash .= X0dash; Ydash .= Y0dash;
    ϕ .= ϕ0bumpy;
    ϕ2 = ϕ*0 .+ ϕ0bumpy; 
    ϕ3 = ϕ*0 .+ ϕ0bumpy;
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    pltA = visQtys(ϕ, ϵsq, X, Y, "Before diffusion eigenfinding", plotE = false); # display(pltA);

    for i = 1:3000; # 3000
        rest = UpdateMorphogen!(ϕ, X, Y, ϵsq);
        (dZZt, projRes, evalEst) = ProjectEvec!(ϕ, X, Y);
        dϕ = dZZt[:,1]; dX = dZZt[:,2]; dY = dZZt[:,3]; 
        dZZtn = inp(dZZt, dZZt);

        UpdateMorphogen!(ϕ2, X, Y, ϵsq);
        (dZZt2, projRes2, evalEst2) = ProjectEvec!(ϕ2, X, Y; 
                dEVs = [dZZt], dEVnorms = [dZZtn]);
        dZZt2n = inp(dZZt2, dZZt2);

        UpdateMorphogen!(ϕ3, X, Y, ϵsq);
        (dZZt3, projRes3, evalEst3) = ProjectEvec!(ϕ3, X, Y; 
                dEVs = [dZZt2, dZZt], dEVnorms = [dZZtn, dZZt2n]);
    end
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    # plt = visQtys(ϕ, ϵsq, X, Y, "After diffusion eigenfinding", plotE = false); display(plt); 
    
    # plot evecs 
    pltB = plot(ξi, ϕ, color = 1, label = L"\varphi_1", lw = 2,
                xlabel = "ξ", ylabel = "Morphogen", title = "After diffusion eigenfinding",
                legend = :topleft, legendfontsize = 10); 
    plot!(ξi, ϕ2, color = 2, label = L"\varphi_2", lw = 2); plot!(ξi, ϕ3, color = 3, label = L"\varphi_3", lw = 2); 
    
    # next, compute hypotheses
    rescale(home, target) = minimum(target) .+ (home .- minimum(home)) .* ((maximum(target) - minimum(target)) / (maximum(home) - minimum(home)));
    e1 = 0*ϕ .+ ϕ[1];
    e2 = rescale(-Sinξ, ϕ2); e3 = rescale(cos.(2*ξi), ϕ3);

    # finally, plot hypotheses 
    plot!(ξi, e1, color = :black, ls = :dot, lw = 1, label = L"Y_{ij}");
    plot!(ξi, e2, color = :black, ls = :dot, lw = 1, label = "");
    plot!(ξi, e3, color = :black, ls = :dot, lw = 1, label = "");

    combined = plot(pltA, pltB, 
                    layout = (1,2), size = (800, 300),  margin = 2mm); 
    display(combined); savefig(combined, "testSHM.png")
end


function TroubleshootICs()
    ## computing initial states as numerical optimisers not analytic optimisers 
    # resultNumShape = optimize(x -> EFwrapped(x, ϕ0), 
    #                     [X0; Y0], 
    #                     # lbounds, ubounds, Fminbox(GradientDescent()),
    #                     method = LBFGS(), 
    #                     autodiff = :forward, iterations = 100000);  # 15000 or 100 000
    # fullXnum = resultNumShape.minimizer;
    # X0num = fullXnum[1:Ndisc]; Y0num = fullXnum[Ndisc+1:end];

    # resultNumPhi = optimize(ϕ -> MorphogenFunctional(ϕ, ϕ0, X0num, Y0num, ϵ0), 
    #                     ϕ0, method = LBFGS(), 
    #                     autodiff = :forward, iterations = 1000); 
    # ϕ0num = resultNumPhi.minimizer;

    # ZZ0 = hcat(ϕ0num, X0num, Y0num);
end