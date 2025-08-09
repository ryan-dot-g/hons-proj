############ -------------------------------------------------- ############
############ --------------- CHECKING FORMULAS ---------------- ############
############ -------------------------------------------------- ############
checkFormulas = false
if checkFormulas 
    X = r * CosS; Y = r * SinS; 
    Xdash = FDX(X); Ydash = FDY(Y); 

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
############ ---------------- TROUBLESHOOTING FNS ------------------ ############
############ -------------------------------------------------- ############
function TroubleshootMorphUpdates()
    # wrapper to test morphogen updates on a nonuniform shape 
    X .= X0test; Y .= Y0test; Xdash .= X0testDash; Ydash .= Y0testDash
    ϕ .= ϕ0;
    plt = visShape(ϕ, X, Y,"Before morphogen updates"); display(plt); savefig(plt,"testMorph-shapeB4.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    plt = visQtys(ϕ, ϵsq, "Before morphogen updates"); display(plt); savefig(plt, "testMorph-morphB4.png")

    # ϵsqtot = ElEn(ϕ, X, Y, Xdash, Ydash)
    # println("Elastic energy: ", ϵsqtot)
    # println("Total volume energy: ", Ω * (VolInt(X, Y, Xdash, Ydash) - V0)^2)

    for i = 1:100;
        res = UpdateMorphogen!(ϕ, X, Y, ϵsq);
    end
    plt = visShape(ϕ, X, Y, "After morphogen updates"); display(plt); savefig(plt,"testMorph-shapeAF.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    plt = visQtys(ϕ, ϵsq, "After morphogen updates"); display(plt); savefig(plt,"testMorph-morphAF.png")
end

function TroubleshootShapeUpdates()
    # code wrapper for when shape updates need to be tested/plotted 
    X .= X0; Y .= Y0; Xdash .= X0dash; Ydash .= Y0dash;
    ϕ .= ϕ0doubleBump;
    plt = visShape(ϕ, X, Y,"Before shape update"); display(plt); savefig(plt,"testShape-shapeB4.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    plt = visQtys(ϕ, ϵsq, "Before shape update"); display(plt); savefig(plt, "testShape-morphB4.png")

    # ϵsqtot = ElEn(ϕ, X, Y, Xdash, Ydash)
    # println("Elastic energy: ", ϵsqtot)
    # println("Total volume energy: ", Ω * (VolInt(X, Y, Xdash, Ydash) - V0)^2)

    res = UpdateShape!(ϕ, X, Y, Xdash, Ydash)
    plt = visShape(ϕ, X, Y, "After shape update"); display(plt); savefig(plt,"testShape-shapeAF.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    plt = visQtys(ϕ, ϵsq, "After shape update"); display(plt); savefig(plt,"testShape-morphAF.png")
end

function TroubleshootSphericalHarmonics()
    # code wrapper to get spherical harmonics from diffusion 
    X .= X0; Y .= Y0; Xdash .= X0dash; Ydash .= Y0dash;
    ϕ .= ϕ0doubleBump;
    plt = visShape(ϕ, X, Y,"Before diffusion eigenfinding"); display(plt); savefig(plt,"testSphHm-shapeB4.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    plt = visQtys(ϕ, ϵsq, "Before diffusion eigenfinding", plotE = false); display(plt); savefig(plt, "testSphHm-morphB4.png")

    for i = 1:1000;
        res = UpdateMorphogen!(ϕ, X, Y, ϵsq);
        (dZZ, projRes) = ProjectEvec!(ϕ, X, Y);
        dϕ = dZZ[:,1]; dX = dZZ[:,2]; dY = dZZ[:,3]; 
    end
    plt = visShape(ϕ, X, Y, "After diffusion eigenfinding"); display(plt); savefig(plt,"testSphHm-shapeAF.png")
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    plt = visQtys(ϕ, ϵsq, "After diffusion eigenfinding", plotE = false, plotTrig = true); display(plt); savefig(plt,"testSphHm-morphAF.png")
end

function TroubleshootSPHM2()
    # better, finds all evals and evecs 
    X .= X0; Y .= Y0; Xdash .= X0dash; Ydash .= Y0dash;
    ϕ .= ϕ0bumpy;
    ϕ2 = ϕ*0 .+ ϕ0bumpy; 
    ϕ3 = ϕ*0 .+ ϕ0bumpy;
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    plt = visQtys(ϕ, ϵsq, "Before diffusion eigenfinding", plotE = false); display(plt);

    for i = 1:3000;
        rest = UpdateMorphogen!(ϕ, X, Y, ϵsq);
        (dZZt, projRes) = ProjectEvec!(ϕ, X, Y);
        dϕ = dZZt[:,1]; dX = dZZt[:,2]; dY = dZZt[:,3]; 

        UpdateMorphogen!(ϕ2, X, Y, ϵsq);
        (dZZt2, projRes2) = ProjectEvec!(ϕ2, X, Y; EVs = [dZZt]);

        UpdateMorphogen!(ϕ3, X, Y, ϵsq);
        (dZZt3, projRes3) = ProjectEvec!(ϕ3, X, Y; EVs = [dZZt2, dZZt]);
    end
    ϵsq = trStrSq(X, Y, Xdash, Ydash); 
    plt = visQtys(ϕ, ϵsq, "After diffusion eigenfinding", plotE = false); display(plt); 

    plt = plot(Si, ϕ); plot!(Si, ϕ2); plot!(Si, ϕ3); display(plt);
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