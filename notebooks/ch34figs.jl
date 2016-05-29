# plotting scripts for the book


# ===========================================================================


# homogeneous triangular lattice
X = domain(10.0, NoDefect)
B = bonds(X)
axis = [-5.3; 5.3; -3.2; +3.2]
ctx = compose( context(axis),
               compose_atoms(X, [0.15], 0.5, axis, MsUtils.blue),
               compose_bonds(X, B, 0.6, MsUtils.blue, axis) )
draw(MsUtils.auto_img(PDF, "ex-trilattice.pdf", axis, 10cm), ctx)
draw(MsUtils.auto_img(SVG, axis, 10cm), ctx)


# ===========================================================================

# point defect in infinite medium
X = domain(10.0, Interstitial)
B = bonds(X)
ee = prestrain(X, B, Interstitial)
ctx = plot_int( X, B, ee, axis=[-ξi[1]-4.0; ξi[1]+5.0; ξi[2]-3; ξi[2]+3],
filename="ex_int_atm.pdf", printwidth=10cm, plotwidth=10cm );


# ===========================================================================

# interstitial strain field
m = AtmModel(R = 20, defect=Interstitial)
u = dof2disp(m, solve(m))
Du = fdiff(u, m.B)
axis = [ξi[1]-7; ξi[1]+7; ξi[2]-5; ξi[2]+5]
MsUtils.plot_strain(m.X, m.B, log(abs(m.e+Du)+2e-4);
                    axis=axis, cmap = colormap("blues"),
                    filename = "apx-strain-int.pdf", printwidth=12cm);

# screw dislocation strain field
m = AtmModel(R = 20, defect=Screw)
u = dof2disp(m, solve(m))
Du = fdiff(u, m.B)
axis = [ξs[1]-7; ξs[1]+7; ξs[2]-5; ξs[2]+5]
MsUtils.plot_strain(m.X, m.B, log(abs(m.e+Du)+3e-2);
                   axis=axis, cmap = colormap("blues"),
                   filename = "apx-strain-screw.pdf", printwidth=12cm);
u_tot = predictor(m.X, Screw) + u
MsUtils.plot_displacement(m.X, m.B, u_tot);

# ===========================================================================

# ENVELOPES of strain fields
R = 50    # computational domain radius
Rp = 20   # radius of domain that we actually plot (to avoid boundary effects)
# --- interstitial
mi = AtmModel(R = R, defect=Interstitial, defcore = ξi_off)
ui = solve(mi)
Dui = fdiff(ui, mi.B)
ri = bonddist(mi.X, mi.B, ξi_off)
I = find(ri .<= Rp)
ri = ri[I]; Dui = Dui[I]
# --- screw
ms = AtmModel(R = R, defect=Screw, defcore = ξs_off)
us = solve(ms)
Dus = fdiff(us, ms.B)
rs = bonddist(ms.X, ms.B, ξs_off)
I = find(rs .<= Rp)
rs = rs[I]; Dus = Dus[I]; es = ms.e[I]
# ---- plot
# to make this a nice plot, it takes a little bit of work
xbins = logspace(0, log10(Rp*0.99), 15)
xi, yi = MsUtils.envelope(ri, abs(Dui), xbins)
xs, ys = MsUtils.envelope(rs, abs(Dus), xbins)
xstot, ystot = MsUtils.envelope(rs, abs(es + Dus), xbins)
p = Axis([
    Plots.Linear(xi, yi, style="thick", legendentry="interstitial");
    Plots.Linear(xs, ys, style="thick", legendentry="screw, corrector");
    Plots.Linear(xstot, ystot, style="green, thick", legendentry="screw, total");
    MsUtils.plot_slope(6.0, 20.0, 0.5, -1.0);
    Plots.Node(L"\sim r_b^{-1}", 12.0, 1e-1);
    MsUtils.plot_slope(6.0, 20.0, 0.07, -2.0);
    Plots.Node(L"\sim r_b^{-2}", 10.0, 3e-3);
        ],
        ymode="log", xmode="log",
        xlabel=L"$r_b$",
        ylabel = "strain (envelope)",
        legendPos="south west" )
# display(p)
save("ex-decay-all.pdf", p)


# plotting the A/C geometry
ac = AcModel(Ra=3.2, Rc=7, defect=Interstitial)
MsUtils.plot_acgeom(ac; filename="acgeom.pdf", printwidth=12cm, plotwidth=nothing);




# ===========================================================================
# strain path integrals

m = AtmModel(R = 20, defect=Screw)
Atri = [[1.0;0.0] [0.5; sqrt(3)/2]]
path1 = Atri * [ [-1;-1] [0;-1] [1;-1] [2;-1] [3;-1] [2;0] [1;1] [0;2] [-1;3] [-1;2] [-1;1] [-1;0] [-1;-1] ]
path2 = Atri * [ [5;-1] [6;-1] [7;-1] [8;-1] [8;0] [7;1]  [6;2]  [5;3]  [4;3] [4;2]  [4;1]  [4;0]  [5;-1] ]

function strain_path(path)
    Ipath = MsUtils.find_index(m.X, path)
    u = predictor([-m.X[1,:]; m.X[2,:]], Screw) + solve(m)
    upath = u[Ipath]
    α = grad2strain(diff(upath))
    Iα = [upath[1]; upath[1] + cumsum(α)] + 0.05
    return upath, Iα
end

Idu1, Iα1 = strain_path(path1)
Idu2, Iα2 = strain_path(path2)
xx = 1:length(Idu1)

p = Axis([
    Plots.Linear(xx, Idu1, legendentry=L"$\sum_{b \in p} Du_b$",
                style="thick, black, dashed, mark=o, mark options={solid, black}");
    Plots.Linear(xx, Iα1, legendentry=L"$\sum_{b \in p} \alpha_b$",
                style="thick, black, mark=*, mark options={solid, black, fill=black}");
    Plots.Linear(xx, Idu2, legendentry=L"$\sum_{b \in p'} Du_b$",
                style="thick, red!70!black, dashed, mark=square, mark options={solid}");
    Plots.Linear(xx, Iα2, legendentry=L"$\sum_{b \in p'} \alpha_b$",
                 style="thick, red!70!black, solid, mark=square*, mark options={fill=red!70!black}" )
        ],
        xlabel="",
        ylabel = "",
        legendPos="south west" )
display(p)
save("bvec_cumsum.svg", p)
