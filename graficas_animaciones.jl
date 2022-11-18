using FileIO
using Plots, Colors
#using CairoMakie, Colors            # Graficas
#CairoMakie.activate!(type = "svg")
gr()
# Directorio
numCarp = "600_3-5";
direc = string(pwd(),"\\gifMagico",numCarp)

# Parametros computacionales
parametros = load( string(direc,"\\parametros.jld"),"nx","ny","nt","dx","dt","und","unt");

# Parametros para los materiales
aux =  load( string(direc,"\\casoMaterialLuneburg.jld"),"caso");

caso = aux;
if caso == "luneburg"        # Lente de Luneburg
    materiales = load( string(direc,"\\casoMaterialLuneburg.jld"),"rad","centr")
    centr = materiales[2];
    th = collect(0:pi/32:2*pi);
    grafObjetos = plot()
    plot!(grafObjetos,parametros[4]*materiales[1]*cos.(th).+centr[2]*parametros[4], materiales[1]*parametros[4]*sin.(th).+centr[1]*parametros[4],label=false,linecolor=:black ,linewidth=1)
    plot(grafObjetos)
    
elseif caso == "luneburgColumna"    # Columna de Luneburgs
    materiales = load( string(direc,"\\casoMaterialLuneburg.jld"),"rad","centr")
    centr = materiales[2]';
    nL = size(centr);
    th = collect(0:pi/32:2*pi);
    grafObjetos = plot()
    for itn = 1:nL[1]
        plot!(grafObjetos,parametros[4]*materiales[1]*cos.(th).+centr[itn,2]*parametros[4], materiales[1]*parametros[4]*sin.(th).+centr[itn,1]*parametros[4],label=false,linecolor=:black ,linewidth=1)
    end
    plot(grafObjetos)

end

# Vista previa de objeto
fsx = (1/parametros[4])*parametros[6];
fsy = (1/parametros[4])*parametros[6];
domx = fsx*collect(1:parametros[1]);
domy = fsy.*collect(1:parametros[2]);

# Intensidad
inte = zeros(parametros[2],parametros[1]);
for itint = 1:1730
        info = load( string(direc,"\\matHz-",parametros[3]-itint,".jld"), "hz");
        info = info./maximum(info);
        inte = inte .+ abs.(info).^2;
end

grafZoom2 = plot(layout=grid(1,1),colorbar=:top);
heatmap!(grafZoom2[1],domx,domy,(1/1740).*inte,colorbar=:top,xlabel="x [mts]",ylabel="y [mts]",
        title="Intensidad",aspect_ratio=:equal,
        xlims=(limX[1],limX[2]),ylims=(limY[1],limY[2]),
        c = :berlin
        );
plot!(grafZoom2[1],fsx*materiales[1]*cos.(th).+centr[2]*fsx, materiales[1]*fsx*sin.(th).+centr[1]*fsx,label=false,
        linecolor=:black ,linewidth=1
        );
savefig(grafZoom2,string("intenHz",numCarp,".pdf"))    


# Skectch final
ref = 1;
itn = parametros[3];#parametros[3];
info = load( string(direc,"\\matHz-",itn,".jld"), "hz");
info = info./maximum(info);
grafFinal = heatmap(domx,domy,info,clims=(-ref,ref),xlabel="x [mts]",ylabel="y [mts]",
                    aspect_ratio=:equal,
                    xlims=(0,fsx*parametros[1]),ylims=(0,fsy*parametros[2]),
                    c = :berlin
                    )
plot!(grafFinal,fsx*materiales[1]*cos.(th).+centr[2]*fsx, materiales[1]*fsx*sin.(th).+centr[1]*fsx,label=false,
                    linecolor=:white ,linewidth=0.5
                    )
annotate!(grafFinal,[(0.026,0.014,(string("Iteración: ",itn),10,:white))])
annotate!(grafFinal,[(0.026,0.001,(string("t= ",round(itn*parametros[5]*parametros[7],digits=16)," s"),10,:white))])

# Zoom en la frontera de los lentes
ref = 1;
itn = parametros[3];#parametros[3];
info = load( string(direc,"\\matHz-",itn,".jld"), "hz");
info = info./maximum(info);
limY = fsx*round.( [centr[1]-1.5materiales[1],centr[1]+1.5materiales[1]],RoundUp);
limX = fsx*round.( [centr[2]-1.5materiales[1],centr[2]+1.5materiales[1]],RoundUp);

grafZoom1 = plot(layout=grid(1,1),colorbar=:top,clims=(-ref,ref))
heatmap!(grafZoom1[1],domx,domy,info,colorbar=:top,xlabel="x [mts]",ylabel="y [mts]",
        title="Campo Hz",aspect_ratio=:equal,
        xlims=(limX[1],limX[2]),ylims=(limY[1],limY[2]),
        c = :berlin
        );
plot!(grafZoom1[1],fsx*materiales[1]*cos.(th).+centr[2]*fsx, materiales[1]*fsx*sin.(th).+centr[1]*fsx,label=false,
        linecolor=:white ,linewidth=1
        );
savefig(grafZoom1,string("campoHz",numCarp,".pdf"))

# Cortes Transversales
ref = 1;
itn = parametros[3];#parametros[3];
refT = round(itn*parametros[5]*parametros[7],digits=16);
info = load( string(direc,"\\matHz-",itn,".jld"), "hz");
info = info./maximum(info);
inte = abs.(info).^2;

grafCT = plot(layout=grid(2,1))
plot!(grafCT[1],domx,info[div(parametros[2],2),:],aspect_ratio=:auto,xlims=(0,parametros[1]*fsx),ylims=(-1,1),
        title="Amplitud Campo Hz",xlabel="x [mts]",ylabel="Amplitud",label = string("t=",refT," s"),
        linecolor=:purple)
plot!(grafCT[2],domx,inte[div(parametros[2],2),:],aspect_ratio=:auto,xlims=(0,parametros[1]*fsx),ylims=(0,1),
        title="Intensidad Campo Hz",xlabel="x [mts]",ylabel="Intensidad",label = string("t=",refT," s"),
        linecolor=:purple)

savefig(grafCT,string("campoHz",numCarp,"_CT.pdf"))

# Creación del gif
ref = 1;
@time begin
anim = @animate for itn = 1:4:parametros[3]
    info = load( string(direc,"\\matHz-",itn,".jld"), "hz");
    info = info./maximum(info);
    grafFinal = heatmap(domx,domy,info,clims=(-ref,ref),xlabel="x [mts]",ylabel="y [mts]",
                        aspect_ratio=:equal,
                        xlims=(0,fsx*parametros[1]),ylims=(0,fsy*parametros[2]),
                        c = :berlin
                        )
    plot!(grafFinal,fsx*materiales[1]*cos.(th).+centr[2]*fsx, materiales[1]*fsx*sin.(th).+centr[1]*fsx,label=false,
                        linecolor=:white ,linewidth=0.5
                        )
    annotate!(grafFinal,[(0.026,0.014,(string("Iteración: ",itn),10,:white))])
    annotate!(grafFinal,[(0.026,0.001,(string("t= ",round(itn*parametros[5]*parametros[7],digits=16)," s"),10,:white))])
end
end

gif(anim,string("anim",numCarp,".gif"), fps = 120)
