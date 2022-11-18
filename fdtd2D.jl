# FDTD en dos dimensiones

#     Program author: Susan C. Hagness
#                     Department of Electrical and Computer Engineering
#                     University of Wisconsin-Madison
#                     1415 Engineering Drive
#                     Madison, WI 53706-1691
#                     hagness@engr.wisc.edu

# Se traduce de lenguaje MatLab a lenguaje Julia
# Modifiaciones para simplificar un poco el código y redefinir cosntantes.

#######################################################################
# Paquetes
#######################################################################

using LinearAlgebra                 # Declaración de variables
using Printf                        # Impresión de datos
using JLD, FileIO                   # Para guardar datos
#######################################################################
# Constantes físicas
#######################################################################

println("Declaración de constantes físicas")

# Redifinición de unidades de longitud y tiempo
und = 1*10^-6;                       # Unidad espacial
unt = 1*10^0;                       # Unidad temporal

epsz = (8.854187817*10^-12)*((und^3)/(unt^4)); 
muz = (4*pi*1e-7)*((unt^2)/und);
etaz=sqrt(muz/epsz);
cc =  sqrt(1/(epsz*muz));

freq=3.3e+14*(unt);                        # Frecuencia de la fuente de campo
lambda=cc/freq;                     # Longitud de onda de la fuente de campo 
omega=2.0*pi*freq;                  # Frecuencia en radiandes

#######################################################################
# Parámetros de la malla numérica espacial
#######################################################################

println("Declaración de parámetros computacionales")

# Diferenciales
CFL = 99e-2;                                  # Constante de Courant
dx = lambda/20;#(100e-3)*und;                             # Diferencial espacial
dt = CFL/(cc*sqrt((1/dx^2)+(1/dx^2)));      # Diferencial temporal

tau = Int(round(1/freq/dt));                # Nodos por periodo
delay = 3*tau;
scl = 60;                                   # Cantidad de periodos

nmax= length(1:scl*tau);         # Nodos temporales

# Ventana numérica
ventNum = round.([30*10^-3 15*10^-3]/und*dx,RoundUp)

ie= Int(ventNum[2]);                                     # nodos en y
je= Int(ventNum[1]);                                     # nodos en x

# Impresión de parámetros
@printf("Velocidad de la luz: %.3e \n",cc);              @printf("Longitud de onda: %.3e \n",lambda*und);
@printf("Unidad de distancia: %.3e \n",und);            @printf("Unidad de tiempo: %.3e \n",unt)
@printf("Diferencial de distancia: %.3e \n",dx);    @printf("Diferencial de tiempo: %.3e \n",dt)
@printf("Cantidad de nodos temporales: %i \n",nmax);      @printf("Cantidad de nodos en x: %i \n",ie);    @printf("Cantidad de nodos en y: %i \n",je)

ib=ie+1;
jb=je+1;

# Ubicación de la fuente del campo
is=div(ie,2)-15:div(ie,2)+15; # ubicacion en y                                 
js=div(je,6);                               # ubicación en x

# Parámetros para PMLs
iebc=8;                                 # espesor de PMLs derecha e izquierda
jebc=8;                                 # espesor de PMLs adelante y atras 
rmax=1.0e-7;                            # R[0]
orderbc=4;                              # Orden de la función de sigmas
ibbc=iebc+1;
jbbc=jebc+1;
iefbc=ie+2*iebc;
jefbc=je+2*jebc;
ibfbc=iefbc+1;
jbfbc=jefbc+1;

########################################################################
# Fuente de campo eléctrico
#######################################################################
#println("Declaración de la fuente")

#source = sin.(omega*(collect(1:scl*tau).-delay)*dt).*exp.(-((collect(1:scl*tau).-delay).^2/tau^2));

# Guardar parámetros espaciales
save(File(format"JLD",string(pwd(),"\\gifMagico\\parametros.jld")),"nx",je,"ny",ie,"nt",nmax,"dx",dx,"dt",dt,"und",und,"unt",unt)

########################################################################
# Inicialización de matrices
#######################################################################

println("Declaración de matrices")
# Campos en la ventana principal
  ex =    zeros(ie,jb);                    
  ey =    zeros(ib,je);
  hz =    zeros(ie,je);

# Campos en PML Frontal
  exbcf = zeros(iefbc,jebc);            
  eybcf = zeros(ibfbc,jebc);
  hzxbcf= zeros(iefbc,jebc);
  hzybcf= zeros(iefbc,jebc);

# Campos en PML Atrás
  exbcb = zeros(iefbc,jbbc);            
  eybcb = zeros(ibfbc,jebc);
  hzxbcb = zeros(iefbc,jebc);
  hzybcb = zeros(iefbc,jebc);

# Campos en PML Izquierda
  exbcl = zeros(iebc,jb);               
  eybcl = zeros(iebc,je);
  hzxbcl= zeros(iebc,je);
  hzybcl= zeros(iebc,je);

# Campos en PML Derecha
  exbcr = zeros(iebc,jb);               
  eybcr = zeros(ibbc,je);
  hzxbcr= zeros(iebc,je);
  hzybcr= zeros(iebc,je);

#

# Las constantes si inicializan en el vacio
# aux = dt*sig[1]/(2.0*epsz*eps[1]);
# caE = (1-aux)/(1+aux);
# cbE = dt/epsz/eps[1]/dx/(1+aux);

# aux = dt*sim[1]/(2.0*muz*mur[1]);
# caH = (1-aux)/(1+aux);
# cbH = dt/muz/mur[1]/dx/(1+aux);

caE = 1;      cbE = dt/epsz/dx;
caH = 1;      cbH = dt/muz/dx;

# Constantes para la venta numérica
  caex=       ones(ie,jb);                  
  cbex=       cbE.*ones(ie,jb);
  caey=       ones(ib,je);
  cbey=       cbE.*ones(ib,je);
  dahz=       ones(ie,je);
  dbhz=       cbH.*ones(ie,je);

# Constantes para PML frontal
  caexbcf=    ones(iefbc,jebc);          
  cbexbcf=    cbE.*ones(iefbc,jebc);
  caexbcl=    ones(iebc,jb);
  cbexbcl=    cbE.*ones(iebc,jb);
  caexbcr=    ones(iebc,jb);
  cbexbcr=    cbE.*ones(iebc,jb);
  caeybcf=    ones(ibfbc,jebc);
  cbeybcf=    cbE.*ones(ibfbc,jebc);
  dahzybcf=   ones(iefbc,jebc);
  dbhzybcf=   cbH.*ones(iefbc,jebc);
  dahzxbcf=   ones(iefbc,jebc);
  dbhzxbcf=   cbH.*ones(iefbc,jebc);

# Constantes para PML atrás
  caexbcb=    ones(iefbc,jbbc);          
  cbexbcb=    cbE.*ones(iefbc,jbbc);
  dahzybcb=   ones(iefbc,jebc);
  dbhzybcb=   cbH.*ones(iefbc,jebc);
  caeybcb=    ones(ibfbc,jebc);
  cbeybcb=    cbE.*ones(ibfbc,jebc);
  dahzxbcb=   ones(iefbc,jebc);
  dbhzxbcb=   cbH.*ones(iefbc,jebc);

# Constantes para PML izquierda
  caeybcl=    ones(iebc,je);             
  cbeybcl=    cbE.*ones(iebc,je);
  dahzxbcl=   ones(iebc,je);
  dbhzxbcl=   cbH.*ones(iebc,je);
  dahzybcl=   ones(iebc,je);
  dbhzybcl=   cbH.*ones(iebc,je);

# Constantes para PML derecha
  caeybcr=    ones(ibbc,je);             
  cbeybcr=    cbE.*ones(ibbc,je);
  dahzxbcr=   ones(ibbc,je);
  dbhzxbcr=   cbH.*ones(ibbc,je);
  dahzybcr=   ones(ibbc,je);
  dbhzybcr=   cbH.*ones(ibbc,je);

#

########################################################################
# Se define la geometría en la ventana numérica principal
#######################################################################

# Se agregan objetos en el medio
casoMaterial = "vacio";
if casoMaterial=="vacio"
  idx = 1:je;
  idy = 1:ie;
  medioProp = epsz;
  save(File(format"JLD",string(pwd(),"\\gifMagico\\casoMaterialLuneburg.jld")),"caso",casoMaterial)
elseif casoMaterial=="luneburg"
  # Lente de Luneburg
    f = 1;
    # En nodos computacionales
    centr = [div(ie,2),div(je,2)];
    rad = Int(round(( 1*lambda )/dx,RoundUp));

    idx = centr[2]-rad:centr[2]+rad;
    idy = centr[1]-rad:centr[1]+rad;

    domx = dx*collect(idx);
    domy = dx*collect(idy);
    aux = length(domx);

    domR = sqrt.( ((ones(aux,1).*domx').-centr[2]*dx).^2 + ((ones(1,aux).*domy).-centr[1]*dx).^2 );
    n = real.((1/f)* sqrt.(Complex.( (1).+(f^2).-((domR./rad*dx).^2) )));
    n[domR.>rad*dx].=1;
    medioProp = epsz*(n.^2);

    save(File(format"JLD",string(pwd(),"\\gifMagico\\casoMaterialLuneburg.jld")),"caso",casoMaterial,"f",f,"rad",rad,"centr",centr)

elseif casoMaterial=="luneburgColumna"
    # Columan de lentes de Luneburg
    nL = 5;
    rad = 60; #div(ie-101,nL);
    f = 1;
    #f = 1.5*rad;
  
    auxCr = (2*rad.*[collect(-div(nL-1,2):1:div(nL-1,2)) zeros(nL,1)])
    ctrA = [div(ie,2) div(je,4)];
    centr = ones(nL,1).*ctrA;
    centr[:,1] = centr[:,1].+auxCr[:,1];
    centr[:,2] = centr[:,2].+auxCr[:,2];
    centr = Int.(centr);
  
    n = zeros(ie,je);
    for itl = 1:nL
      idx = centr[itl,2]-rad-1:centr[itl,2]+rad-1;
      idy = centr[itl,1]-rad-1:centr[itl,1]+rad-1;
      
      domx = dx*collect(idx);
      domy = dx*collect(idy);
      aux = length(domx);
      
      domR = sqrt.( ((ones(aux,1).*domx').-centr[itl,2]*dx).^2 + ((ones(1,aux).*domy).-centr[itl,1]*dx).^2 );
      nAux = real.((1/f)* sqrt.(Complex.( (1).+(f^2).-((domR./rad*dx).^2) )));
      #alf = 0.74;
      #nAux = real.((1/f)* sqrt.(Complex.( (rad^2).+(f^2).-alf*(domR.^2) )));
      n[idy,idx] = n[idy,idx].+nAux.*(domR.<rad*dx);
    end
    n[n.==0].=1;
    idx = 1:je;
    idy = 1:ie;
    medioProp = epsz*(n.^2);

    save(File(format"JLD",string(pwd(),"\\gifMagico\\casoMaterialLuneburg.jld")),"caso",casoMaterial,"f",f,"rad",rad,"centr",centr)

end

# Cambio en las constantes de propagación
  auxObj = 0; #dt*sig[1]./(2.0*medioProp);
  caEObj = ((1).-auxObj)./((1).+auxObj);
  cbEObj = dt./medioProp./dx./((1).+auxObj);

  auxObj = 0; #dt*sim[1]/(2.0*muz*mur[1]);
  caHObj = ((1).-auxObj)/((1).+auxObj);
  cbHObj = dt/muz/1/dx/((1).+auxObj);

# Constantes para la ventana principal
  caex[idy,idx].=caEObj;                 
  cbex[idy,idx].=cbEObj;
  caey[idy,idx].=caEObj;
  cbey[idy,idx].=cbEObj;
  dahz[idy,idx].=caHObj;
  dbhz[idy,idx].=cbHObj;
#

########################################################################
# Se inicializan las PMLs 
#######################################################################

# Definición de las densidades de corriente
delbc=iebc*dx;
sigmam=-log(rmax)*(orderbc+1)/(2*etaz*delbc);   # Eq. 7.62 in text
bcfactor=sigmam/(dx*(orderbc+1)*(delbc^orderbc));  

#
# Enfrente
  # Campo eléctrico
    idx = 2:jebc;         idy = 1:iefbc;
    caexbcf[idy,1].=1.0;
    cbexbcf[idy,1].=0.0;

    dAux1 = (dx.*(jebc.-collect(idx).+1.5)).^(orderbc+1);
    dAux2 = (dx.*(jebc.-collect(idx).+0.5)).^(orderbc+1);
    sigmaPML = bcfactor.*(dAux1.-dAux2);
    aux1PML = exp.((-1).*sigmaPML.*(dt/epsz));
    aux2PML = ((1).-aux1PML)./(sigmaPML.*dx);
    caexbcf[idy,idx].= ones(iefbc,1).*aux1PML';
    cbexbcf[idy,idx].= ones(iefbc,1).*aux2PML';

    sigmaAux = bcfactor*(0.5*dx)^(orderbc+1);
    ca1=exp(-sigmaAux*dt/epsz);   cb1=(1-ca1)/(sigmaAux*dx);
    caex[1:ie,1].=ca1;
    cbex[1:ie,1].=cb1;
    caexbcl[1:iebc,1].=ca1;
    cbexbcl[1:iebc,1].=cb1;
    caexbcr[1:iebc,1].=ca1;
    cbexbcr[1:iebc,1].=cb1;

  # Campo magnético
    idx = 1:jebc;         idy = 1:iefbc;
    dAux1 = (dx.*(jebc.-collect(idx).+1)).^(orderbc+1);
    dAux2 = (dx.*(jebc.-collect(idx))).^(orderbc+1);
    sigmaPML = bcfactor.*(muz/epsz).*(dAux1.-dAux2);
    aux1PML = exp.((-1).*sigmaPML.*(dt/muz));
    aux2PML = ((1).-aux1PML)./(sigmaPML.*dx);
    dahzybcf[idy,idx].= ones(iefbc,1).*aux1PML';
    dbhzybcf[idy,idx].= ones(iefbc,1).*aux2PML';

# Atras
  # Campo eléctrico
  idx = 2:jebc;         idy = 1:iefbc;
  caexbcb[idy,jbbc].=1.0;
  cbexbcb[idy,jbbc].=0.0;

  dAux1 = (dx.*(collect(idx).-0.5)).^(orderbc+1);
  dAux2 = (dx.*(collect(idx).-1.5)).^(orderbc+1);
  sigmaPML = bcfactor.*(dAux1.-dAux2);
  aux1PML = exp.((-1).*sigmaPML.*(dt/epsz));
  aux2PML = ((1).-aux1PML)./(sigmaPML.*dx);
  caexbcb[idy,idx].= (ones(iefbc,1).*aux1PML');
  cbexbcb[idy,idx].= (ones(iefbc,1).*aux2PML');
  caexbcb[idy,1].=0;    
  cbexbcb[idy,1].=0;

  sigmaAux = bcfactor*(0.5*dx)^(orderbc+1);
  ca1=exp(-sigmaAux*dt/epsz);   cb1=(1-ca1)/(sigmaAux*dx);
  caex[1:ie,jb].=ca1;
  cbex[1:ie,jb].=cb1;
  caexbcl[1:iebc,jb].=ca1;
  cbexbcl[1:iebc,jb].=cb1;
  caexbcr[1:iebc,jb].=ca1;
  cbexbcr[1:iebc,jb].=cb1;

  # Campo magnético
  idx = 1:jebc;         idy = 1:iefbc;
  dAux1 = (dx.*(collect(idx))).^(orderbc+1);
  dAux2 = (dx.*(collect(idx).-1)).^(orderbc+1);
  sigmaPML = bcfactor.*(muz/epsz).*(dAux1.-dAux2);
  aux1PML = exp.((-1).*sigmaPML.*(dt/muz));
  aux2PML = ((1).-aux1PML)./(sigmaPML.*dx);
  dahzybcb[idy,idx].= ones(iefbc,1).*aux1PML';
  dbhzybcb[idy,idx].= ones(iefbc,1).*aux2PML';

# Izquierda
  # Campo eléctrico
  idx = 1:jebc;         idy = 2:iebc;
  caexbcl[1,1:je].=1.0;
  cbexbcl[1,1:je].=0.0;

  dAux1 = (dx.*(iebc.-collect(idy).+1.5)).^(orderbc+1);
  dAux2 = (dx.*(iebc.-collect(idy).+0.5)).^(orderbc+1);
  sigmaPML = bcfactor.*(dAux1.-dAux2);
  aux1PML = exp.((-1).*sigmaPML.*(dt/epsz));
  aux2PML = ((1).-aux1PML)./(sigmaPML.*dx);
  caeybcl[idy,1:je].= ones(1,je).*aux1PML;
  cbeybcl[idy,1:je].= ones(1,je).*aux2PML;
  caeybcf[idy,idx].= ones(1,jebc).*aux1PML;
  cbeybcf[idy,idx].= ones(1,jebc).*aux2PML;
  caeybcb[idy,idx].= ones(1,jebc).*aux1PML;
  cbeybcb[idy,idx].= ones(1,jebc).*aux2PML;

  sigmaAux = bcfactor*(0.5*dx)^(orderbc+1);
  ca1=exp(-sigmaAux*dt/epsz);   cb1=(1-ca1)/(sigmaAux*dx);
  caey[1,1:je].=ca1;
  cbey[1,1:je].=cb1;
  caexbcf[iebc+1,idx].=ca1;
  cbeybcf[iebc+1,idx].=cb1;
  caeybcb[iebc+1,idx].=ca1;
  cbeybcb[iebc+1,idx].=cb1;

  # Campo magnético
      idx = 1:jebc;         idy = 1:iebc;
    dAux1 = (dx.*(iebc.-collect(idx).+1)).^(orderbc+1);
    dAux2 = (dx.*(iebc.-collect(idx))).^(orderbc+1);
    sigmaPML = bcfactor.*(muz/epsz).*(dAux1.-dAux2);
    aux1PML = exp.((-1).*sigmaPML.*(dt/muz));
    aux2PML = ((1).-aux1PML)./(sigmaPML.*dx);
    dahzxbcl[idy,1:je].= ones(1,je).*aux1PML;
    dbhzxbcl[idy,1:je].= ones(1,je).*aux2PML;
    dahzxbcf[idy,idx].= ones(1,jebc).*aux1PML;
    dbhzxbcf[idy,idx].= ones(1,jebc).*aux2PML;
    dahzxbcb[idy,idx].= ones(1,jebc).*aux1PML;
    dbhzxbcb[idy,idx].= ones(1,jebc).*aux2PML;

# Derecha
  # Campo eléctrico
  idx = 1:jebc;         idy = 2+iebc+ie:iebc+iebc+ie;
  caeybcr[ibbc,1:je].=1.0;
  cbeybcr[ibbc,1:je].=0.0;

  dAux1 = (dx.*(collect(2:iebc).-0.5)).^(orderbc+1);
  dAux2 = (dx.*(collect(2:iebc).-1.5)).^(orderbc+1);
  sigmaPML = bcfactor.*(dAux1.-dAux2);
  aux1PML = exp.((-1).*sigmaPML.*(dt/epsz));
  aux2PML = ((1).-aux1PML)./(sigmaPML.*dx);
  caeybcr[2:iebc,1:je].= ones(1,je).*aux1PML;
  cbeybcr[2:iebc,1:je].= ones(1,je).*aux2PML;
  caeybcf[idy,idx].= ones(1,jebc).*aux1PML;
  cbeybcf[idy,idx].= ones(1,jebc).*aux2PML;
  caeybcb[idy,idx].= ones(1,jebc).*aux1PML;
  cbeybcb[idy,idx].= ones(1,jebc).*aux2PML;

  sigmaAux = bcfactor*(0.5*dx)^(orderbc+1);
  ca1=exp(-sigmaAux*dt/epsz);   cb1=(1-ca1)/(sigmaAux*dx);
  caey[1,1:je].=ca1;
  cbey[1,1:je].=cb1;
  caeybcf[iebc+ib,idx].=ca1;
  cbeybcf[iebc+ib,idx].=cb1;
  caeybcb[iebc+ib,idx].=ca1;
  cbeybcb[iebc+ib,idx].=cb1;

  # Campo magnético
    idx = 1:jebc;         idy = 1+iebc+ie:iebc+iebc+ie;
    dAux1 = (dx.*(collect(idx))).^(orderbc+1);
    dAux2 = (dx.*(collect(idx).-1)).^(orderbc+1);
    sigmaPML = bcfactor.*(muz/epsz).*(dAux1.-dAux2);
    aux1PML = exp.((-1).*sigmaPML.*(dt/muz));
    aux2PML = ((1).-aux1PML)./(sigmaPML.*dx);
    dahzxbcr[1:iebc,1:je].= ones(1,je).*aux1PML;
    dbhzxbcr[1:iebc,1:je].= ones(1,je).*aux2PML;
    dahzxbcf[idy,idx].= ones(1,jebc).*aux1PML;
    dbhzxbcf[idy,idx].= ones(1,jebc).*aux2PML;
    dahzxbcb[idy,idx].= ones(1,jebc).*aux1PML;
    dbhzxbcb[idy,idx].= ones(1,jebc).*aux2PML;

#

########################################################################
# Ecuaciones de propagación 
#######################################################################
println("Inicia Propagación")
@time begin

for nit=1:nmax

  # Campo eléctrico
  ex[:,2:je]=caex[:,2:je].*ex[:,2:je]+cbex[:,2:je].*(hz[:,2:je]-hz[:,1:je-1]);
  ey[2:ie,:]=caey[2:ie,:].*ey[2:ie,:]+cbey[2:ie,:].*(hz[1:ie-1,:]-hz[2:ie,:]);


  # Acutalización en la frontera Ex
  #     FRONT
  exbcf[:,2:jebc]=    caexbcf[:,2:jebc].*exbcf[:,2:jebc].-cbexbcf[:,2:jebc].*(hzxbcf[:,1:jebc-1]+hzybcf[:,1:jebc-1]-hzxbcf[:,2:jebc]-hzybcf[:,2:jebc]);
  ex[1:ie,1]=         caex[1:ie,1].*ex[1:ie,1].-cbex[1:ie,1].*(hzxbcf[ibbc:iebc+ie,jebc]+hzybcf[ibbc:iebc+ie,jebc]-hz[1:ie,1]);
  
  #     BACK
  exbcb[:,2:jebc-1]=  caexbcb[:,2:jebc-1].*exbcb[:,2:jebc-1].-cbexbcb[:,2:jebc-1].*(hzxbcb[:,1:jebc-2]+hzybcb[:,1:jebc-2]-hzxbcb[:,2:jebc-1]-hzybcb[:,2:jebc-1]);
  ex[1:ie,jb]=        caex[1:ie,jb].*ex[1:ie,jb].-cbex[1:ie,jb].*(hz[1:ie,jb-1]-hzxbcb[ibbc:iebc+ie,1]-hzybcb[ibbc:iebc+ie,1]);
  
  #     LEFT
  exbcl[:,2:je]=      caexbcl[:,2:je].*exbcl[:,2:je].-cbexbcl[:,2:je].*(hzxbcl[:,1:je-1]+hzybcl[:,1:je-1]-hzxbcl[:,2:je]-hzybcl[:,2:je]);
  exbcl[:,1]=         caexbcl[:,1].*exbcl[:,1].-cbexbcl[:,1].*(hzxbcf[1:iebc,jebc]+hzybcf[1:iebc,jebc]-hzxbcl[:,1]-hzybcl[:,1]);
  exbcl[:,jb]=        caexbcl[:,jb].*exbcl[:,jb].-cbexbcl[:,jb].*(hzxbcl[:,je]+hzybcl[:,je]-hzxbcb[1:iebc,1]-hzybcb[1:iebc,1]);
  
  #     RIGHT
  exbcr[:,2:je]=      caexbcr[:,2:je].*exbcr[:,2:je].-cbexbcr[:,2:je].*(hzxbcr[:,1:je-1]+hzybcr[:,1:je-1]-hzxbcr[:,2:je]-hzybcr[:,2:je]);
  exbcr[:,1]=         caexbcr[:,1].*exbcr[:,1].-cbexbcr[:,1].*(hzxbcf[1+iebc+ie:iefbc,jebc]+hzybcf[1+iebc+ie:iefbc,jebc]-hzxbcr[:,1]-hzybcr[:,1]);
  exbcr[:,jb]=        caexbcr[:,jb].*exbcr[:,jb].-cbexbcr[:,jb].*(hzxbcr[:,je]+hzybcr[:,je]-hzxbcb[1+iebc+ie:iefbc,1]-hzybcb[1+iebc+ie:iefbc,1]);
  

  # Acutalización en la frontera Ey
  #     FRONT
  eybcf[2:iefbc,:]=   caeybcf[2:iefbc,:].*eybcf[2:iefbc,:].-cbeybcf[2:iefbc,:].*(hzxbcf[2:iefbc,:]+hzybcf[2:iefbc,:]-hzxbcf[1:iefbc-1,:]-hzybcf[1:iefbc-1,:]);
  
  #     BACK
  eybcb[2:iefbc,:]=   caeybcb[2:iefbc,:].*eybcb[2:iefbc,:].-cbeybcb[2:iefbc,:].*(hzxbcb[2:iefbc,:]+hzybcb[2:iefbc,:]-hzxbcb[1:iefbc-1,:]-hzybcb[1:iefbc-1,:]);
  
  #     LEFT
  eybcl[2:iebc,:]=    caeybcl[2:iebc,:].*eybcl[2:iebc,:].-cbeybcl[2:iebc,:].*(hzxbcl[2:iebc,:]+hzybcl[2:iebc,:]-hzxbcl[1:iebc-1,:]-hzybcl[1:iebc-1,:]);
  ey[1,:]=            caey[1,:].*ey[1,:].-cbey[1,:].*(hz[1,:]-hzxbcl[iebc,:]-hzybcl[iebc,:]);
  
  #     RIGHT
  eybcr[2:iebc,:]=    caeybcr[2:iebc,:].*eybcr[2:iebc,:].-cbeybcr[2:iebc,:].*(hzxbcr[2:iebc,:]+hzybcr[2:iebc,:]-hzxbcr[1:iebc-1,:]-hzybcr[1:iebc-1,:]);
  ey[ib,:]=           caey[ib,:].*ey[ib,:].-cbey[ib,:].*(hzxbcr[1,:]+hzybcr[1,:]- hz[ie,:]);


  # Campo Magnético Hz
  hz[1:ie,1:je]=      dahz[1:ie,1:je].*hz[1:ie,1:je].+dbhz[1:ie,1:je].*(ex[1:ie,2:jb]-ex[1:ie,1:je]+ey[1:ie,1:je]-ey[2:ib,1:je]);
  #hz[is,js].=         hz[is,js].+sin.(omega*(nit.-delay)*dt).*exp.(-((nit.-delay).^2/tau^2));
  #hz[is,js].=         hz[is,js].+exp.(-((nit.-delay).^2/tau^2));
  hz[is,js].=         hz[is,js].+sin.(omega*(nit.-delay)*dt)

  # Actualización en la mfrontera Hzx
  #     FRONT
  hzxbcf[1:iefbc,:]=  dahzxbcf[1:iefbc,:].*hzxbcf[1:iefbc,:].-dbhzxbcf[1:iefbc,:].*(eybcf[2:ibfbc,:]-eybcf[1:iefbc,:]);

  #     BACK
  hzxbcb[1:iefbc,:]=  dahzxbcb[1:iefbc,:].*hzxbcb[1:iefbc,:].-dbhzxbcb[1:iefbc,:].*(eybcb[2:ibfbc,:]-eybcb[1:iefbc,:]);
 
  #     LEFT
  hzxbcl[1:iebc-1,:]= dahzxbcl[1:iebc-1,:].*hzxbcl[1:iebc-1,:].-dbhzxbcl[1:iebc-1,:].*(eybcl[2:iebc,:]-eybcl[1:iebc-1,:]);
  hzxbcl[iebc,:]=     dahzxbcl[iebc,:].*hzxbcl[iebc,:].-dbhzxbcl[iebc,:].*(ey[1,:]-eybcl[iebc,:]);
 
  #     RIGHT
  hzxbcr[2:iebc,:]=   dahzxbcr[2:iebc,:].*hzxbcr[2:iebc,:].-dbhzxbcr[2:iebc,:].*(eybcr[3:ibbc,:]-eybcr[2:iebc,:]);
  hzxbcr[1,:]=        dahzxbcr[1,:].*hzxbcr[1,:].-dbhzxbcr[1,:].*(eybcr[2,:]-ey[ib,:]);


  # Actualización en la frontera Hzy
  
  #     FRONT
  hzybcf[:,1:jebc-1]= dahzybcf[:,1:jebc-1].*hzybcf[:,1:jebc-1].-dbhzybcf[:,1:jebc-1].*(exbcf[:,1:jebc-1]-exbcf[:,2:jebc]);
  hzybcf[1:iebc,jebc]=dahzybcf[1:iebc,jebc].*hzybcf[1:iebc,jebc].-dbhzybcf[1:iebc,jebc].*(exbcf[1:iebc,jebc]-exbcl[1:iebc,1]);
  hzybcf[iebc+1:iebc+ie,jebc]=dahzybcf[iebc+1:iebc+ie,jebc].*hzybcf[iebc+1:iebc+ie,jebc].-dbhzybcf[iebc+1:iebc+ie,jebc].*(exbcf[iebc+1:iebc+ie,jebc]-ex[1:ie,1]);
  hzybcf[iebc+ie+1:iefbc,jebc]=dahzybcf[iebc+ie+1:iefbc,jebc].*hzybcf[iebc+ie+1:iefbc,jebc].-dbhzybcf[iebc+ie+1:iefbc,jebc].*(exbcf[iebc+ie+1:iefbc,jebc]-exbcr[1:iebc,1]);

  #     BACK
  hzybcb[1:iefbc,2:jebc]=dahzybcb[1:iefbc,2:jebc].*hzybcb[1:iefbc,2:jebc].-dbhzybcb[1:iefbc,2:jebc].*(exbcb[1:iefbc,2:jebc]-exbcb[1:iefbc,3:jbbc]);
  hzybcb[1:iebc,1]=dahzybcb[1:iebc,1].*hzybcb[1:iebc,1].-dbhzybcb[1:iebc,1].*(exbcl[1:iebc,jb]-exbcb[1:iebc,2]);
  hzybcb[iebc+1:iebc+ie,1]=dahzybcb[iebc+1:iebc+ie,1].*hzybcb[iebc+1:iebc+ie,1].-dbhzybcb[iebc+1:iebc+ie,1].*(ex[1:ie,jb]-exbcb[iebc+1:iebc+ie,2]);
  hzybcb[iebc+ie+1:iefbc,1]=dahzybcb[iebc+ie+1:iefbc,1].*hzybcb[iebc+ie+1:iefbc,1].-dbhzybcb[iebc+ie+1:iefbc,1].*(exbcr[1:iebc,jb]-exbcb[iebc+ie+1:iefbc,2]);
  
  #     LEFT
  hzybcl[:,1:je]=dahzybcl[:,1:je].*hzybcl[:,1:je]-dbhzybcl[:,1:je].*(exbcl[:,1:je]-exbcl[:,2:jb]);
  
  #     RIGHT
  hzybcr[:,1:je]=dahzybcr[2:ibbc,1:je].*hzybcr[:,1:je].-dbhzybcr[2:ibbc,1:je].*(exbcr[:,1:je]-exbcr[:,2:jb]);

  # Guardar matrices de infromación
  save(File(format"JLD",string(pwd(),"\\gifMagico\\matHz-",nit,".jld")),"hz",hz)

end

end
