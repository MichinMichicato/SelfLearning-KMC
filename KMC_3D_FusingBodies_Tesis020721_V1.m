%Algorithm Fortran 95 - Last Version to New Version Matlab
%Revisar Algoritmo de creacion de tamaño de lista y llenado de eventos
%cuando hay celulas acercandose a las fronteras.
%-/-/-/USES FusingSpheresAdjust, All3DEval, All2DEvalAdHoc /-/-/-/-/-/-/-/
%-/-/-/USES EvalMorphPoint, EvalUniform /-/-/-/-/-/-/-/

global E;

%Parametros Principales
ta0=double(1); %Ta0 as the inverse of the characteristic frequency, it can be interpreted to have a value of 1 in t0 units [its value does not affect the overall algorithm nor the relative probs but does affect single rates of occurrance - let as 1 so that tao of fusion process can be found through best fitted analytical equation // DO NOT CONFUSE WITH THE TIME UNITS OF THE SIMULATION     "t0"]
Ycm=double(0.5);
t0=7.40490629219698;%Equivalent value of a t0 unit in minutes [per transivity it'd be the value of ta0]
redhot=0;

%Declaracion de otras variables
time=double(0);
treal = double(0);%In MINUTES
Even=double.empty;%Even = zeros (3,3)%zeros(:,8));% - Lista de eventos
%Posi = int32(2,2) - Lista coordenadas de celulas en el sistema
%nevent = double(0);
%numcel= double(0);
tamanolista= int32(0);
ncm= double(0);
nsn=double(26);%numero de vecinos significativos
Rtot=double(0);
%Pacum=double(0);
E1=double(0);
cotainf = int32(1);%Cota inf para busqueda logaritmica en lista de eventos (Occurring event selection)
cotasup = int32(1);%Cota sup para busqueda logaritmica en lista de eventos (Occurring event selection)
eval= int32(1);%Evaluation Index over the list of events for the logarithmic search (Occurring event selection)
K_esimo= int32(1); %Result Index  for the logarithmic search over 'Even' list (Occurring event selection)
dif=0;%Parameter criterium for the convergence of logarithmic search over 'Even'
iter=0;%Counter for the MCS number of current iteration
numr=1;%Index for registration of Fusing_Spheres output parameters
count=0;%Auxiliar Parameter for registration of Fusing_Spheres function to decide when to register
register=0;%Auxliliar Parameter to determine when to register (expressed in minutes);
Parameters=double.empty;
Results=double.empty;
SystemRange=double.empty; %Min, Max, CM, Sig position for Validating System
SystemUniformity=double.empty; %Uniformity Indeces for Validating System
ValidatingPoint=zeros(5,3);%Matrix containing the x,y,z, of the 5 points along the symmetry axis that are used for mapping structure properties
SystemVar=double.empty;
timepicture=0;
PictureE=double.empty;
%-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

%%CREACION DEL SISTEMA (DEFINICION INICIAL DEL ESPACIO)

%system Style define what system is created, follow the list to pick one
%1: Fusing Spheres (Cell Aggregates) - Typically for calibration
%21: Torus arrangement of spherical aggregates with hand defined coordinates for the aggregates - For when
%validating a specific system
%22: Cyllinder arrangement of spherical aggregates with hand defined rows positions - For when validating a
%specific system
%31: Analitical defined cyllinder
%32: Analitical defined Torus
%33: Torus arrangement of spherical aggregates with parametrized (tunable)
%position for the aggregates

SystemStyle=33;

switch( SystemStyle)


    case 1
%% 1.- Fusing Spheres - Calibration // creando sistema (ubicación de células)

tic
%global E;

%N tamaño del espacio, tsimulacion tiempo max a simular
N=60; %N shall be at least equal to sum of radiuses +4 (due to frontier clearance condition)
tsimulacion=0.7;

%Creando Matriz del espacio
E = single(int32(zeros(N,N,N)));

%Cuadrado de los Radios de los N agregados
RadLeft=10;
RadRight=10;
Rad1=RadLeft^2;
Rad2=RadRight^2;

%Ubicación Centros de las N esferas (agregados)
Cx1 = int32(N/2)-(Rad1^0.5)-0.25;
Cy1 = int32(N/2) ;%-Rad1)
Cz1 = int32((N/2));
Cx2 = int32(N/2)+(Rad2^0.5)+0.25;
Cy2 = int32((N/2)) ;%+Rad2
Cz2 = int32((N/2));


for i=1:N
    for o = 1:N
        for z=1:N
            if ((i-Cy1)^2+(o-Cx1)^2+(z-Cz1)^2)<= Rad1 || ((i-Cy2)^2+(o-Cx2)^2+(z-Cz2)^2)<= Rad2
                E(i,o,z)=1;
            end
        end
    end
end
%}


%% 2.- %% Algorithm  for creating the system that will validate our model that is the ring of spherical aggregates and the cylinder of aggregates that shall fuse into a toroid and a pipe
    case 21
%      2.1 PART 1 - TOROIDAL RING

%
tic
%global E;

%N tamaño del espacio, tsimulacion tiempo max a simular
N=200; %N shall be at least equal to sum of radiuses +4 (due to frontier clearance condition)
tsimulacion=0.7;
NumSpheres=4;

%Coordinates over a plane for the CM of the spheres (Aggregates)
CentX=double.empty;
CentY=double.empty;



%Creando Matriz del espacio
E = single(int32(zeros(N,N,N)));

%Cuadrado del Radio de los NumSpheres agregados / Assumes one single raidus
%for all

Rad=15;
Rad1=Rad^2;


%Ubicación Centros de las N esferas (agregados)


CentX(1)= 40;
CentY(1)= 40;

CentX(2)= 70;
CentY(2)= 40;

CentX(3)= 40;
CentY(3)= 70;

CentX(4)= 70;
CentY(4)= 70;


%Defines coordinates of the plane on which all aggregates' centers lie

CentZ=int32(N/2);

%creando sistema (ubicación de células)

    %To better narrow the system and fasten the process we reduce the
    %limits into which the system creation loop evaluates the lattices ->
    %Taking the most "excentric" center in that coordinate and
    %adding/substracting the spheres' radius+1
    
    
    %Time went from 31 sec to 51 sec for 4 spheres of rad=10, but could
    %enhance time for larger systems / keep it and check later when we have
    %actual system coordinates
     %{
    minI = min(CentY)-(Rad+1);
    maxI = max(CentY)+(Rad+1);
    minO = min(CentX)-(Rad+1);
    maxO = max(CentX)+(Rad+1);
    minZ = CentZ-(Rad+1);
    maxZ = CentZ+(Rad+1);
    %}
    %Ok this should make sense at least
    minZ = CentZ-(Rad+1);
    maxZ = CentZ+(Rad+1);
    
for i=1:N
    for o = 1:N
        for z=minZ:maxZ
            for counterSpheres=1:NumSpheres
                if ((i-CentY(counterSpheres))^2+(o-CentX(counterSpheres))^2+(z-CentZ)^2)<= Rad1 
                    E(i,o,z)=1;
            end
        end
    end
end

end
%Generando Sistema

%Registro sistema inicial en una matriz dentro del programa (para acceso
%rapido)

%EINI=E;

%2.1
%{
%ValidatingPoint(1,3)=(N/2)-Rad+2;
%ValidatingPoint(2,3)=(N/2)-(Rad/2);
%ValidatingPoint(3,3)=(N/2);
%ValidatingPoint(4,3)=(N/2)+(Rad/2);
%ValidatingPoint(5,3)=(N/2)+(Rad/2)-2;

%}

    case 22
%      2.2 PART 2 - CYLLINDER

tic
%global E;

%Here we will assume the spheres are aligned in parallels around the axis
%of the cylinder (i.e. no X-degree twisted arrangement)

%Rationale, the cylinder as a tower thus it has levels and are all
%identical, one defines the coordinates only in x & y for the main columns

%N tamaño del espacio, tsimulacion tiempo max a simular
N=400; %N shall be at least equal to sum of radiuses +4 (due to frontier clearance condition)
tsimulacion=0.7;
NumSpheres=16; %Total in the system
NumAggreg=4;  %Per level

NumLevels=NumSpheres/NumAggreg; %Number of levels

%Coordinates over a plane for the CM of the spheres (Aggregates)
CentX=double.empty;
CentY=double.empty;



%Creando Matriz del espacio
E = single(int32(zeros(N,N,N)));

%Cuadrado del Radio de los NumSpheres agregados / Assumes one single raidus
%for all

Rad=15;
Rad1=Rad^2;


%Ubicación Centros de las N esferas (agregados)

    %Only define coordinates in x and y for rows conforming the cylinder,
    %the z coordinate is given as z0+(2*rad) for the next slide

CentX(1)= 40;
CentY(1)= 40;
CentX(2)= 70;
CentY(2)= 40;
CentX(3)= 40;
CentY(3)= 70;
CentX(4)= 70;
CentY(4)= 70;


%Defines coordinates of the  first plane of aggregates conforming the
%cyllinder

CentZini=(N/2)-((NumLevels/2)*(2*Rad))+Rad;
CentZ=double.empty;
CentZ(1)=CentZini;
for i=2:NumLevels
    CentZ(i)= CentZ(i-1)+(2*Rad);
end


%creando sistema (ubicación de células)

    %To better narrow the system and fasten the process we reduce the
    %limits into which the system creation loop evaluates the lattices ->
    %Taking the most "excentric" center in that coordinate and
    %adding/substracting the spheres' radius+1
    
    
    %Time went from 31 sec to 51 sec for 4 spheres of rad=10, but could
    %enhance time for larger systems / keep it and check later when we have
    %actual system coordinates
     %{
    minI = min(CentY)-(Rad+1);
    maxI = max(CentY)+(Rad+1);
    minO = min(CentX)-(Rad+1);
    maxO = max(CentX)+(Rad+1);
    minZ = CentZ-(Rad+1);
    maxZ = CentZ+(Rad+1);
    %}
    %Ok this should make sense at least
   % minZ = CentZ-(Rad+1);
   % maxZ = CentZ+(Rad+1);
    
for i=1:N
    for o = 1:N
        for z=1:N
            for counterSpheres=1:NumAggreg
                for level=1:NumLevels
                if ((i-CentY(counterSpheres))^2+(o-CentX(counterSpheres))^2+(z-CentZ(level))^2)<= Rad1 
                    E(i,o,z)=1;
                end
            end
        end
    end
end

end
%Generando Sistema

%Registro sistema inicial en una matriz dentro del programa (para acceso
%rapido)

%EINI=E;

%2.2
%{
%ValidatingPoint(1,3)=(N/2)-((NumLevels/2)*(2*Rad))+2;
%ValidatingPoint(2,3)=(N/2)-((NumLevels/2)*(Rad));
%ValidatingPoint(3,3)=(N/2);
%ValidatingPoint(4,3)=(N/2)+(NumLevels*Rad/2);
%ValidatingPoint(5,3)=(N/2)+((NumLevels/2)*(2*Rad))-2;

%}

%% 3.- %%Analytical (Customizable Torus w/ shperoids) System
%Creates a solid cylinder (Part 1) & Toroid (Part 2) from analytical equations

    %This system is to be used to test out sampling functions and be able
    %to interprete their results in the validation simulations

    case 31
%       3.1 PART 1   - CYLLINDER
    
    tic
%global E;

%Here we will assume the spheres are aligned in parallel column lines around the main axis
%of the cylinder (i.e. no X-degree twisted arrangement)
    %Rationale, the cylinder as a tower thus it has levels and are all
    %identical, one defines the coordinates only in x & y for the main columns

%N tamaño del espacio, tsimulacion tiempo max a simular
N=300; %N shall be at least equal to sum of radiuses +4 (due to frontier clearance condition)
NumLevels=90; %Number of levels (expressed in cell diameters)

%Defines the inner and the outter radius
RadOut=100; %Outter radius
RadIn=85; %Inner radius
Rad1=RadOut^2;
Rad2=RadIn^2;

%Coordinates on x-y for the axis of the cylinder
CentX=(N/2);
CentY=(N/2);
CentZini=(N/2)-(NumLevels/2);
CentZfinal=CentZini+NumLevels;

%Creando Matriz del espacio
E = single(int32(zeros(N,N,N)));

%creando sistema (ubicación de células)

    %To better narrow the system and fasten the process we reduce the
    %limits into which the system creation loop evaluates the lattices ->
    %Taking the most "excentric" center in that coordinate and
    %adding/substracting the spheres' radius+1
    
    
    %Time went from 31 sec to 51 sec for 4 spheres of rad=10, but could
    %enhance time for larger systems / keep it and check later when we have
    %actual system coordinates
     %{
    minI = min(CentY)-(Rad+1);
    maxI = max(CentY)+(Rad+1);
    minO = min(CentX)-(Rad+1);
    maxO = max(CentX)+(Rad+1);
    minZ = CentZ-(Rad+1);
    maxZ = CentZ+(Rad+1);
    %}
    %Ok this should make sense at least
   % minZ = CentZ-(Rad+1);
   % maxZ = CentZ+(Rad+1);
    
for i=1:N
    for o = 1:N
        for z=CentZini:CentZfinal
                if ((i-CentY)^2+(o-CentX)^2)<= Rad1 && ((i-CentY)^2+(o-CentX)^2)>= Rad2 
                    E(i,o,z)=1;
        end
    end
end

end
%Generando Sistema

%Registro sistema inicial en una matriz dentro del programa (para acceso
%rapido)

%EINI=E;

%3.1
%{
%ValidatingPoint(1,3)=(N/2)-(NumLevels/2)+2;
%ValidatingPoint(2,3)=(N/2)-(NumLevels/4);
%ValidatingPoint(3,3)=(N/2);
%ValidatingPoint(4,3)=(N/2)+(NumLevels/4);
%ValidatingPoint(5,3)=(N/2)+((NumLevels/2)-2;

    %}
    

    case 32
%       3.2 PART 2 - TORUS


tic
%global E;
%N tamaño del espacio, tsimulacion tiempo max a simular
N=50; %N shall be at least equal to sum of radiuses +4 (due to frontier clearance condition)

%Defines the inner and the outter radius
Rad=20; %Radius of Torus (from pipe center to torus center)
RadIn=9; %Inner radius of pipe (radius of circle)
Rad1=Rad^2;
Rad2=RadIn^2;
AspectRatio=Rad/RadIn;

%Coordinates on x-y-z for the center of the Torus
CentX=N/2;
CentY=N/2;
CentZ=N/2;

%Creando Matriz del espacio
E = single(int32(zeros(N,N,N)));

%creando sistema (ubicación de células)

    %To better narrow the system and fasten the process we reduce the
    %limits into which the system creation loop evaluates the lattices ->
    %Taking the most "excentric" center in that coordinate and
    %adding/substracting the spheres' radius+1
    
    
    %Time went from 31 sec to 51 sec for 4 spheres of rad=10, but could
    %enhance time for larger systems / keep it and check later when we have
    %actual system coordinates
     %{
    minI = min(CentY)-(Rad+1);
    maxI = max(CentY)+(Rad+1);
    minO = min(CentX)-(Rad+1);
    maxO = max(CentX)+(Rad+1);
    minZ = CentZ-(Rad+1);
    maxZ = CentZ+(Rad+1);
    %}
    %Ok this should make sense at least
   % minZ = CentZ-(Rad+1);
   % maxZ = CentZ+(Rad+1);

   delta=RadIn;
   
   
for i=1:N
    for o = 1:N
        for z=1:N
                if ((((i-CentY)^2+(o-CentX)^2)^0.5)-Rad)^2+((z-CentZ)^2)-Rad2<=0
                    E(i,o,z)=1;
        end
    end
end
end
 
%3.2
%{
%ValidatingPoint(1,3)=(N/2)-Rad+2;
%ValidatingPoint(2,3)=(N/2)-(Rad/2);
%ValidatingPoint(3,3)=(N/2);
%ValidatingPoint(4,3)=(N/2)+(Rad/2);
%ValidatingPoint(5,3)=(N/2)+(Rad/2)-2;

   

%}


    case 33
%   3.3    PART 2 - TORUS Customizable Spheroids


tic

%for rep = 1:1

%Algorithm Fortran 95 - Last Version to New Version Matlab
%Revisar Algoritmo de creacion de tamaño de lista y llenado de eventos
%cuando hay celulas acercandose a las fronteras.


%N tamaño del espacio, tsimulacion tiempo max a simular
N=200;

%Creando Matriz del espacio
%global E;
E = int32(zeros(N,N,N));
NAg=10; %Numero de agregados en el anillo
RR=55.39;%Radio Anillo sobre el que estan centros de los agregados (Radio mayor Toroide)


%Cuadrado de los Radios de los N agregados
Rad=16.67;%Radio tipo 1 de Agregados
%Rad2=12;%Radio tipo 2 de Agregados
CAg=(0:0:0); %Lista de coordenadas dx & dy de centros de agregados esfericos respecto a centro del macro-anillo


%Ubicación Centros de Agregados

%Centro macro-anillo 
Cx=(N/2);
Cy=(N/2);
Cz=(N/2);
%{
i=int32(Cx);
j=int32(Cy);
o=int32(Cz);
Cx=double(i);
Cy=double(j);
Cz=double(o);
%}
    
%Calculo coordenadas centro de agregados

dAngle=2*pi()/NAg;
i=1;
Angle=(pi()/2);
if Rad > dAngle*RR; %Ithink it would be the other way aroun tho' if rad > arch lenght then there is such risk
disp ('Possible Overlapping Aggregates!!!');
end

while Angle < 5*(pi()/2)

    CAg(i,1)= int32(RR*cos(Angle));
    CAg(i,2)=int32(RR*sin(Angle));
    CAg(i,3)=0;
    Angle=Angle+dAngle;
        i=i+1;
end

%Correccion-Traslacion de las coordenadas respecto al centro del espacio
%computacional (pues el ciclo while anterior define las coordenadas
%respecto al centro que se asume en el origen, pero el origen lo queremos
%en el centro del espacio reticular (x,y,z)=N/2).
CAg(:,1)=CAg(:,1)+Cx;
CAg(:,2)=CAg(:,2)+Cy;
CAg(:,3)=CAg(:,3)+Cz;


NAgEf=NAg; %Puede cambiarse por la longitud (num filas de CAg) para definir un cilindro.


%creando sistema (ubicación de células para cada agregado celular)
for cont=1:NAgEf
for i=1:N%CAg(cont,1)-2*Rad:CAg(cont,1)+2*Rad %This +- 2*Rad is made to "fasten the process yet avoid missing any cell due to too narrow space definition (in theory could be just 1 rad)
    for o = 1:N%CAg(cont,2)-2*Rad:CAg(cont,2)+2*Rad
        for z=1:N
            if ((i-CAg(cont,2))^2+(o-CAg(cont,1))^2+(z-CAg(cont,3))^2)<= (Rad)^2 %|| ((i-Cy2)^2+(o-Cx2)^2+(z-Cz2)^2)<= Rad2
                E(i,o,z)=1;
            end
        end
    end
end
end

%3.3
%{
%ValidatingPoint(1,3)=(N/2)-Rad+2;
%ValidatingPoint(2,3)=(N/2)-(Rad/2);
%ValidatingPoint(3,3)=(N/2);
%ValidatingPoint(4,3)=(N/2)+(Rad/2);
%ValidatingPoint(5,3)=(N/2)+(Rad/2)-2;

%}

end


ValidatingPoint(:,1:2)=int32(N/2);
switch (SystemStyle)
    case 1
        ValidatingPoint(:,3)=int32(N/2);
    case 21
        ValidatingPoint(1,3)=int32((N/2)-Rad+2);
        ValidatingPoint(2,3)=int32((N/2)-(Rad/2));
        ValidatingPoint(3,3)=int32((N/2));
        ValidatingPoint(4,3)=int32((N/2)+(Rad/2));
        ValidatingPoint(5,3)=int32((N/2)+(Rad/2)-2);
    case 22
        ValidatingPoint(1,3)= int32((N/2)-((NumLevels/2)*(2*Rad))+2);
        ValidatingPoint(2,3)=int32((N/2)-((NumLevels/2)*(Rad)));
        ValidatingPoint(3,3)=int32((N/2));
        ValidatingPoint(4,3)=int32((N/2)+(NumLevels*Rad/2));
        ValidatingPoint(5,3)=int32((N/2)+((NumLevels/2)*(2*Rad))-2);
    case 31
        ValidatingPoint(1,3)=int32((N/2)-(NumLevels/2)+2);
        ValidatingPoint(2,3)=int32((N/2)-(NumLevels/4));
        ValidatingPoint(3,3)=int32((N/2));
        ValidatingPoint(4,3)=int32((N/2)+(NumLevels/4));
        ValidatingPoint(5,3)=int32((N/2)+(NumLevels/2)-2);
    case 32
        ValidatingPoint(1,3)=int32((N/2)-Rad+2);
        ValidatingPoint(2,3)=int32((N/2)-(Rad/2));
        ValidatingPoint(3,3)=int32((N/2));
        ValidatingPoint(4,3)=int32((N/2)+(Rad/2));
        ValidatingPoint(5,3)=int32((N/2)+(Rad/2)-2);
    case 33
        ValidatingPoint(1,3)=int32((N/2)-Rad+2);
        ValidatingPoint(2,3)=int32((N/2)-(Rad/2));
        ValidatingPoint(3,3)=int32((N/2));
        ValidatingPoint(4,3)= int32((N/2)+(Rad/2));
        ValidatingPoint(5,3)=int32((N/2)+(Rad/2)-2);
end

%%
%Registro sistema inicial en una matriz dentro del programa (para acceso
%rapido)

EINI=E;
PictureE(1,:,:,:)=E;
%Registrando el espacio (sistema) en un archivo txt (long ago was
%dissactivated)
%{
    name = sprintf('Espa_inicial%d.txt',0);
    fileID=fopen(name,'w');
for p = 1:N
for k = 1:N
for h = 1:N

    if k ~=N
        if h~= N
            formatSpec = '%d ';
            fprintf(fileID,formatSpec,E(k,h,p));
        else
            formatSpec = '%d \r\n';
            fprintf(fileID,formatSpec,E(k,h,p));
        end
        
    else
        
        if h~= N
            formatSpec = '%d ';
            fprintf(fileID,formatSpec,E(k,h,p));
        else
            formatSpec = '%d \r\n\v';
            fprintf(fileID,formatSpec,E(k,h,p));
        end
    end
    
end
end
end
    fclose(fileID);
%}

%-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

% Blanqueo de fronteras pre-evaluacion
E(1,:,:)=0;
E(2,:,:)=0;
E(:,1,:)=0;
E(:,2,:)=0;
E(:,:,1)=0;
E(:,:,2)=0;
E(N,:,:)=0;
E(N-1,:,:)=0;
E(:,N,:)=0;
E(:,N-1,:)=0;
E(:,:,N)=0;
E(:,:,N-1)=0;

%-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

    numcel=0;
    nevent=0;
    Rtot=0;
    
   %Conteo numero total de celulas y determinacion de num total de eventos
   %totales posibles (can be inefficient algorithm) para tamño de listas
   
        %Probably it would be more efficient to turn this part down, allow
        %Eve=double.void;(same for Posi), and instead create a loop of
        %feasibility (once an event is identified we can proceed, no need
        %to further evaluate the whole space)
        
%{   
   for i = 2:N-1
       for o = 2:N-1
           for z=2:N-1
               
               if E(i,o,z) == 1
                   numcel = numcel+1;
               end
               
               for h = -1:1
                   for k = -1:1
                       for p=-1:1
                           if E(i,o,z)~=E(i+h,o+k,z+p)
                               nevent=nevent+1;
                           end
                       end
                   end
               end
           end
       end
   end
tamanolista=int32(nevent/2);


%Si lista vacia salir del ciclo
if tamanolista==0
    disp ('All cells have dissappeared at');
    treal;
	iter=0;
end


posi = int32 (zeros(numcel,3));
Even = single (zeros(tamanolista,10));

%}
 %Ojo tenemos un barrido que está ignorando los bordes (primera y segunda
 %linea de borde).
 
 %-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
 
Even(1,10)=0;
 %Identificacion eventos
     %gposi=1;
    %Creation of List of Events "Eve" where position of involved nodes, ncm, rates (ri), Rtot,
    %etc... It also fills Posi (position of cells) [Disabled for versions 2020]
  for i=3:(N-2)
     for o=3:(N-2)
         for z=3:(N-2)
         
             
        %{
        Register coordinates of cell in the list (vector) posi, index used is gposi     
        if E(i,o,z)==1
            posi(gposi,1)=i;
            posi(gposi,2)=o;
            posi(gposi,3)=z;
            gposi=gposi+1;
         end
         %}
             
        for h=-1:1
        for k=-1:1
        for p=-1:1
            %Evento Identificado
            if E(i,o,z) ~= E(i+h,o+k,z+p)
               rep=0; %Initializes counter for loop verifying whether the event has already been registered
               g=1; %Initializes index of the vector where events are registered
               
               %Busqueda en lista y escrita segun caso
               while rep == 0
                 if g==size(Even,1)%(Even(g,1)==0 && Even(g,2)==0 && Even(g,3)==0 && Even(g,4)==0&& Even(g,5)==0&& Even(g,6)==0) %Busqueda lleva a espacio vacio, se resgitra evento en lista
                   Even(g+1,1)=i;
                   Even(g+1,2)=o;
                   Even(g+1,3)=z;
                   Even(g+1,4)=i+h;
                   Even(g+1,5)=o+k;
                   Even(g+1,6)=z+p;
                           %disp ('Evento encontrado')
                           %eventonum=eventonum
                           %eventonum=eventonum+1;
                   %Registro de Numero de vecinos no homologos (homo-clase)
                   ncm=0;
                   for t=-1:1
                       for e=-1:1
                        for d=-1:1   %Ciclo para conteo de vecindario de particulas no homo-clase

                            if E(i,o,z)~=E(i+t,o+e,z+d)
                                ncm=ncm+1;
                            end
                           
                            if E(i+h,o+k,z+p)~= E(i+h+t,o+e+k,z+p+d)
                                ncm=ncm+1;
                            end
                            
                        end   
                       end
                   end
                   Even(g+1,7)=ncm/2; %Registro prom vecindario no homo-clase
                   Even(g+1,8)=(exp(-(Ycm*(nsn-(2.0*Even(g+1,7))+1.0))))/ta0; %Velocidad de evento (rate of occurrance)
                   Rtot=Rtot+Even(g+1,8); %Conteo para (sum over each r(i))Rtotal.
                   rep=1; %Exits the loop of verification
                
                 else
                   
                        if (((i == Even(g,1) && o== Even(g,2) && z== Even(g,3)) && (i+h == Even(g,4) && o+k==Even(g,5) && z+p==Even(g,6))) || ((i+h==Even(g,1) && o+k==Even(g,2)&& z+p==Even(g,3)) &&(i==Even(g,4) && o==Even(g,5) && z==Even(g,6))))
                             rep=1; %Exits the loop of verification when the two involved nodes are already present in an event (in any order)
                        end
                             g=g+1; %Advances the index over the list of events 'Even', to check next registered event row
                  
                 end
               end 
           end
        end
        end
        end 
        end
     end
 end
  %disp (Even)
 if (Even(1,:)==0)
   Even(1,:)=[]; 
 end
  
 
 %Filling Probability Columns in 'Eve' for the events
 Pacum=0;%Cummulative Prob
 for g=1:size(Even)
     Even(g,9)=Even(g,8)/Rtot; %Prob. for each event
     Even(g,10)=Pacum+Even(g,9); %Cummulative prob in the list until that event
     Pacum=Pacum+Even(g,9);
         %disp('Llenando lista')
         %disp(g)
         %disp(Rtot)
         %disp ( Even(g,7))
         %disp(Even(g,8))
 end
 

 %Occurring Event Selection
  E1=rand; %Generación de número aleatorio
 cotainf=1;
 cotasup=double(size(Even,1));
 K_esimo=cotasup;
 eval=int32(cotainf+((cotasup-cotainf)/2));
 dif=0;
  while dif == 0
    if (cotasup-cotainf) >1
        if E1== Even((eval),10)
                K_esimo=eval;
                dif=1;
               % disp('igualito')
               % eval=int32((cotasup-cotainf)/2)+cotainf           
        else
               if E1 > Even(eval,10)
                    cotainf = eval;
                    eval=int32(cotainf+((cotasup-cotainf)/2));
                    dif=0;
                 %   disp('mayor')
                 %   eval=cotainf + int32((cotasup-cotainf)/2)
               else
                    if E1 < Even(eval,10)
                         cotasup = eval;
                         eval=int32(cotainf+((cotasup-cotainf)/2));
                         dif=0;
                  %       disp('menor')
                  %       eval=cotainf+ int32((cotasup-cotainf)/2)
                    end
               end
        end
    else
      K_esimo=cotasup;
      %disp('igualito bueno no, pero casi')
      %K_esimo
      dif=1;
    end
  end
 
 
 coordinatesofevent=Even(K_esimo,1:6);%saves the coordinates of the involved particles in the vector 'coordinates ofevent'
 %Ejecución del evento
 spin=E(Even(K_esimo,1),Even(K_esimo,2),Even(K_esimo,3));
 E(Even(K_esimo,1),Even(K_esimo,2),Even(K_esimo,3))=E(Even(K_esimo,4),Even(K_esimo,5),Even(K_esimo,6));
 E(Even(K_esimo,4),Even(K_esimo,5),Even(K_esimo,6))= spin;

 
 %Avance Temporal
 %Revisar UNIDADES en este caso ... (para el artículo)
 E2=rand;
 dT=-(log(E2)/Rtot); %Time leap for the event (in ta0 units)
 time=time+dT; %Current time expressed in units of ta0 (charact time)
 treal=(t0*time); %Conversion to real time equivalent (minutes) of current time
 iter=iter+1;
 register=register+treal;
 timepicture=timepicture+(dT*t0);

%-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


%Algoritmo Principal
while treal<=10080 %treal in minutes %10080 min = 168 h

    
 
 % Blanqueo de fronteras pre-evaluacion
E(1,:,:)=0;
E(2,:,:)=0;
E(:,1,:)=0;
E(:,2,:)=0;
E(:,:,1)=0;
E(:,:,2)=0;
E(N,:,:)=0;
E(N-1,:,:)=0;
E(:,N,:)=0;
E(:,N-1,:)=0;
E(:,:,N)=0;
E(:,:,N-1)=0;

%-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
            
%Make Sure there are events to happen on the list (i.e. the system has not
%desintegrated yet)

if size(Even,1)<1
   break
end
             

   %Evaluar Espacio modificado & Actualizar Lista 'Even'
  
   %Deletes the events from the list that involved any of the particles
   %that took part in the event executed
   g=1;
    endofroad=0;
    ficlos=0;
   while endofroad ~= 1
       
       if g > size(Even,1)
           endofroad=1;
       else
                 ficlos=ficlos+1;     
          if (((Even(g,1)==coordinatesofevent(1) && Even(g,2)==coordinatesofevent(2) && Even(g,3)==coordinatesofevent(3)))||(Even(g,1)==coordinatesofevent(4) && Even(g,2)==coordinatesofevent(5) && Even(g,3)==coordinatesofevent(6)))||((Even(g,4)==coordinatesofevent(1) && Even(g,5)==coordinatesofevent(2) && Even(g,6)==coordinatesofevent(3))||(Even(g,4)==coordinatesofevent(4) && Even(g,5)==coordinatesofevent(5) && Even(g,6)==coordinatesofevent(6))) %If any (other) event that involves any of the two selected (inolved) particles in the picked and executed event
                Even(g,:)=[];
            else
                g=g+1;
           end
               
            
      end
       
   end 
    
 
%Evaluate the local neighborhood in E to write the new events after occurrance at the end of the
%list 'Even'

                        %for i=min(coordintaesofevent(1),coordintaesofevent(4))-1:max(coordintaesofevent(1),coordintaesofevent(4))+1 %Evaluate in the cube formed by most extreme points resulting from evaluating y coordinates for both particles -1 and most extreme +1 (same for subsequent dimensions x,z)
                         %    for o=min(coordintaesofevent(2),coordintaesofevent(5))-1:max(coordintaesofevent(2),coordintaesofevent(5))+1
                          %       for z=min(coordintaesofevent(3),coordintaesofevent(6))-1:max(coordintaesofevent(3),coordintaesofevent(6))+1
deadE=(size(Even,1));
deadend=int64(deadE);
  for particle =0:1
    
    if particle ==0
        i=coordinatesofevent(1);
        o=coordinatesofevent(2);
        z=coordinatesofevent(3);
    else
        i=coordinatesofevent(4);
        o=coordinatesofevent(5);
        z=coordinatesofevent(6);
    end
    
        for h=-1:1
        for k=-1:1
        for p=-1:1
            %Evento Identificado
            if E(i,o,z) ~= E(i+h,o+k,z+p)
               rep=0; %Initializes counter for loop verifying whether the event has already been registered
               g=deadend; %Initializes index of the vector where events are registered
               
               %Busqueda en lista y escrita segun caso
               while rep == 0
                   
                 if  g==size(Even,1) %Busqueda lleva a espacio vacio, se resgitra evento en lista
                   Even(g+1,1)=i;
                   Even(g+1,2)=o;
                   Even(g+1,3)=z;
                   Even(g+1,4)=i+h;
                   Even(g+1,5)=o+k;
                   Even(g+1,6)=z+p;                           
                   %Registro de Numero de vecinos no homologos (homo-clase)
                   ncm=0;
                   for t=-1:1
                       for e=-1:1
                        for d=-1:1   %Ciclo para conteo de vecindario de particulas no homo-clase

                            if E(i,o,z)~=E(i+t,o+e,z+d)
                                ncm=ncm+1;
                            end
                           
                            if E(i+h,o+k,z+p)~= E(i+h+t,o+e+k,z+p+d)
                                ncm=ncm+1;
                            end
                            
                        end   
                       end
                   end
                   Even(g+1,7)=ncm/2; %Registro prom vecindario no homo-clase
                   Even(g+1,8)=(exp(-(Ycm*(nsn-(2.0*Even(g+1,7))+1.0))))/ta0; %Velocidad de evento (rate of occurrance)
                   %Rtot=Rtot+Even(g,8); %Conteo para (sum over each r(i))Rtotal.
                   rep=1; %Exits the loop of verification
                
                 else
                                        
                    if (((i == Even(g,1) && o== Even(g,2) && z== Even(g,3)) && (i+h == Even(g,4) && o+k==Even(g,5) && z+p==Even(g,6))) || ((i+h==Even(g,1) && o+k==Even(g,2)&& z+p==Even(g,3)) &&(i==Even(g,4) && o==Even(g,5) && z==Even(g,6))))
                         rep=1; %Exits the loop of verification when the two involved nodes are already present in an event (in any order)
                    end
                         g=g+1; %Advances the index over the list of events 'Even', to check next registered event row

                 end 
               end
            end
        end
        end 
        end
  
  end

 
%Filling Probability Columns in 'Eve' for the events
Rtot=0; 
for g=1:size(Even)
Rtot=Rtot+Even(g,8);
 end

 Pacum=0;%Cummulative Prob
 for g=1:size(Even)
     Even(g,9)=Even(g,8)/Rtot; %Prob. for each event
     Even(g,10)=Pacum+Even(g,9); %Cummulative prob in the list until that event
     Pacum=Pacum+Even(g,9);
 end
 

 
 
 %Occurring Event Selection
 E1=rand; %Generación de número aleatorio
 cotainf=1;
 cotasup=double(size(Even,1));
 K_esimo=cotasup;
 eval=int32(cotainf+((cotasup-cotainf)/2));
 dif=0;
 while dif == 0
    if (cotasup-cotainf) >1
        if E1== Even((eval),10)
                K_esimo=eval;
                dif=1;
               % disp('igualito')
               % eval=int32((cotasup-cotainf)/2)+cotainf           
        else
               if E1 > Even(eval,10)
                    cotainf = eval;
                    eval=int32(cotainf+((cotasup-cotainf)/2));
                    dif=0;
                 %   disp('mayor')
                 %   eval=cotainf + int32((cotasup-cotainf)/2)
               else
                    if E1 < Even(eval,10)
                         cotasup = eval;
                         eval=int32(cotainf+((cotasup-cotainf)/2));
                         dif=0;
                  %       disp('menor')
                  %       eval=cotainf+ int32((cotasup-cotainf)/2)
                    end
               end
        end
    else
      K_esimo=cotasup;
      %disp('igualito bueno no, pero casi')
      %K_esimo
      dif=1;
    end
  end
 
 coordinatesofevent=Even(K_esimo,1:6);%saves the coordinates of the involved particles in the vector 'coordinates ofevent'
 
 
 
 
 %Ejecución del evento
 spin=E(Even(K_esimo,1),Even(K_esimo,2),Even(K_esimo,3));
 E(Even(K_esimo,1),Even(K_esimo,2),Even(K_esimo,3))=E(Even(K_esimo,4),Even(K_esimo,5),Even(K_esimo,6));
 E(Even(K_esimo,4),Even(K_esimo,5),Even(K_esimo,6))= spin;

 
 %Avance Temporal
 %Revisar UNIDADES en este caso ... (para el artículo)
 E2=rand;
 dT=-(log(E2)/Rtot); %Time leap for the event (in t0 units)
 time=time+dT; %Current time expressed in units of ta0 (charact time)
 treal=(t0*time); %Conversion to real time equivalent (minutes) of current time

 
 %Blanqueo de fronteras post suceso
%{
E(1,:,:)=0;
E(2,:,:)=0;
E(:,1,:)=0;
E(:,2,:)=0;
E(:,:,1)=0;
E(:,:,2)=0;
E(N,:,:)=0;
E(N-1,:,:)=0;
E(:,N,:)=0;
E(:,N-1,:)=0;
E(:,:,N)=0;
E(:,:,N-1)=0; 
%}

 
 %{
%Registrar CM (Hay que definir cada cuantos pasos se va a registrar)
%Podria omitir el conteo de M de nuevo y basarme solo en numcel que fueron
%contadas anteriormente.

%Grabar espacio cada intervalo fijo
%{
%lap=lap+1;
%if (lap==333|lap==30|lap==100|lap==300|lap==1000|lap==3000|lap==10000|lap==30000|lap==33333)
       
   % name= sprintf('Centros De Masa_%d',rep);
    %fileID=fopen(name,'w')
%CMy=0;
%CMx=0;
%M=0;
%for i=1:N
%   for o=1:N
  %     if E(i,o)==1
   %         M=M+1;
   %         CMy=CMy+i;
 %           CMx=CMx+o;
  %     end
  % end
%end
%CMy=CMy/M;
%CMx=CMx/M;
                                %formatSpec='%d';
                                %dlmwrite(fileID,CMy,CMx,'-apend')
                                %fclose(fileID);
%         name = sprintf('Espa_%d_%d.txt',rep,lap);
%         fileID=fopen(name,'w');
%   for k = 1:N
%     for h = 1:N
%         if h~= N
%             formatSpec = '%d ';
%             fprintf(fileID,formatSpec,E(k,h));
%         else
%             formatSpec = '%d \r\n';
%             fprintf(fileID,formatSpec,E(k,h));
%         end
%     end
%   end
% fclose(fileID);

%Function to record on a logarithmic basis the events (having about the same amount of points per magnitude order (i.e per log10 unit)(

%fileID=fopen('Parameters','w');


%f=floor((log10(iter)+1)^(log10(iter)+1));
%}
%}

 
 
%RE- ACTIVAR PARA REGISTRAR PARAMETROS DE SIMULACION




switch (SystemStyle)
%Parameters Registration for Fusing Spheres - Produce two vectors that lead
%to XLS files (Parameters & Results).

    case 1

f=floor(log10(iter));
if f<6
if count<1
	count=count+(1/(10^f));
	else
    [OptRadR,OptRadL,Rneck,RneckB,RneckC,RneckD,TetaA,TetaB,TetaC,Xinter,DistCM,VoidR,VoidL] = FusingSpheresAdjust(N);
	%treal,OptRadR,Rneck,TetaB
    count=0;
    Parameters(numr,:)= [iter,dT,time,treal,OptRadR,OptRadL,Rneck,RneckB,RneckC,RneckD,TetaA,TetaB,TetaC,DistCM, VoidR, VoidL];
    Results(numr,:)=[iter,dT,time,treal,OptRadR,OptRadL,Rneck,RneckB,RneckC,DistCM];
    %formatSpec='%d,';
    %dlmwrite(fileID,iter,treal,OptRadL,Rneck,TetaB,'-apend');
    numr=numr+1;
end

else
    if count<1
	count=count+(1/(10^(f-1)));
	else
    [OptRadR,OptRadL,Rneck,RneckB,RneckC,RneckD,TetaA,TetaB,TetaC,Xinter,DistCM,VoidR,VoidL] = FusingSpheresAdjust(N);
	%treal,OptRadR,Rneck,TetaB
    count=0;
    Parameters(numr,:)= [iter,dT,time,treal,OptRadR,OptRadL,Rneck,RneckB,RneckC,RneckD,TetaA,TetaB,TetaC,DistCM, VoidR, VoidL];
    Results(numr,:)=[iter,dT,time,treal,OptRadR,OptRadL,Rneck,RneckB,RneckC,DistCM];
    %formatSpec='%d,';
    %dlmwrite(fileID,iter,treal,OptRadL,Rneck,TetaB,'-apend');
    numr=numr+1;
end
    
end
%{
%fclose(fileID);    

    [OptRadR,OptRadL,Rneck,RneckB,RneckC,RneckD,TetaA,TetaB,TetaC,Xinter,DistCM,VoidR,VoidL] = FusingSpheresAdjust(N);
	Parameters(numr,:)= [iter,dT,time,treal,OptRadR,OptRadL,Rneck,RneckB,RneckC,RneckD,TetaA,TetaB,TetaC,DistCM];
    Results(numr,:)=[iter,dT,time,treal,OptRadR,OptRadL,Rneck,RneckB,RneckC,DistCM];
    numr=numr+1;

%}



        
    otherwise
%Parameters Registration for Validating System
%It is necessary to select the coordinates for evaluation and keep them
%fixed -> Better if defined after system creation and then referenced here
%(5 points are suggested)

if register >120 || iter==1 %Time in Minutes 
    %Determine the time interval in minutes for the measurement/scanning (determine how many minutes between measurements) 
             
    [Xmin_1,Ymin_1,Xmax_1,Ymax_1,XCM_1,YCM_1,Sig_1]= EvalMorphPoint (ValidatingPoint(1,1),ValidatingPoint(1,2),ValidatingPoint(1,3),N);
    [Xmin_2,Ymin_2,Xmax_2,Ymax_2,XCM_2,YCM_2,Sig_2]= EvalMorphPoint (ValidatingPoint(2,1),ValidatingPoint(2,2),ValidatingPoint(2,3),N);
    [Xmin_3,Ymin_3,Xmax_3,Ymax_3,XCM_3,YCM_3,Sig_3]= EvalMorphPoint (ValidatingPoint(3,1),ValidatingPoint(3,2),ValidatingPoint(3,3),N);
    [Xmin_4,Ymin_4,Xmax_4,Ymax_4,XCM_4,YCM_4,Sig_4]= EvalMorphPoint (ValidatingPoint(4,1),ValidatingPoint(4,2),ValidatingPoint(4,3),N); 
    [Xmin_5,Ymin_5,Xmax_5,Ymax_5,XCM_5,YCM_5,Sig_5]= EvalMorphPoint (ValidatingPoint(5,1),ValidatingPoint(5,2),ValidatingPoint(5,3),N);     
   % [Xmin_6,Ymin_6,Xmax_6,Ymax_6,XCM_6,YCM_6,Sig_6]= EvalMorphPoint (ValidatingPoint(6,1),ValidatingPoint(6,2),ValidatingPoint(6,3),N);
   % [Xmin_7,Ymin_7,Xmax_7,Ymax_7,XCM_7,YCM_7,Sig_7]= EvalMorphPoint (ValidatingPoint(7,1),ValidatingPoint(7,2),ValidatingPoint(7,3),N);
   % [Xmin_8,Ymin_8,Xmax_8,Ymax_8,XCM_8,YCM_8,Sig_8]= EvalMorphPoint (ValidatingPoint(8,1),ValidatingPoint(8,2),ValidatingPoint(8,3),N);
    SystemRange(numr,:)= [Xmin_1,Ymin_1,Xmax_1,Ymax_1,XCM_1,YCM_1,Xmin_2,Ymin_2,Xmax_2,Ymax_2,XCM_2,YCM_2,Xmin_3,Ymin_3,Xmax_3,Ymax_3,XCM_3,YCM_3,Xmin_4,Ymin_4,Xmax_4,Ymax_4,XCM_4,YCM_4,Xmin_5,Ymin_5,Xmax_5,Ymax_5,XCM_5,YCM_5];%,Xmin_6,Ymin_6,Xmax_6,Ymax_6,XCM_6,YCM_6,Sig_6,Xmin_7,Ymin_7,Xmax_7,Ymax_7,XCM_7,YCM_7,Sig_7,Xmin_8,Ymin_8,Xmax_8,Ymax_8,XCM_8,YCM_8,Sig_8];    
    SystemVar(numr,:,:)=[Sig_1,Sig_2,Sig_3,Sig_4,Sig_5];
    [UG_1, UCM_1] = EvalUniform (Xmin_1,Ymin_1,Xmax_1,Ymax_1,XCM_1,YCM_1,Sig_1,ValidatingPoint(1,1),ValidatingPoint(1,2),ValidatingPoint(1,3));
    [UG_2, UCM_2] = EvalUniform (Xmin_2,Ymin_2,Xmax_2,Ymax_2,XCM_2,YCM_2,Sig_2,ValidatingPoint(2,1),ValidatingPoint(2,2),ValidatingPoint(2,3));
    [UG_3, UCM_3] = EvalUniform (Xmin_3,Ymin_3,Xmax_3,Ymax_3,XCM_3,YCM_3,Sig_3,ValidatingPoint(3,1),ValidatingPoint(3,2),ValidatingPoint(3,3));
    [UG_4, UCM_4] = EvalUniform (Xmin_4,Ymin_4,Xmax_4,Ymax_4,XCM_4,YCM_4,Sig_4,ValidatingPoint(4,1),ValidatingPoint(4,2),ValidatingPoint(4,3));
    [UG_5, UCM_5] = EvalUniform (Xmin_5,Ymin_5,Xmax_5,Ymax_5,XCM_5,YCM_5,Sig_5,ValidatingPoint(5,1),ValidatingPoint(5,2),ValidatingPoint(5,3));
   % [UG_6, UCM_6] = EvalUniform (Xmin_6,Ymin_6,Xmax_6,Ymax_6,XCM_6,YCM_6,Sig_6,ValidatingPoint(6,1),ValidatingPoint(6,2),ValidatingPoint(6,3));
   % [UG_7, UCM_7] = EvalUniform (Xmin_7,Ymin_7,Xmax_7,Ymax_7,XCM_7,YCM_7,Sig_7,ValidatingPoint(7,1),ValidatingPoint(7,2),ValidatingPoint(7,3));
   % [UG_8, UCM_8] = EvalUniform (Xmin_8,Ymin_8,Xmax_8,Ymax_8,XCM_8,YCM_8,Sig_8,ValidatingPoint(8,1),ValidatingPoint(8,2),ValidatingPoint(8,3));
    SystemUniformity(numr,:)= [UG_1,UCM_1,UG_2,UCM_2,UG_3,UCM_3,UG_4,UCM_4,UG_5,UCM_5];%,UG_6,UCM_6,UG_7,UCM_7,UG_8,UCM_8];
    register=0;
    numr=numr+1;
    %redhot=redhot+1
end

register=register+(dT*t0);


%{


f=floor(log10(iter));
if f<6
if count<1
	count=count+(1/(10^f));
else
     
            
    [Xmin_1,Ymin_1,Xmax_1,Ymax_1,XCM_1,YCM_1,Sig_1]= EvalMorphPoint (ValidatingPoint(1,1),ValidatingPoint(1,2),ValidatingPoint(1,3),N);
    [Xmin_2,Ymin_2,Xmax_2,Ymax_2,XCM_2,YCM_2,Sig_2]= EvalMorphPoint (ValidatingPoint(2,1),ValidatingPoint(2,2),ValidatingPoint(2,3),N);
    [Xmin_3,Ymin_3,Xmax_3,Ymax_3,XCM_3,YCM_3,Sig_3]= EvalMorphPoint (ValidatingPoint(3,1),ValidatingPoint(3,2),ValidatingPoint(3,3),N);
    [Xmin_4,Ymin_4,Xmax_4,Ymax_4,XCM_4,YCM_4,Sig_4]= EvalMorphPoint (ValidatingPoint(4,1),ValidatingPoint(4,2),ValidatingPoint(4,3),N); 
    [Xmin_5,Ymin_5,Xmax_5,Ymax_5,XCM_5,YCM_5,Sig_5]= EvalMorphPoint (ValidatingPoint(5,1),ValidatingPoint(5,2),ValidatingPoint(5,3),N);     
   % [Xmin_6,Ymin_6,Xmax_6,Ymax_6,XCM_6,YCM_6,Sig_6]= EvalMorphPoint (ValidatingPoint(6,1),ValidatingPoint(6,2),ValidatingPoint(6,3),N);
   % [Xmin_7,Ymin_7,Xmax_7,Ymax_7,XCM_7,YCM_7,Sig_7]= EvalMorphPoint (ValidatingPoint(7,1),ValidatingPoint(7,2),ValidatingPoint(7,3),N);
   % [Xmin_8,Ymin_8,Xmax_8,Ymax_8,XCM_8,YCM_8,Sig_8]= EvalMorphPoint (ValidatingPoint(8,1),ValidatingPoint(8,2),ValidatingPoint(8,3),N);
    SystemRange(numr,:)= [Xmin_1,Ymin_1,Xmax_1,Ymax_1,XCM_1,YCM_1,Xmin_2,Ymin_2,Xmax_2,Ymax_2,XCM_2,YCM_2,Xmin_3,Ymin_3,Xmax_3,Ymax_3,XCM_3,YCM_3,Xmin_4,Ymin_4,Xmax_4,Ymax_4,XCM_4,YCM_4,Xmin_5,Ymin_5,Xmax_5,Ymax_5,XCM_5,YCM_5];%,Xmin_6,Ymin_6,Xmax_6,Ymax_6,XCM_6,YCM_6,Sig_6,Xmin_7,Ymin_7,Xmax_7,Ymax_7,XCM_7,YCM_7,Sig_7,Xmin_8,Ymin_8,Xmax_8,Ymax_8,XCM_8,YCM_8,Sig_8];    
    SystemVar(numr,:,:)=[Sig_1,Sig_2,Sig_3,Sig_4,Sig_5];
    [UG_1, UCM_1] = EvalUniform (Xmin_1,Ymin_1,Xmax_1,Ymax_1,XCM_1,YCM_1,Sig_1,ValidatingPoint(1,1),ValidatingPoint(1,2),ValidatingPoint(1,3));
    [UG_2, UCM_2] = EvalUniform (Xmin_2,Ymin_2,Xmax_2,Ymax_2,XCM_2,YCM_2,Sig_2,ValidatingPoint(2,1),ValidatingPoint(2,2),ValidatingPoint(2,3));
    [UG_3, UCM_3] = EvalUniform (Xmin_3,Ymin_3,Xmax_3,Ymax_3,XCM_3,YCM_3,Sig_3,ValidatingPoint(3,1),ValidatingPoint(3,2),ValidatingPoint(3,3));
    [UG_4, UCM_4] = EvalUniform (Xmin_4,Ymin_4,Xmax_4,Ymax_4,XCM_4,YCM_4,Sig_4,ValidatingPoint(4,1),ValidatingPoint(4,2),ValidatingPoint(4,3));
    [UG_5, UCM_5] = EvalUniform (Xmin_5,Ymin_5,Xmax_5,Ymax_5,XCM_5,YCM_5,Sig_5,ValidatingPoint(5,1),ValidatingPoint(5,2),ValidatingPoint(5,3));
   % [UG_6, UCM_6] = EvalUniform (Xmin_6,Ymin_6,Xmax_6,Ymax_6,XCM_6,YCM_6,Sig_6,ValidatingPoint(6,1),ValidatingPoint(6,2),ValidatingPoint(6,3));
   % [UG_7, UCM_7] = EvalUniform (Xmin_7,Ymin_7,Xmax_7,Ymax_7,XCM_7,YCM_7,Sig_7,ValidatingPoint(7,1),ValidatingPoint(7,2),ValidatingPoint(7,3));
   % [UG_8, UCM_8] = EvalUniform (Xmin_8,Ymin_8,Xmax_8,Ymax_8,XCM_8,YCM_8,Sig_8,ValidatingPoint(8,1),ValidatingPoint(8,2),ValidatingPoint(8,3));
    SystemUniformity(numr,:)= [UG_1,UCM_1,UG_2,UCM_2,UG_3,UCM_3,UG_4,UCM_4,UG_5,UCM_5];%,UG_6,UCM_6,UG_7,UCM_7,UG_8,UCM_8];

    count=0;
    numr=numr+1;
       
end

else
    if count<1
	count=count+(1/(10^(f-1)));
    else
        
    [Xmin_1,Ymin_1,Xmax_1,Ymax_1,XCM_1,YCM_1,Sig_1]= EvalMorphPoint (ValidatingPoint(1,1),ValidatingPoint(1,2),ValidatingPoint(1,3),N);
    [Xmin_2,Ymin_2,Xmax_2,Ymax_2,XCM_2,YCM_2,Sig_2]= EvalMorphPoint (ValidatingPoint(2,1),ValidatingPoint(2,2),ValidatingPoint(2,3),N);
    [Xmin_3,Ymin_3,Xmax_3,Ymax_3,XCM_3,YCM_3,Sig_3]= EvalMorphPoint (ValidatingPoint(3,1),ValidatingPoint(3,2),ValidatingPoint(3,3),N);
    [Xmin_4,Ymin_4,Xmax_4,Ymax_4,XCM_4,YCM_4,Sig_4]= EvalMorphPoint (ValidatingPoint(4,1),ValidatingPoint(4,2),ValidatingPoint(4,3),N); 
    [Xmin_5,Ymin_5,Xmax_5,Ymax_5,XCM_5,YCM_5,Sig_5]= EvalMorphPoint (ValidatingPoint(5,1),ValidatingPoint(5,2),ValidatingPoint(5,3),N);
   % [Xmin_6,Ymin_6,Xmax_6,Ymax_6,XCM_6,YCM_6,Sig_6]= EvalMorphPoint (ValidatingPoint(6,1),ValidatingPoint(6,2),ValidatingPoint(6,3),N);
   % [Xmin_7,Ymin_7,Xmax_7,Ymax_7,XCM_7,YCM_7,Sig_7]= EvalMorphPoint (ValidatingPoint(7,1),ValidatingPoint(7,2),ValidatingPoint(7,3),N);
   % [Xmin_8,Ymin_8,Xmax_8,Ymax_8,XCM_8,YCM_8,Sig_8]= EvalMorphPoint (ValidatingPoint(8,1),ValidatingPoint(8,2),ValidatingPoint(8,3),N);
    SystemRange(numr,:)= [Xmin_1,Ymin_1,Xmax_1,Ymax_1,XCM_1,YCM_1,Xmin_2,Ymin_2,Xmax_2,Ymax_2,XCM_2,YCM_2,Xmin_3,Ymin_3,Xmax_3,Ymax_3,XCM_3,YCM_3,Xmin_4,Ymin_4,Xmax_4,Ymax_4,XCM_4,YCM_4,Xmin_5,Ymin_5,Xmax_5,Ymax_5,XCM_5,YCM_5];
    SystemVar(numr,:,:)=[Sig_1,Sig_2,Sig_3,Sig_4,Sig_5];
    
    [UG_1, UCM_1] = EvalUniform (Xmin_1,Ymin_1,Xmax_1,Ymax_1,XCM_1,YCM_1,Sig_1,ValidatingPoint(1,1),ValidatingPoint(1,2),ValidatingPoint(1,3));
    [UG_2, UCM_2] = EvalUniform (Xmin_2,Ymin_2,Xmax_2,Ymax_2,XCM_2,YCM_2,Sig_2,ValidatingPoint(2,1),ValidatingPoint(2,2),ValidatingPoint(2,3));
    [UG_3, UCM_3] = EvalUniform (Xmin_3,Ymin_3,Xmax_3,Ymax_3,XCM_3,YCM_3,Sig_3,ValidatingPoint(3,1),ValidatingPoint(3,2),ValidatingPoint(3,3));
    [UG_4, UCM_4] = EvalUniform (Xmin_4,Ymin_4,Xmax_4,Ymax_4,XCM_4,YCM_4,Sig_4,ValidatingPoint(4,1),ValidatingPoint(4,2),ValidatingPoint(4,3));
    [UG_5, UCM_5] = EvalUniform (Xmin_5,Ymin_5,Xmax_5,Ymax_5,XCM_5,YCM_5,Sig_5,ValidatingPoint(5,1),ValidatingPoint(5,2),ValidatingPoint(5,3));
   % [UG_6, UCM_6] = EvalUniform (Xmin_6,Ymin_6,Xmax_6,Ymax_6,XCM_6,YCM_6,Sig_6,ValidatingPoint(6,1),ValidatingPoint(6,2),ValidatingPoint(6,3));
   % [UG_7, UCM_7] = EvalUniform (Xmin_7,Ymin_7,Xmax_7,Ymax_7,XCM_7,YCM_7,Sig_7,ValidatingPoint(7,1),ValidatingPoint(7,2),ValidatingPoint(7,3));
   % [UG_8, UCM_8] = EvalUniform (Xmin_8,Ymin_8,Xmax_8,Ymax_8,XCM_8,YCM_8,Sig_8,ValidatingPoint(8,1),ValidatingPoint(8,2),ValidatingPoint(8,3));
    SystemUniformity(numr,:)= [UG_1,UCM_1,UG_2,UCM_2,UG_3,UCM_3,UG_4,UCM_4,UG_5,UCM_5];%UG_6,UCM_6,UG_7,UCM_7,UG_8,UCM_8];        
    
    count=0;
    numr=numr+1;
end
    
end
%}




end


timepicture=timepicture+(dT*t0);

if timepicture > 120
    i=size(PictureE,1);
    PictureE(i+1,:,:,:)=E; 
    timepicture=0;
end


%Registering copies of the system at given iterations
%{
if iter ==500
E500=E;
t500=time;
end

if iter ==1000
E1k=E;
t1k=time;
end
 
if iter ==5000
E5k=E;
t5k=time;
end 

if iter ==10000
E10k=E;
t10k=time;
end

if iter ==15000
E15k=E;
t15k=time;
end


if iter ==20000
E20k=E;
t20k=time;
end

if iter ==50000
E50k=E;
t50k=time;
end 
 
if iter ==100000
E100k=E;
t100k=time;
end
 
 
if iter ==200000
E200k=E;
t200k=time;
end


if iter ==500000
E500k=E;
t500k=time;
end


if iter ==1000000
E1M=E;
t1M=time;
end

if iter ==2500000
E2_5M=E;
t2_5M=time;
end

if iter ==5000000
E5M=E;
t5M=time;
end
%}
  


iter=iter+1;
end




switch (SystemStyle)
    case 1
%Writing Vector of parameters for Fusing Spheres
    xlswrite('Parameters_Fusion_Y=1.1_Tao=1_Dini=0.5_v4.xls',Parameters);
    xlswrite('Results_Fusion_Y=1.1_Tao=1_Dini=0.5_v4.xls',Results);
    save ('Y=0.5_Tao=1_Dini=0.5_v4.mat');
%writematrix(Parameters);
%writematrix(Results);
%}
    otherwise
%Writing Vector of parameters for Validating System
    xlswrite('SystemRange_Case01_Tao=1_v1.xls',SystemRange);
    xlswrite('Uniformity_Case01_Tao=1_v1.xls',SystemUniformity);
    save ('Y=0.5_Tao=1_Case01_v1.mat');
end
toc 
 

    
    
   







