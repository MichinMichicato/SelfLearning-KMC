%Algorithm Fortran 95 - Last Version to New Version Matlab
%Revisar Algoritmo de creacion de tamaño de lista y llenado de eventos
%cuando hay celulas acercandose a las fronteras.
%-/-/-/USES FusingSpheresAdjust, All3DEval, All2DEvalAdHoc /-/-/-/-/-/-/-/
tic
global E;


%N tamaño del espacio, tsimulacion tiempo max a simular
N=60; %N shall be at least equal to sum of radiuses +4 (due to frontier clearance condition)
tsimulacion=7000;

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

%Parametros Principales
ta0=double(1);
Ycm=double(1.0);

%Declaracion de otras variables
time=double(0);
treal = double(0);
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
Parameters=double.empty;
Results=double.empty;





%creando sistema (ubicación de células)
for i=1:N
    for o = 1:N
        for z=1:N
            if ((i-Cy1)^2+(o-Cx1)^2+(z-Cz1)^2)<= Rad1 || ((i-Cy2)^2+(o-Cx2)^2+(z-Cz2)^2)<= Rad2
                E(i,o,z)=1;
            end
        end
    end
end
%Generando Sistema

%Registro sistema inicial en una matriz dentro del programa (para acceso
%rapido)

EINI=E;

%Registrando el espacio (sistema) en un archivo txt
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
 treal=(ta0*time); %Conversion to real time equivalent (sec) of current time
 iter=iter+1;


%-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


%Algoritmo Principal
while iter<=10000000 %treal < tsimulacion

    
 
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
 dT=-(log(E2)/Rtot); %Time leap for the event (in ta0 units)
 time=time+dT; %Current time expressed in units of ta0 (charact time)
 treal=(ta0*time); %Conversion to real time equivalent (sec) of current time

 
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


%Registrar CM (Hay que definir cada cuantos pasos se va a registrar)
%Podria omitir el conteo de M de nuevo y basarme solo en numcel que fueron
%contadas anteriormente.

%Grabar espacio cada intervalo fijo
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

%RE- ACTIVAR PARA REGISTRAR PARAMETROS DE SIMULACION



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
	count=count+(1/(5^(f-1)));
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

%fclose(fileID);    
%{
    [OptRadR,OptRadL,Rneck,RneckB,RneckC,RneckD,TetaA,TetaB,TetaC,Xinter,DistCM,VoidR,VoidL] = FusingSpheresAdjust(N);
	Parameters(numr,:)= [iter,dT,time,treal,OptRadR,OptRadL,Rneck,RneckB,RneckC,RneckD,TetaA,TetaB,TetaC,DistCM];
    Results(numr,:)=[iter,dT,time,treal,OptRadR,OptRadL,Rneck,RneckB,RneckC,DistCM];
    numr=numr+1;

%}

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

if iter ==7500000
E7_5M=E;
t7_5M=time;
end

if iter ==10000000
E10M=E;
t10M=time;
end
    
iter=iter+1;
end

xlswrite('Parameters_Fusion_Y=1.0_Tao=1_Dini=0.5_v5_10Mn.xls',Parameters);
xlswrite('Results_Fusion_Y=1.0_Tao=1_Dini=0.5_v5_10Mn.xls',Results);
save ('Y=1.0_Tao=1_Dini=0.5_v5_10Mn.mat');
%writematrix(Parameters);
%writematrix(Results);







toc 
 

    
    
   







