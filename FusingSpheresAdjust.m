function [OptRadR,OptRadL,Rneck,RneckB,RneckC,RneckD,TetaA,TetaB,TetaC,Xinter,DistCM,VoidR,VoidL] = FusingSpheresAdjust(N)

global E;

OptRadR=0;
OptRadL=0;


%This approach was not used, yet was left in notes to help the reader to
%think about alternative implementations of the current script.

%%First we're gonna read the data from a space named E (or EINI) provided
%%by any other script routine
%%Then we're gonna run the routine to adjust the geometries of two spehres
%%to the system so that we can find the geometrical parameters analyzed
%%when calibrating the model using the analogy with two fusing droplets

%%ESP = zeros(:,:,:);
%fid = fopen('Test_For_Reading_Flat_File_Into_Matrix.txt', 'r');
%ESP = textscan(fid, '%d');


%%COPIED FROM KMC_3D_FUSING_SPHERES_TESIS031019
%-------------------------------------------------------------------------    


%N tamaño del espacio, tsimulacion tiempo max a simular
%N=50;

%Creando Matriz del espacio

%E = int32(zeros(N,N,N));

%Cuadrado de los Radios de los N agregados
%Rad1=100;
%Rad2=100;

%Ubicación Centros de las N esferas (agregados)
%{
Cx1 = int32(N/2)-(Rad1^0.5)-1;
Cy1 = int32(N/2) ;
Cz1 = int32((N/2));
Cx2 = int32(N/2)+(Rad2^0.5)+1;
Cy2 = int32((N/2)) ;
Cz2 = int32((N/2));
%}

%creando sistema (ubicación de células)
%{
for i=1:N
    for o = 1:N
        for z=1:N
            ru=rand;
            if ((i-Cy1)^2+(o-Cx1)^2+(z-Cz1)^2)<= Rad1 || ((i-Cy2)^2+(o-Cx2)^2+(z-Cz2)^2)<= Rad2
             if ru>0.7
                    E(i,o,z)=1;
             end
            end
        end
    end
end
%}

%-------------------------------------------------------------------------    

%Call the center of mass coordinates for the system, either Whole (X,Y,Z),
%and for the right and left part (i.e. divinding on the X axis in halfs)
%returning XR,YR,LR and XL,YL,ZL respectively.


[X,Y,Z,XR,YR,ZR,XL,YL,ZL] = All3DEval(N);


%Initial Guess for Radius Search (not necessary for fminbnd
%RadRguess=(N/2);

%In addition to fminbnd, one can use: fminsearch, fmin.

%_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

%Notice the variable MaxPosR has as objective reduce the time of evaluation
%every time the function preCond is called (reduces the domain over which
%it is evaluated)

%MaxR MinR (see the Feasilibility evaluating module) may hamper the
%performance of the function for very irregular distributions where the
%derivative (of Error_R_R function) may change of sign more than once before reaching the global
%minimum


%Example of how to use PreCond
MaxPosR=ceil(0.87*N);

%{


global SizeRangeEval;
global RangeEval;
Division =N/200;
SizeRangeEval=ceil(MaxPosR / Division);
RangeEval = linspace (0,MaxPosR,SizeRangeEval);

[RANGR,RAGR]=PreCond(XR,YR,ZR,N,1);
[RANGL,RAGL]=PreCond(XL,YL,ZL,N,0);

negR=0;
posR=0;
negL=0;
posL=0;

%minR=0;
%maxR=0;
%minL=0;
%maxL=0;


%Module Evaluating the feasibility of the sphere adjustment

%Right
ind=1;
while negR ==0 && ind<=size(RANGR,2)
 if (RANGR(ind)<0)
    % minR=ind;
     negR=1;
 end
 ind=ind+1;
end
while posR ==0 && ind<=size(RANGR,2)
 if (RANGR(ind)>0)
    % maxR=ind+1;
     posR=1;
 end
 ind=ind+1;
end

%Left
ind=1;
while negL ==0 && ind<=size(RANGL,2)
 if (RANGL(ind)<0)
    % minL=ind;
     negL=1;
 end
 ind=ind+1;
end
while posL ==0 && ind<=size(RANGL,2)
 if (RANGL(ind)>0)
    % maxL=ind+1;
     posL=1;
 end
 ind=ind+1;
end
%}

%Radius of Sphere fitting

%%Note: In case of non homogeneous spheres the minR and maxR may cause a
%%false positive for convergence as Error_R_R may not have a soft behavior
%%and the value of derivatives change of sign many times before function
%%reaches its minimum

%TolX is intended to adjuts the tolerance under which the function fminbnd
%can be stopped as it is considered that it has converged

%{
%This part was used at first without the use of the precond function, in
case of instability abilitate this and disactivate the conditional based on
precond function

OptRadR = fminbnd(@(RadR) Error_R_R (XR,YR,ZR,N,RadR,1),minL,MaxPosR);%,optimset('TolX',0.00025));
OptRadL = fminbnd(@(RadL) Error_R_R (XL,YL,ZL,N,RadL,0),minL,MaxPosR);%,optimset('TolX',0.00025));
disp (OptRadR);
disp (OptRadL);
%}


%if (negR==1 && posR==1)&& (negL==1 && posL==1)
    
    OptRadR = fminbnd(@(RadR) Error_R_R (XR,YR,ZR,N,RadR,1),0,MaxPosR,optimset('TolX',0.1));
    OptRadL = fminbnd(@(RadL) Error_R_R (XL,YL,ZL,N,RadL,0),0,MaxPosR,optimset('TolX',0.1));
 %{
    %OptRad1=  fmincon(@(RadR) Error_R_R (XR,YR,ZR,N,RadR),RadRguess,[],[],[],[],0,(N^2));
    %disp (OptRadR);
    %disp (OptRadL);
    %}
%else
%    if (negL==1 && posL==1)
%        disp('The system does not adjust an spehere on the Right - check manually');
%    else
%        disp('The system does not adjust an spehere on the Left - check manually');
%    end
%    OptRadR=0;
%    OptRadL=0;
%end

%{
if (OptRadR==0 || OptRadL==0)
    
    DistCM=(((XR-XL)^2)+((YR-YL)^2)+((ZR-ZL)^2))^0.5;
    Rneck=0;
    Xinter=1;
    TetaA=0;
    TetaB=0;
    
else
%}   

    DistCM = (((XR-XL)^2)+((YR-YL)^2)+((ZR-ZL)^2))^0.5;
    a1= 2*(XR-XL);
    da=a1/2;
    Sa=(XR+XL);
   % b1=2*(YR-YL);
   % db=b1/2;
    %Sb =(YR+YL);
   % c1= 2*(ZR-ZL);
    %dc=c1/2;
    %Sc=(ZR+ZL);
    D1=Sa*da+((OptRadL)^2)-((OptRadR)^2);%+Sb*db+Sc*dc
  %{
    %coeff1= -2*(b1/a1)*(XL);
    %coeff2= -2*(c1/a1)*(XL);
    %coeff3= 2*(b1*D1/(a1^2));
    %coeff4= 2*(c1*D1/(a1^2));
    %coeff5= -(2*(b1*c1/(a1^2)));
    %coeff6= -((b1/a1)^2);
    %coeff7= -((c1/a1)^2);
    %coeffTotal= coeff1+coeff2+coeff3+coeff4+coeff5+coeff6+coeff7;
    %}
    
    if OptRadR+OptRadL<DistCM
        Rneck=0;
        Xinter=1;
        TetaA=0;
        TetaB=0;
    else
        Rneck=((OptRadL^2)-(XL^2)+(2*D1*(XL)/a1)-((D1/a1)^2))^0.5;
        Xinter=D1/a1;
        TetaA=(360/(2*pi()))*acos((DistCM/2)/OptRadL);
        TetaB=(360/(2*pi()))*asin((Rneck)/OptRadL);
    end
    
        
%// - Alternative computations for the neck Radius ('RneckX')
    
Xintersec=int64(Xinter);
[CMzinter,CMyinter]=All2DEvalAdHoc(N,E,Xintersec);
CMzinter;
CMyinter;
%Rneck as resulting from best fitting circle at X=Xinter plane taking as
%center the CM for that plane
    RneckB = fminbnd(@(RneckB) Error_R_Plane (CMyinter,CMzinter,Xintersec,N,RneckB),0,N,optimset('TolX',0.1));
       
%Rneck as resulting from best fitting circle at X=Xinter plane taking as
%center the average of CM of both spheres
CMyinter2= ((YR+YL)/2);
CMzinter2= ((ZR+ZL)/2);
    RneckC = fminbnd(@(RneckC) Error_R_Plane (CMyinter2,CMzinter2,Xintersec,N,RneckC),0,N,optimset('TolX',0.1));
       
%Rneck as resulting from best fitting circle at the plane in the Mid-way of the vector connecting the CM of both spheres 
CMxinter3= int64((XR+XL)/2);
CMyinter3= ((YR+YL)/2);
CMzinter3= ((ZR+ZL)/2);
CMxinter3;
CMyinter3;
CMzinter3;

    RneckD = fminbnd(@(RneckD) Error_R_Plane (CMyinter3,CMzinter3,CMxinter3,N,RneckD),0,N,optimset('TolX',0.1));

%Alternatice Calculation of Tetca
RneckAvg=(RneckB+RneckC+RneckD)/3;
RadSphereAvg=(OptRadL+OptRadR)/2;
TetaC=(360/(2*pi()))*asin((RneckAvg)/RadSphereAvg);
    

%%-_--_--_--_--___--_---_--_-_-__-_-_-_-_--_-_-____---_--_--_--_--___--_---
%CALCULATING THE VOID FRACTION FOR EACH SPHERE IN THEIR RESPECTIVE HALF OF
%THE SPACE

%RIGHT

CelTotR =0;
CelINR=0;

for i=floor(YR-OptRadR):ceil(YR+OptRadR)
   for o = floor(XR-OptRadR): ceil(XR+OptRadR) 
       for z = floor(ZR-OptRadR): ceil(ZR+OptRadR)
           if (((i-YR)^2)+((o-XR)^2)+((z-ZR)^2)<=(OptRadR^2))&&(o>=(N/2))&& (i>=1 && i<=N) && (o>=1 && o<=N)&& (z>=1 && z<=N)
                CelTotR=CelTotR+1;
            if E(i,o,z)==1
                CelINR=CelINR+1;
            end
           end
           
       end
   end
end


%LEFT

CelTotL=0;
CelINL=0;

for i=floor(YL-OptRadL):ceil(YL+OptRadL)
   for o = floor(XL-OptRadL): ceil(XL+OptRadL) 
       for z = floor(ZL-OptRadL): ceil(ZL+OptRadL)
           if (((i-YL)^2)+((o-XL)^2)+((z-ZL)^2)<=(OptRadL^2))&&(o<=(N/2))&& (i>=1 && i<=N) && (o>=1 && o<=N)&& (z>=1 && z<=N)
            CelTotL=CelTotL+1;
            if E(i,o,z)==1
                CelINL=CelINL+1;
            end
           end
           
       end
   end
end

%Void Fraction calculation

VoidL=1-(CelINL/CelTotL);
%disp (VoidL);
VoidR=1-(CelINR/CelTotR);
%disp (VoidR);

%{
%THIS CYCLE ONLY MAKES SENSE WHEN THE SYSTEM TO BE EVALUATED IS RANDOMLY
GENERATED EACH TIME THE EVALUATION CYCLE IS RUN (I.E. AS IN THE INITIAL
PROGRAMMING PHASE OF THE SCRIPT) BUT IT DOES NOT MAKES SENSE WHEN IT IS
RECEIVED AS A FIXED INPUT SINCE SEVERAL EVALUATIONS LEAD TO THE SAME
ESTIMATED RADIUS

if (VoidL>0.5 || VoidR>0.5) ||(CelTotR==0 || CelTotL==0)
      
    if counter <=5
       
        GrabaRR(counter)=OptRadR;
        GrabaRL(counter)=OptRadL;
        RadAproxR=RadAproxR+OptRadR;
       RadAproxL=RadAproxL+OptRadL;
       counter=counter+1;
    else
        RadAproxR=RadAproxR/counter;
        RadAproxL=RadAproxL/counter;
        OptRadR=RadAproxR;
        OptRadL=RadAproxL;
        
       convergence=1;
       
       %Calculation of the fusing sphere parametric parameters in cas of
       %aprox average radius (high void fraction)
        DistCM = (((XR-XL)^2)+((YR-YL)^2)+((ZR-ZL)^2))^0.5;
        a1= 2*(XR-XL);
        da=a1/2;
        Sa=(XR+XL);
        b1=2*(YR-YL);
        db=b1/2;
        Sb =(YR+YL);
        c1= 2*(ZR-ZL);
        dc=c1/2;
        Sc=(ZR+ZL);
        D1=Sa*da+Sb*db+Sc*dc+((RadAproxL)^2)-((RadAproxR)^2);
        Rneck=((RadAproxL^2)-(XL^2)+(2*D1*(XL)/a1)-((D1/a1)^2))^0.5;
        Xinter=D1/a1;
        TetaA=(360/(2*pi()))*acos((DistCM/2)/RadAproxL);
        TetaB=(360/(2*pi()))*asin((Rneck)/RadAproxL);
       
    end
else
    convergence=1;
end
%}

%-_-_-_-_--_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_--_-_-_-_-_-_-_

%FUNCTIONS TO VISUALIZE STABILITY AND THE SYSTEM

%Options to monitor the behavior of Error_R_R given the System conditions
%fun=@(RadL) Error_R_R (XL,YL,ZL,N,RadL,0);
%fplot(@(L)fun (L),[0,12]);
%fun=@(RadR) Error_R_R (XR,YR,ZR,N,RadR,1);
%fplot(@(R)fun (R),[-0.5,10]);
%Using the comman "isosurface([1:N],[1:N],[1:N],E,0.5)" will generate a 3D
%plot of the Space for a "interpolating surcafe among i=1 and i=0.

end


%-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

%Suport Sub-Functions


function [RANG,RAG] = PreCond(X,Y,Z,N,side)
%global MaxPosR;
global SizeRangeEval;
global RangeEval;
VAL=double.empty;
for k=1: SizeRangeEval
   VAL(k)=Error_R_R(X,Y,Z,N,RangeEval(k),side);
end
RAG=VAL;
RANG=diff(VAL);
end


%On June 22nd a modification to Error_R_R was made, here a heavier (TE=3)
%contribution to the error when the spin is 0 (medium) and is within the
%spehere is set. In contrast spin 1 (cells) out of the sphere contribute to
%error in the same magnitude as before (TE=1).

function [TER] = Error_R_R(X,Y,Z,N,R,side)
TE=0;
global E;

if side == 0

    for i=1:N
        for o=1:N/2
            for z=1:N
            
            if (E(i,o,z)==1 && ((i-Y)^2+(o-X)^2+(z-Z)^2)> (R^2))||(E(i,o,z)==0 && ((i-Y)^2+(o-X)^2+(z-Z)^2)<= (R^2))
                if E(i,o,z)==0
                    TE=TE+1;
                else %E(i,o,z)==1
                    TE=TE+1;
                end
            end
                
            end
            
        end
    end

    TER=TE;

end

if side == 1
    
    for i=1:N
        for o=N/2:N
            for z=1:N
            
            if (E(i,o,z)==1 && ((i-Y)^2+(o-X)^2+(z-Z)^2)> (R^2))||(E(i,o,z)==0 && ((i-Y)^2+(o-X)^2+(z-Z)^2)<= (R^2))
                
                if E(i,o,z)==0
                    TE=TE+1;
                else %E(i,o,z)==1
                    TE=TE+1;
                end
                
            end
                
            end
            
        end
    end

    TER=TE;
end
end


function [TER] = Error_R_Plane(Y,Z,o,N,R)
TE=0;
global E;
%o is the Xinter (x at which both speheres intersect)

for i=1:N      
  for z=1:N
            if (E(i,o,z)==1 && ((i-Y)^2+(z-Z)^2)> (R^2))||(E(i,o,z)==0 && ((i-Y)^2+(z-Z)^2)<= (R^2))
                if E(i,o,z)==0
                    TE=TE+0.5;
                else %E(i,o,z)==1
                    TE=TE+1;
                end
            end
  end
end

    TER=TE;

end