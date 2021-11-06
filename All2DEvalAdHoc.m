function [CMz,CMy]=All2DEvalAdHoc(N,E,Xinter)
%Evalua los centros de masa global y del lado izquierdo y derecho de un

%XInter is the plane on X (over cannonic axes) on which All 2DEvalAdHoc Evaluates
%the system (system's slice)

%espacio 2D
%
%[CMxL,CMyL]=LeftSide(N);
%[CMxR,CMyR]=RightSide(N);
[CMz,CMy]=Whole(N,E,Xinter);
end
%{
function [CMxL,CMyL] = LeftSide (N)
global E;
Rx=0;
Ry=0;
M=0;
for i=1:N
for j=1:N/2
        if E(i,j)==1
            M=M+1;
        Rx=Rx+j;
        Ry=Ry+i;
        end
    
   end
end
CMxL=Rx/M;
CMyL=Ry/M;
end
function [CMxR,CMyR] = RightSide (N)
global E;
Rx=0;
Ry=0;
M=0;
for i=1:N
for j=N/2:N
        if E(i,j)==1
            M=M+1;
        Rx=Rx+j;
        Ry=Ry+i;
        end
    
   end
end
CMxR=Rx/M;
CMyR=Ry/M;

end
%}
function [CMz,CMy]=Whole(N,E,Xinter)

Rz=0;
Ry=0;
M=0;
for i=1:N
for z=1:N
        if E(i,Xinter,z)==1
        M=M+1;
        Rz=Rz+z;
        Ry=Ry+i;
        end
    
   end
end
CMz=Rz/M;
CMy=Ry/M;

end