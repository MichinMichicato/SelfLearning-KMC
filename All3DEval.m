function [CMx,CMy,CMz,CMxR,CMyR, CMzR,CMxL, CMyL, CMzL]=All3DEval(N)
%Evalua los centros de masa global y del lado izquierdo y derecho de un
%espacio 3D
global E;
[CMxL,CMyL,CMzL]=LeftSide(N);
[CMxR,CMyR,CMzR]=RightSide(N);
[CMx,CMy,CMz]=Whole(N);

end

function [CMxL,CMyL,CMzL] = LeftSide (N)
global E;
Rx=0;
Ry=0;
Rz=0;
M=0;
for i=1:N
for j=1:N/2
for z=1:N    
  
        if E(i,j,z)==1
        M=M+1;
        Rx=Rx+j;
        Ry=Ry+i;
        Rz=Rz+z;
        end
 
end
end
end
CMxL=Rx/M;
CMyL=Ry/M;
CMzL=Rz/M;


end
function [CMxR,CMyR,CMzR] = RightSide (N)
global E;
Rx=0;
Ry=0;
Rz=0;
M=0;
for i=1:N
for j=N/2:N
for z=1:N
        if E(i,j,z)==1
        M=M+1;
        Rx=Rx+j;
        Ry=Ry+i;
        Rz=Rz+z;
        end
        
end    
end
end
CMxR=Rx/M;
CMyR=Ry/M;
CMzR=Rz/M;

end

function [CMx,CMy,CMz]=Whole(N)
global E;
Rx=0;
Ry=0;
Rz=0;
M=0;
for i=1:N
for j=1:N
for z=1:N
        if E(i,j,z)==1
        M=M+1;
        Rx=Rx+j;
        Ry=Ry+i;
        Rz=Rz+z;
        end
    
end
end   
end
CMx=Rx/M;
CMy=Ry/M;
CMz=Rz/M;
end