function [UG, UCM] = EvalUniform (Xmin,Ymin,Xmax,Ymax,XCM,YCM,Sig,i,j,z)
%The min and max and CM on each direction are found in EvaMorphPoint or Subroutine_Indexes_140619
global E;
UG=[];
UCM=[];
k=2.5 %About 85% of population in a normal distrib (See Z-score normal distrib tables)
UG(1)=evalung1(Ymin(1),Ymax(1),i,j,z);
UG(2)=evalung2(Xmin(2),Ymin(2),Xmax(2),Ymax(2),i,j,z);
UG(3)=evalung3(Xmin(3),Xmax(3),i,j,z);
UG(4)=evalung4(Xmin(4),Ymin(4),Xmax(4),Ymax(4),i,j,z);
UG(5)=evalung5(Ymin(5),Ymax(5),i,j,z);
UG(6)=evalung6(Xmin(6),Ymin(6),Xmax(6),Ymax(6),i,j,z);
UG(7)=evalung7(Xmin(7),Xmax(7),i,j,z);
UG(8)=evalung8(Xmin(8),Ymin(8),Xmax(8),Ymax(8),i,j,z);
UG(9)=(sum(UG(1:8)))/8;

UCM(1)=evaluncm1 (YCM(1),Sig,i,j,z,k);
UCM(2)=evaluncm2 (XCM(2),YCM(2),Sig,i,j,z,k);
UCM(3)=evaluncm3 (XCM(3),Sig,i,j,z,k);
UCM(4)=evaluncm4 (XCM(4),YCM(4),Sig,i,j,z,k);
UCM(5)=evaluncm5 (YCM(5),Sig,i,j,z,k);
UCM(6)=evaluncm6 (XCM(6),YCM(6),Sig,i,j,z,k);
UCM(7)=evaluncm7 (XCM(7),Sig,i,j,z,k);
UCM(8)=evaluncm8 (XCM(8),YCM(8),Sig,i,j,z,k);
UCM(9)=sum(UCM(1:8))/8;

end

%Calculates Global - Uniformity indexes (From closest to furthest)

function [Ux1]= evalung1 (Ymin,Ymax1,i,j,z)
global E;
M1=0;
i=Ymin;
Ux1=0;

if i == 0
Ux1=0;
else
while i>=Ymax1
   for h=-1:1
       for o =-1:1
           M1=M1+1;
           if E(i,j+h,z+o)==1
              Ux1=Ux1+1;
           end
    
   end
   end
   i=i-1;
end
Ux1=Ux1/M1;
end
end

function [Ux2] = evalung2 (Xmin,Ymin,Xmax,Ymax,i,j,z)
global E;
M2=0;
i=Ymin;
j=Xmin;
Ux2=0;

if i == 0 || j==0
Ux2=0;
else
while i>=Ymax && j<=Xmax
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
              M2=M2+1;
              if E(i,j,z+o)==1
              Ux2=Ux2+1;        
              end
            end
        case 2
            for o=-1:1
              M2=M2+1;
              if E(i,j+1,z+o)==1
              Ux2=Ux2+1;        
              end
            end
        case 3
            
            for o=-1:1
              M2=M2+1;
              if E(i-1,j,z+o)==1
              Ux2=Ux2+1;        
              end
            end
    end
    end
   i=i-1;
   j=j+1;
end
Ux2=Ux2/M2;
end
end

function [Ux3] = evalung3 (Xmin, Xmax,i,j,z)
global E;
M3=0;
j=Xmin;
Ux3=0;

if j == 0
Ux3=0;
else
while j<=Xmax
   for q=-1:1
       for o =-1:1
            M3=M3+1;
           if E(i+q,j,z+o)==1  
              Ux3=Ux3+1;
           end
           
       end
   end
   j=j+1;
end 
Ux3=Ux3/M3;
end
end

function [Ux4] = evalung4 (Xmin,Ymin,Xmax,Ymax,i,j,z)
global E;
M4=0;
i=Ymin;
j=Xmin;
Ux4=0;

if i == 0 || j==0
Ux4=0;
else

while i<=Ymax && j<=Ymax
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
                M4=M4+1;
              if E(i,j,z+o)==1
              Ux4=Ux4+1;        
              end
            end
        case 2
            for o=-1:1
                M4=M4+1;
              if E(i,j+1,z+o)==1
              Ux4=Ux4+1;        
              end
            end
        case 3
            
            for o=-1:1
                M4=M4+1;
              if E(i+1,j,z+o)==1
              Ux4=Ux4+1;        
              end
            end
     end
   end
   i=i+1;
   j=j+1;
end
Ux4=Ux4/M4;
end
end

function [Ux5] = evalung5 (Ymin,Ymax,i,j,z)
global E;
M5=0;
i=Ymin;
Ux5=0;

if i == 0
Ux5=0;
else
while i<=Ymax
   for h=-1:1
       for o =-1:1
            M5=M5+1;
           if E(i,j+h,z+o)==1
              Ux5=Ux5+1;
           end
           
       end
   end
   i=i+1;
end 
Ux5=Ux5/M5;
end
end

function [Ux6] = evalung6 (Xmin,Ymin,Xmax,Ymax,i,j,z)
global E;
M6=0;
i=Ymin;
j=Xmin;
Ux6=0;

if i == 0 || j==0
Ux6=0;
else
while i<=Ymax && j>=Xmax
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
                M6=M6+1;
              if E(i,j,z+o)==1
              Ux6=Ux6+1;        
              end
            end
        case 2
            for o=-1:1
                M6=M6+1;
              if E(i,j-1,z+o)==1
              Ux6=Ux6+1;        
              end
            end
        case 3
            
            for o=-1:1
                M6=M6+1;
              if E(i+1,j,z+o)==1
              Ux6=Ux6+1;        
              end
            end
       end
   end
   i=i+1;
   j=j-1;
end
Ux6=Ux6/M6;
end
end

function [Ux7] = evalung7 (Xmin, Xmax,i,j,z)
global E;
M7=0;
j=Xmin;
Ux7=0;

if j == 0
Ux7=0;
else
while j>=Xmax
   for q=-1:1
       for o =-1:1
            M7=M7+1;
            if E(i+q,j,z+o)==1
              Ux7=Ux7+1;
            end
           
       end
   end
   j=j-1;
end
Ux7=Ux7/M7;
%Ux7
end
end

function [Ux8] = evalung8 (Xmin,Ymin,Xmax,Ymax,i,j,z)
global E;
M8=0;
i=Ymin;
j=Ymin;
Ux8=0;

if i == 0 || j==0
Ux8=0;
else
while i>=Ymax && j>=Xmax
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
                M8=M8+1;
              if E(i,j,z+o)==1
              Ux8=Ux8+1;        
              end
            end
        case 2
            for o=-1:1
                M8=M8+1;
              if E(i,j-1,z+o)==1
              Ux8=Ux8+1;        
              end
            end
        case 3
            for o=-1:1
                M8=M8+1;
              if E(i-1,j,z+o)==1
              Ux8=Ux8+1;        
              end
            end
       end
   end
   i=i-1;
   j=j-1;
end
Ux8=Ux8/M8;
end
end

%Calculates de CM - Uniformity indexes (from CM +- k std deviations)
%Sig : Std Deviations -> vector of 8 (radial directions) lines and two
%columns (for X and Y coordinates)

function [Ucm1]= evaluncm1 (YCM,Sig,I,j,z,k)
global E;
M1=0;
i=ceil(YCM+k*Sig(1,2));
Ucm1=0;

if i==0
  Ucm1=0; 
    
else
while i>=floor(YCM-k*Sig(1,2))
   % YCM
    %k
    %Sig(1,2)
   for h=-1:1
       for o =-1:1
           M1=M1+1;
          % i
         %  j+h
          % z+o
           if E(i,j+h,z+o)==1
              Ucm1=Ucm1+1;
           end
    
   end
   end
   i=i-1;
end
Ucm1=Ucm1/M1;

end


end

function [Ucm2] = evaluncm2 (XCM,YCM,Sig,I,J,z,k)
global E;
M2=0;
i=ceil(YCM+k*Sig(2,2));
j=floor(XCM-k*Sig(2,1));
Ucm2=0;

if i==0 || j==0
  Ucm2=0; 
else
while i>=floor(YCM-k*Sig(2,2)) && j<=ceil(XCM+k*Sig(2,1))
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
              M2=M2+1;
              if E(i,j,z+o)==1
              Ucm2=Ucm2+1;        
              end
            end
        case 2
            for o=-1:1
              M2=M2+1;
              if E(i,j+1,z+o)==1
              Ucm2=Ucm2+1;        
              end
            end
        case 3
            
            for o=-1:1
              M2=M2+1;
              if E(i-1,j,z+o)==1
              Ucm2=Ucm2+1;        
              end
            end
    end
    end
   i=i-1;
   j=j+1;
end
Ucm2=Ucm2/M2;
end

end

function [Ucm3] = evaluncm3 (XCM,Sig,i,J,z,k)
global E;
M3=0;
j=floor(XCM-k*Sig(3,1));
Ucm3=0;
if j==0
  Ucm3=0;
else
while j<=ceil(XCM+k*Sig(3,1))
   for q=-1:1
       for o =-1:1
            M3=M3+1;
           if E(i+q,j,z+o)==1  
              Ucm3=Ucm3+1;
           end
           
       end
   end
   j=j+1;
end 
Ucm3=Ucm3/M3;
end
end

function [Ucm4] = evaluncm4 (XCM,YCM,Sig,I,J,z,k)
global E;
M4=0;
i=floor(YCM-k*Sig(4,2));
j=floor(XCM-k*Sig(4,1));
Ucm4=0;

if i==0 || j==0
  Ucm4=0; 
else

while i<=ceil(YCM+k*Sig(4,2)) && j<=ceil(YCM+k*Sig(4,1))
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
                M4=M4+1;
              if E(i,j,z+o)==1
              Ucm4=Ucm4+1;        
              end
            end
        case 2
            for o=-1:1
                M4=M4+1;
              if E(i,j+1,z+o)==1
              Ucm4=Ucm4+1;        
              end
            end
        case 3
            
            for o=-1:1
                M4=M4+1;
              if E(i+1,j,z+o)==1
              Ucm4=Ucm4+1;        
              end
            end
     end
   end
   i=i+1;
   j=j+1;
end
Ucm4=Ucm4/M4;
end
end

function [Ucm5] = evaluncm5 (YCM,Sig,I,j,z,k)
global E;
M5=0;
i=floor(YCM-k*Sig(5,2));
Ucm5=0;

if i==0
  Ucm5=0; 
else
while i<=ceil(YCM+k*Sig(5,2))
   for h=-1:1
       for o =-1:1
            M5=M5+1;
           if E(i,j+h,z+o)==1
              Ucm5=Ucm5+1;
           end
           
       end
   end
   i=i+1;
end 
Ucm5=Ucm5/M5;
end
end

function [Ucm6] = evaluncm6 (XCM,YCM,Sig,I,J,z,k)
global E;
M6=0;
i=floor(YCM-k*Sig(6,2));
j=ceil(XCM+k*Sig(6,1));
Ucm6=0;
if i==0 || j==0
  Ucm6=0; 
else
while i<=ceil(YCM+k*Sig(6,2)) && j>=floor(XCM-k*Sig(6,1))
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
                M6=M6+1;
              if E(i,j,z+o)==1
              Ucm6=Ucm6+1;        
              end
            end
        case 2
            for o=-1:1
                M6=M6+1;
              if E(i,j-1,z+o)==1
              Ucm6=Ucm6+1;        
              end
            end
        case 3
            
            for o=-1:1
                M6=M6+1;
              if E(i+1,j,z+o)==1
              Ucm6=Ucm6+1;        
              end
            end
       end
   end
   i=i+1;
   j=j-1;
end
Ucm6=Ucm6/M6;
end
end

function [Ucm7] = evaluncm7 (XCM,Sig,i,J,z,k)
global E;
M7=0;
j=ceil(XCM+k*Sig(7,1));
Ucm7=0;
if j==0
  Ucm7=0; 
else
while j>=floor(XCM-k*Sig(7,1))
   for q=-1:1
       for o =-1:1
            M7=M7+1;
            if E(i+q,j,z+o)==1
              Ucm7=Ucm7+1;
            end
           
       end
   end
   j=j-1;
end
Ucm7=Ucm7/M7;
end
end

function [Ucm8] = evaluncm8 (XCM,YCM,Sig,I,J,z,k)
global E;
M8=0;
i=ceil(YCM+k*Sig(8,2));
j=ceil(YCM+k*Sig(8,1));
Ucm8=0;
if i==0 || j==0;
  Ucm8=0; 
else
while i>=floor(YCM-k*Sig(8,2)) && j>=floor(XCM-k*Sig(8,1))
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
                M8=M8+1;
              if E(i,j,z+o)==1
              Ucm8=Ucm8+1;        
              end
            end
        case 2
            for o=-1:1
                M8=M8+1;
              if E(i,j-1,z+o)==1
              Ucm8=Ucm8+1;        
              end
            end
        case 3
            for o=-1:1
                M8=M8+1;
              if E(i-1,j,z+o)==1
              Ucm8=Ucm8+1;        
              end
            end
       end
   end
   i=i-1;
   j=j-1;
end
Ucm8=Ucm8/M8;
end
end


