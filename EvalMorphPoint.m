function [Xmin,Ymin,Xmax,Ymax,XCM,YCM,Sig]= EvalMorphPoint (i,j,z,N)
%Evalua los parametros de evolucion morfologica dado un punto de evaluacion (i,j,z) en E.
%Determina la celula más cercana y la más lejana en cada dirección de
%evaluación & determina la ubicación de cada CM dada una de las
%direcciones.

%EvalMorphPoint2 is a back up, only missing to includeSigma in the output of
%functions (and completing function 8)
global E;

Xmin=[];
Ymin=[];
Xmax=[];
Ymax=[];
XCM=[];
YCM=[];
Sig=[];

[Xmin(1), Ymin(1), Xmax(1), Ymax(1),  XCM(1), YCM(1),Sig(1,1),Sig(1,2)]=eval1(i,j,z,N);
[Xmin(2), Ymin(2), Xmax(2), Ymax(2),  XCM(2), YCM(2),Sig(2,1),Sig(2,2)]=eval2(i,j,z,N);
[Xmin(3), Ymin(3), Xmax(3), Ymax(3),  XCM(3), YCM(3),Sig(3,1),Sig(3,2)]=eval3(i,j,z,N);
[Xmin(4), Ymin(4), Xmax(4), Ymax(4),  XCM(4), YCM(4),Sig(4,1),Sig(4,2)]=eval4(i,j,z,N);
[Xmin(5), Ymin(5), Xmax(5), Ymax(5),  XCM(5), YCM(5),Sig(5,1),Sig(5,2)]=eval5(i,j,z,N);
[Xmin(6), Ymin(6), Xmax(6), Ymax(6),  XCM(6), YCM(6),Sig(6,1),Sig(6,2)]=eval6(i,j,z,N);
[Xmin(7), Ymin(7), Xmax(7), Ymax(7),  XCM(7), YCM(7),Sig(7,1),Sig(7,2)]=eval7(i,j,z,N);
[Xmin(8), Ymin(8), Xmax(8), Ymax(8),  XCM(8), YCM(8),Sig(8,1),Sig(8,2)]=eval8(i,j,z,N);


end



function [min1x, min1y, max1x, max1y, CMx1, CMy1,Sig11,Sig12]= eval1 (i,j,z,N)
global E;
%size(E)
iInitial=i;
jInitial=j;
min1=0;
max1=0;
min1y=0;
min1x=0;
max1y=0;
max1x=0;
meanx=0;
meany=0;
Rx1=0;
Ry1=0;
M1=0;
difx=0;
dify=0;
emptyspace=0;


while i>=2
   for h=-1:1
       for o =-1:1  
           if E(i,j+h,z+o)==1
              M1=M1+1;
              Rx1=Rx1+j+h;
              Ry1=Ry1+i;
              meanx=meanx+j+h;
              meany=meany+i;
           end
    
           if(min1==0)
           if(E(i,j+h,z+o)==1)
                min1x=j+h;
                min1y=i;
                min1z=z+o;
                min1=1;
           end
           end
   
           if(E(i,j+h,z+o)==1)
               if(max1y ~= i)
                max1x=j+h;
                max1y=i;
                max1z=z+o;
               end
           end
       end
   end
   i=i-1;
end

if M1==0
CMx1=0;
CMy1=0;
xm=0;
ym=0;
emptyspace=1;
else
CMx1=Rx1/M1;
CMy1=Ry1/M1; 
xm=meanx/M1;
ym=meany/M1;
end

i=iInitial;
while i>=2
   for h=-1:1
       for o =-1:1
           if E(i,j+h,z+o)==1
              difx=difx+(j+h-xm)^2;
              dify=dify+(i-ym)^2;
           end
       end
   end
   i=i-1;
end

if emptyspace ==0
Sig11=(difx/M1)^0.5;
Sig12=(dify/M1)^0.5;
else
Sig11=0;
Sig12=0;
end

end


%%%Falta actualizar calculo de la matriz Sig (desviacion estandar) para las
%%%demàs subfunciones e incluirlo como output de la funcion principal


function [min2x, min2y, max2x, max2y, CMx2, CMy2, Sig21, Sig22]= eval2 (i,j,z,N)
global E;
min2=0;
max2=0;
M2=0;
iInitial=i;
jInitial=j;
meanx=0;
meany=0;
Rx2=0;
Ry2=0;
difx=0;
dify=0;
min2y=0;
min2x=0;
max2y=0;
max2x=0;
emptyspace=0;

while i>=2 && j<=N-1
   for count=1:3
    switch (count)
        case 1
            for o=-1:1
            
           if E(i,j,z+o)==1
                  M2=M2+1;
                  Rx2=Rx2+j;
                  Ry2=Ry2+i;
                  meanx=meanx+j;
                  meany=meany+i;
           end
              
           if(min2==0)
           if(E(i,j,z+o)==1)
                min2x=j;
                min2y=i;
                min2z=z+o;
                min2=1;
           end
           end
           
           if(E(i,j,z+o)==1)
                max2x=j;
                max2y=i;
                max2z=z+o;
           end
           
            end
        case 2
            for o=-1:1
                if E(i,j+1,z+o)==1
                    M2=M2+1;
                    Rx2=Rx2+j+1;
                    Ry2=Ry2+i;
                    meanx=meanx+j+1;
                    meany=meany+i;
                end
                
                if(min2==0)
                    if(E(i,j+1,z+o)==1)
                        min2x=j+1;
                        min2y=i;
                        min2z=z+o;
                        min2=1;
                    end
                end
                
                if(E(i,j+1,z+o)==1)
                    max2x=j+1;
                    max2y=i;
                    max2z=z+o;
                end
           end
        case 3
            for o=-1:1
              
               if E(i-1,j,z+o)==1
                    M2=M2+1;
                    Rx2=Rx2+j;
                    Ry2=Ry2+i-1;
                    meanx=meanx+j;
                    meany=meany+i-1;
                end
                
                if(min2==0)
                    if(E(i-1,j,z+o)==1)
                        min2x=j;
                        min2y=i-1;
                        min2z=z+o;
                        min2=1;
                    end
                end
                
                if(E(i-1,j,z+o)==1)
                    max2x=j;
                    max2y=i-1;
                    max2z=z+o;
                end
              
            end
    end
    end
   i=i-1;
   j=j+1;
end

if M2==0
CMx2=0;
CMy2=0;
xm=0;
ym=0;
emptyspace=1;
else
CMx2=Rx2/M2;
CMy2=Ry2/M2; 
xm=meanx/M2;
ym=meany/M2;
end

i=iInitial;
j=jInitial;

while i>=2 && j<=N-1
   for count=1:3
    switch (count)
        case 1
            for o=-1:1
            
           if E(i,j,z+o)==1 
                  difx=difx+(j-xm)^2;
                  dify=dify+(i-ym)^2;
                %  i
                %  j
           end
            end
        case 2
            for o=-1:1
                if E(i,j+1,z+o)==1
                  difx=difx+(j+1-xm)^2;
                  dify=dify+(i-ym)^2;
                 % i
                 %j+1
                end
           end
        case 3
            for o=-1:1
              
               if E(i-1,j,z+o)==1
                  difx=difx+(j-xm)^2;
                  dify=dify+(i-1-ym)^2;
                 % i-1
                 % j
               end              
            end
    end
    end
   i=i-1;
   j=j+1;
end

if emptyspace ==0
Sig21=(difx/M2)^0.5;
Sig22=(dify/M2)^0.5;
else
Sig21=0;
Sig22=0;
end

end

function [min3x, min3y, max3x, max3y, CMx3, CMy3, Sig31, Sig32]= eval3 (i,j,z,N)
global E;
min3=0;
max3=0;
M3=0;
iInitial=i;
jInitial=j;
meanx=0;
meany=0;
Rx3=0;
Ry3=0;
difx=0;
dify=0;
min3y=0;
min3x=0;
max3y=0;
max3x=0;
emptyspace=0;

while j<=N-1
   for q=-1:1
       for o =-1:1
    
           if E(i+q,j,z+o)==1
              M3=M3+1;
              Rx3=Rx3+j;
              Ry3=Ry3+i+q;
              meanx=meanx+j;
              meany=meany+i+q;
           end
           
           if(min3==0)
           if(E(i+q,j,z+o)==1)
                min3x=j;
                min3y=i+q;
                min3z=z+o;
                min3=1;
           end
           end
           
           
           if(E(i+q,j,z+o)==1)
               if(max3x ~= j)
                max3x=j;
                max3y=i+q;
                max3z=z+o;
               end
           end
       end
   end
   j=j+1;
end 
if M3==0
CMx3=0;
CMy3=0;
xm=0;
ym=0;
emptyspace=1;
else
CMx3=Rx3/M3;
CMy3=Ry3/M3; 
xm=meanx/M3;
ym=meany/M3;
end

i=iInitial;
j=jInitial;

while j<=N-1
   for q=-1:1
       for o =-1:1
           if E(i+q,j,z+o)==1
              difx=difx+(j-xm)^2;
              dify=dify+(i+q-ym)^2;
           end
       end
   end
   j=j+1;
end

if emptyspace ==0
Sig31=(difx/M3)^0.5;
Sig32=(dify/M3)^0.5;
else
Sig31=0;
Sig32=0;
end

end

function [min4x, min4y, max4x, max4y, CMx4, CMy4, Sig41, Sig42]= eval4 (i,j,z,N)
global E;
min4=0;
max4=0;
M4=0;
iInitial=i;
jInitial=j;
meanx=0;
meany=0;
Rx4=0;
Ry4=0;
difx=0;
dify=0;
min4y=0;
min4x=0;
max4y=0;
max4x=0;
emptyspace=0;

while i<=N-1 && j<=N-1
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
            
                if E(i,j,z+o)==1
                    M4=M4+1;
                    Rx4=Rx4+j;
                    Ry4=Ry4+i;
                    meanx=meanx+j;
                    meany=meanx+i;
                end
                
                if(min4==0)
                    if(E(i,j,z+o)==1)
                        min4x=j;
                        min4y=i;
                        min4z=z+o;
                        min4=1;
                    end
                end
                
                if(E(i,j,z+o)==1)
                    max4x=j;
                    max4y=i;
                    max4z=z+o;
                end
                
            end
        case 2
            for o=-1:1
                if E(i,j+1,z+o)==1
                    M4=M4+1;
                    Rx4=Rx4+j+1;
                    Ry4=Ry4+i;
                    meanx=meanx+j+1;
                    meany=meanx+i;
                end
                
                if(min4==0)
                    if(E(i,j+1,z+o)==1)
                        min4x=j+1;
                        min4y=i;
                        min4z=z+o;
                        min4=1;
                    end
                end
                
                if(E(i,j+1,z+o)==1)
                    max4x=j+1;
                    max4y=i;
                    max4z=z+o;
                end
            end
        case 3
            for o=-1:1
                
                if E(i+1,j,z+o)==1
                    M4=M4+1;
                    Rx4=Rx4+j;
                    Ry4=Ry4+i+1;
                    meanx=meanx+j;
                    meany=meanx+i+1;
                end
                
                if(min4==0)
                    if(E(i+1,j,z+o)==1)
                        min4x=j;
                        min4y=i+1;
                        min4z=z+o;
                        min4=1;
                    end
                end
                
                if(E(i+1,j,z+o)==1)
                    max4x=j;
                    max4y=i+1;
                    max4z=z+o;
                end
                
            end
    end
    end
i=i+1;
j=j+1;
end
if M4==0
CMx4=0;
CMy4=0;
xm=0;
ym=0;
emptyspace=1;
else
CMx4=Rx4/M4;
CMy4=Ry4/M4; 
xm=meanx/M4;
ym=meany/M4;
end

i=iInitial;
j=jInitial;

while i<=N-1 && j<=N-1
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
            
                if E(i,j,z+o)==1
                  difx=difx+(j-xm)^2;
                  dify=dify+(i-ym)^2;
                end
                
            end
        case 2
            for o=-1:1
                if E(i,j+1,z+o)==1
                  difx=difx+(j+1-xm)^2;
                  dify=dify+(i-ym)^2;
                end
            end
        case 3
            for o=-1:1
                
                if E(i+1,j,z+o)==1
                  difx=difx+(j-xm)^2;
                  dify=dify+(i+1-ym)^2;
                end
            end
    end
    end
i=i+1;
j=j+1;
end

if emptyspace ==0
Sig41=(difx/M4)^0.5;
Sig42=(dify/M4)^0.5;
else
Sig41=0;
Sig42=0;
end

end

function [min5x, min5y, max5x, max5y, CMx5, CMy5, Sig51, Sig52]= eval5 (i,j,z,N)
global E;
min5=0;
max5=0;
M5=0;
iInitial=i;
jInitial=j;
meanx=0;
meany=0;
Rx5=0;
Ry5=0;
difx=0;
dify=0;
min5y=0;
min5x=0;
max5y=0;
max5x=0;
emptyspace=0;

while i<=N-1
   for h=-1:1
       for o =-1:1
     
           if E(i,j+h,z+o)==1
              M5=M5+1;
              Rx5=Rx5+j+h;
              Ry5=Ry5+i;
              meanx=meanx+j+h;
              meany=meany+i;
           end
           
           if(min5==0)
           if(E(i,j+h,z+o)==1)
                min5x=j+h;
                min5y=i;
                min5z=z+o;
                min5=1;
           end
           end
   
           if(E(i,j+h,z+o)==1)
               if(max5y ~= i)
                max5x=j+h;
                max5y=i;
                max5z=z+o;
               end
           end
       end
   end
   i=i+1;
end

if M5==0
CMx5=0;
CMy5=0;
xm=0;
ym=0;
emptyspace=1;
else
CMx5=Rx5/M5;
CMy5=Ry5/M5; 
xm=meanx/M5;
ym=meany/M5;
end


i=iInitial;
j=jInitial;

while i<=N-1
   for h=-1:1
       for o =-1:1    
           if E(i,j+h,z+o)==1
              difx=difx+(j+h-xm)^2;
              dify=dify+(i-ym)^2;
           end
       end
   end
   i=i+1;
end

if emptyspace ==0
Sig51=(difx/M5)^0.5;
Sig52=(dify/M5)^0.5;
else
Sig51=0;
Sig52=0;
end

end

function [min6x, min6y, max6x, max6y, CMx6, CMy6, Sig61, Sig62]= eval6 (i,j,z,N)
global E;
min6=0;
max6=0;
M6=0;
iInitial=i;
jInitial=j;
meanx=0;
meany=0;
Rx6=0;
Ry6=0;
difx=0;
dify=0;
min6y=0;
min6x=0;
max6y=0;
max6x=0;
emptyspace=0;

while i<=N-1 && j>=2
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
                
                if E(i,j,z+o)==1
                    M6=M6+1;
                    Rx6=Rx6+j;
                    Ry6=Ry6+i;
                    meanx=meanx+j;
                    meany=meany+i;
                    
                end
                
                if(min6==0)
                    if(E(i,j,z+o)==1)
                        min6x=j;
                        min6y=i;
                        min6z=z+o;
                        min6=1;
                    end
                end
                
                if(E(i,j,z+o)==1)
                    max6x=j;
                    max6y=i;
                    max6z=z+o;
                end
                
            end
        case 2
            for o=-1:1
                if E(i,j-1,z+o)==1
                    M6=M6+1;
                    Rx6=Rx6+j-1;
                    Ry6=Ry6+i;
                    meanx=meanx+j-1;
                    meany=meany+i;
                    
                end
                
                if(min6==0)
                    if(E(i,j-1,z+o)==1)
                        min6x=j-1;
                        min6y=i;
                        min6z=z+o;
                        min6=1;
                    end
                end
                
                if(E(i,j-1,z+o)==1)
                    max6x=j-1;
                    max6y=i;
                    max6z=z+o;
                end
            end
        case 3
            for o=-1:1
                
                if E(i+1,j,z+o)==1
                    M6=M6+1;
                    Rx6=Rx6+j;
                    Ry6=Ry6+i+1;
                    meanx=meanx+j;
                    meany=meany+i+1;
                    
                end
                
                if(min6==0)
                    if(E(i+1,j,z+o)==1)
                        min6x=j;
                        min6y=i+1;
                        min6z=z+o;
                        min6=1;
                    end
                end
                
                if(E(i+1,j,z+o)==1)
                    max6x=j;
                    max6y=i+1;
                    max6z=z+o;
                end
                
            end
    end
    end
   i=i+1;
   j=j-1;
end

if M6==0
CMx6=0;
CMy6=0;
xm=0;
ym=0;
emptyspace=1;
else
CMx6=Rx6/M6;
CMy6=Ry6/M6; 
xm=meanx/M6;
ym=meany/M6;
end

i=iInitial;
j=jInitial;

while i<=N-1 && j>=2
    for count=1:3
    switch (count)
        case 1
            for o=-1:1
                
                if E(i,j,z+o)==1
                  difx=difx+(j-xm)^2;
                  dify=dify+(i-ym)^2;
                end               
            end
        case 2
            for o=-1:1
                if E(i,j-1,z+o)==1
                  difx=difx+(j-1-xm)^2;
                  dify=dify+(i-ym)^2;
                end
            end
        case 3
            for o=-1:1               
                if E(i+1,j,z+o)==1
                  difx=difx+(j-xm)^2;
                  dify=dify+(i+1-ym)^2;
                end               
            end
    end
    end
   i=i+1;
   j=j-1;
end

if emptyspace ==0
Sig61=(difx/M6)^0.5;
Sig62=(dify/M6)^0.5;
else
Sig61=0;
Sig62=0;
end

end

function [min7x, min7y, max7x, max7y, CMx7, CMy7,  Sig71, Sig72]= eval7 (i,j,z,N)
global E;
min7=0;
max7=0;
M7=0;
iInitial=i;
jInitial=j;
meanx=0;
meany=0;
Rx7=0;
Ry7=0;
difx=0;
dify=0;
min7y=0;
min7x=0;
max7y=0;
max7x=0;
emptyspace=0;

while j>=2
   for q=-1:1
       for o =-1:1
    
            if E(i+q,j,z+o)==1
              M7=M7+1;
              Rx7=Rx7+j;
              Ry7=Ry7+i+q;
              meanx=meanx+j;
              meany=meany+i+q;
           end
           
           
           if(min7==0)
           if(E(i+q,j,z+o)==1)
                min7x=j;
                min7y=i+q;
                min7z=z+o;
                min7=1;
           end
           end
   
           if(E(i+q,j,z+o)==1)
               if(max7x ~= j)
                max7x=j;
                max7y=i+q;
                max7z=z+o;
               end
           end
       end
   end
   j=j-1;
end

if M7==0
CMx7=0;
CMy7=0;
xm=0;
ym=0;
emptyspace=1;
else
CMx7=Rx7/M7;
CMy7=Ry7/M7; 
xm=meanx/M7;
ym=meany/M7;
end

i=iInitial;
j=jInitial;
while j>=2
   for q=-1:1
       for o =-1:1
            if E(i+q,j,z+o)==1
              difx=difx+(j-xm)^2;
              dify=dify+(i+q-ym)^2;
            end
       end
   end
   j=j-1;
end

if emptyspace ==0
Sig71=(difx/M7)^0.5;
Sig72=(dify/M7)^0.5;
else
Sig71=0;
Sig72=0;
end

end

function [min8x, min8y, max8x, max8y, CMx8, CMy8,  Sig81, Sig82]= eval8 (i,j,z,N)
global E;
min8=0;
max8=0;
M8=0;
iInitial=i;
jInitial=j;
meanx=0;
meany=0;
Rx8=0;
Ry8=0;
difx=0;
dify=0;
min8y=0;
min8x=0;
max8y=0;
max8x=0;
emptyspace=0;

while i>=2 && j>=2
    for count=1:3
        switch (count)
            case 1
                for o=-1:1
                    
                    if E(i,j,z+o)==1
                        M8=M8+1;
                        Rx8=Rx8+j;
                        Ry8=Ry8+i;
                        meanx=meanx+j;
                        meany=meanx+i;
                    end
                    
                    if(min8==0)
                        if(E(i,j,z+o)==1)
                            min8x=j;
                            min8y=i;
                            min8z=z+o;
                            min8=1;
                        end
                    end
                    
                    if(E(i,j,z+o)==1)
                        max8x=j;
                        max8y=i;
                        max8z=z+o;
                    end
                    
                end
            case 2
                for o=-1:1
                    if E(i,j-1,z+o)==1
                        M8=M8+1;
                        Rx8=Rx8+j-1;
                        Ry8=Ry8+i;
                        meanx=meanx+j-1;
                        meany=meanx+i;
                        
                    end
                    
                    if(min8==0)
                        if(E(i,j-1,z+o)==1)
                            min8x=j-1;
                            min8y=i;
                            min8z=z+o;
                            min8=1;
                        end
                    end
                    
                    if(E(i,j-1,z+o)==1)
                        max8x=j-1;
                        max8y=i;
                        max8z=z+o;
                    end
                end
            case 3
                for o=-1:1
                    
                    if E(i-1,j,z+o)==1
                        M8=M8+1;
                        Rx8=Rx8+j;
                        Ry8=Ry8+i-1;
                        meanx=meanx+j;
                        meany=meanx+i-1;
                        
                    end
                    
                    if(min8==0)
                        if(E(i-1,j,z+o)==1)
                            min8x=j;
                            min8y=i-1;
                            min8z=z+o;
                            min8=1;
                        end
                    end
                    
                    if(E(i-1,j,z+o)==1)
                        max8x=j;
                        max8y=i-1;
                        max8z=z+o;
                    end
                    
                end
        end
    end
   i=i-1;
   j=j-1;
end

if M8==0
CMx8=0;
CMy8=0;
xm=0;
ym=0;
emptyspace=1;
else
CMx8=Rx8/M8;
CMy8=Ry8/M8; 
xm=meanx/M8;
ym=meany/M8;
end

i=iInitial;
j=jInitial;

while i>=2 && j>=2
    for count=1:3
        switch (count)
            case 1
                for o=-1:1   
                    if E(i,j,z+o)==1
                    difx=difx+(j-xm)^2;
                    dify=dify+(i-ym)^2;
                    end
                end
            case 2
                for o=-1:1
                    if E(i,j-1,z+o)==1
                      difx=difx+(j-1-xm)^2;
                      dify=dify+(i-ym)^2;
                    end

                end
            case 3
                for o=-1:1
                    
                    if E(i-1,j,z+o)==1
                     difx=difx+(j-xm)^2;
                     dify=dify+(i-1-ym)^2;
                        
                    end     
                end
        end
    end
   i=i-1;
   j=j-1;
end
    
if emptyspace ==0
Sig81=(difx/M8)^0.5;
Sig82=(dify/M8)^0.5;
else
Sig81=0;
Sig82=0;
end

end

