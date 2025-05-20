function UWM=UnWrapMatBest(MAT,m,n,tolMIN,tolMAX,disDIR)
% function UWM=UnWrapMatBest(MAT,m,n,tolMIN,tolMAX,disDIR)
% Unwraps phase images (one-directional) from a given starting point.
%
% INPUT:
% MAT = 2D phase image to unwrap
% m = starting coordinate x
% n = starting coordinate y
% tolMIN = minimum tolerance in degree
% tolMAX = maximum tolerance in degree
% disDIR = unwrapping direction (1=x, 2=y)
%
% OUTPUT:
% UWM = unwrapped 2D matrix

m=round(mod(abs(m)-1,size(MAT,1))+1);
n=round(mod(abs(n)-1,size(MAT,2))+1);

UWM=mod(MAT+180,360)-180;

if disDIR==1
    for j=n+1:size(UWM,2)
        if abs(UWM(m,j)-UWM(m,j-1))>=abs(tolMAX)
            UWM(m,j)=UWM(m,j)-360*sign(UWM(m,j)-UWM(m,j-1))*round(abs(UWM(m,j)-UWM(m,j-1))/360);
        elseif abs(UWM(m,j)-UWM(m,j-1))>=abs(tolMIN)
            UWM(m,j)=NaN;
        end
    end
    if n>1
        for j=n-1:-1:1
            if abs(UWM(m,j)-UWM(m,j+1))>=abs(tolMAX)
                UWM(m,j)=UWM(m,j)-360*sign(UWM(m,j)-UWM(m,j+1))*round(abs(UWM(m,j)-UWM(m,j+1))/360);
            elseif abs(UWM(m,j)-UWM(m,j+1))>=abs(tolMIN)
                UWM(m,j)=NaN;
            end
        end
    end
    
    for j=1:size(UWM,2)
        for i=m+1:size(UWM,1)
            if abs(UWM(i,j)-UWM(i-1,j))>=abs(tolMAX)
                UWM(i,j)=UWM(i,j)-360*sign(UWM(i,j)-UWM(i-1,j))*round(abs(UWM(i,j)-UWM(i-1,j))/360);
            elseif abs(UWM(i,j)-UWM(i-1,j))>=abs(tolMIN)
                UWM(i,j)=NaN;
            end
        end
        for i=m-1:-1:1
            if abs(UWM(i,j)-UWM(i+1,j))>=abs(tolMAX)
                UWM(i,j)=UWM(i,j)-360*sign(UWM(i,j)-UWM(i+1,j))*round(abs(UWM(i,j)-UWM(i+1,j))/360);
            elseif abs(UWM(i,j)-UWM(i+1,j))>=abs(tolMIN)
                UWM(i,j)=NaN;
            end
        end
    end
          
elseif disDIR==2
    
    for i=m+1:size(UWM,1)
        if abs(UWM(i,n)-UWM(i-1,n))>=abs(tolMAX)
            UWM(i,n)=UWM(i,n)-360*sign(UWM(i,n)-UWM(i-1,n))*round(abs(UWM(i,n)-UWM(i-1,n))/360);
        elseif abs(UWM(i,n)-UWM(i-1,n))>=abs(tolMIN)
            UWM(i,n)=NaN;
        end
    end
    if m>1
        for i=m-1:-1:1
            if abs(UWM(i,n)-UWM(i+1,n))>=abs(tolMAX)
                UWM(i,n)=UWM(i,n)-360*sign(UWM(i,n)-UWM(i+1,n))*round(abs(UWM(i,n)-UWM(i+1,n))/360);
            elseif abs(UWM(i,n)-UWM(i+1,n))>=abs(tolMIN)
                UWM(i,n)=NaN;
            end
        end
    end
         
    for i=1:size(UWM,1)
        for j=n+1:size(UWM,2)
            if abs(UWM(i,j)-UWM(i,j-1))>=abs(tolMAX)
                UWM(i,j)=UWM(i,j)-360*sign(UWM(i,j)-UWM(i,j-1))*round(abs(UWM(i,j)-UWM(i,j-1))/360);
            elseif abs(UWM(i,j)-UWM(i,j-1))>=abs(tolMIN)
                UWM(i,j)=NaN;
            end
        end
        for j=n-1:-1:1
            if abs(UWM(i,j)-UWM(i,j+1))>=abs(tolMAX)
                UWM(i,j)=UWM(i,j)-360*sign(UWM(i,j)-UWM(i,j+1))*round(abs(UWM(i,j)-UWM(i,j+1))/360);
            elseif abs(UWM(i,j)-UWM(i,j+1))>=abs(tolMIN)
                UWM(i,j)=NaN;
            end
        end
    end
       
else
    
    error('Direction selected not allowed (DIR=1 or DIR=2)')
    
end

end

