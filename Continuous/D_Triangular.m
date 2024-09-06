function [ M ] = D_Triangular (k,t,b,row,column)
    M=zeros(row,column);
    for i=1:row
        for j=1:column

            M(i,j)=D_Tri_real(k,t,b);
            
        end
    end
end

function [ val ] = D_Tri_real(k,t,b)
    m=randi([1 10]);
    a=(t-k)/10;
    b=(b-t)/10;
    
    val=unifrnd((t-m*a),(t+m*b),1);
    
end
