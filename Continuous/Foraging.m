function y = Foraging (Position,nghk,upperBound,lowerBound,PSize)
   
    r = nghk * PSize;
    nVar = numel (Position);
    k = randi([1 nVar]);
    y = Position;

    if numel(PSize)>1
        y(k) = y(k)+ unifrnd(-r(k),r(k),1);
        if y(k)>upperBound(k)
            y(k)=upperBound(k);
        end
        if y(k)<lowerBound(k)
            y(k)=lowerBound(k);
        end
    else
        y(k) = y(k)+ unifrnd(-r,r,1);
        y(y>upperBound) = upperBound;
        y(y<lowerBound) = lowerBound;       
    end
    
end
