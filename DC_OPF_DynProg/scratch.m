function xest = scratch(A,x,b)
err = 2;
for n = 1:2000;
    for idx = 1:3
        err = b -  A*x;
        delta =  .00000000013*err(idx).*A(idx,:)';
        x = x +delta; 
    end
end

xest = A*x;
end