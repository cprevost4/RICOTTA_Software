function rmse = rmse(ref,est)

D = size(ref,1); R = size(ref,2);

rmse = 0;
temp = 0;
for r=1:R
    for d=1:D
        temp = temp + (ref(d,r)-est(d,r))^2;
    end
    temp = (1/D)*temp;
    temp = sqrt(temp);
    rmse = rmse + temp;
end

rmse = (1/R)*rmse;
end

