function A_fill = mat_fill_svd(A)
% This function makes a Singular Value Decomposition of an incomplete
% matrix and returns the estimated full matrix value

%History            
%Version    Date        Who     Summary
%1          10/10/2017  JesseB  Initial Version


zero_map = (A == 0);
val_map = ~zero_map;

num_A_val = sum(val_map);
sum_A_col = sum(A);

mean_A = sum_A_col./num_A_val;
std_A = zeros(size(mean_A));

for s_idx = 1:size(A,2)
    A_col = A(val_map(:,s_idx));
    std_A(s_idx) =  std(A_col);
end

A_mean_fill = A - mean_A.*val_map;
A_stand = A_mean_fill./std_A;


[U,S,V] = svd(A_stand);
S(:,3:end) = 0;
A_stand_fill = U*S*V';

A_fill = (A_stand_fill.*std_A)+ mean_A(ones(size(A,1),1),:);


end