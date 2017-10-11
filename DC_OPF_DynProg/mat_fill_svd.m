function A_fill = mat_fill_svd(A,T)
% This function makes a Singular Value Decomposition of an incomplete
% matrix and returns the estimated full matrix value

%History            
%Version    Date        Who     Summary
%1          10/10/2017  JesseB  Initial Version


[U,S,V] = svd(A);

S = S(1:9,1:9);
Sless = S;
Sless(4:end,:)=0;
SV = Sless*V';
x = zeros(9,1);
T_est = zeros(size(T));
for t_idx = 1:length(T)
    T_est(t_idx,:) = scratch(SV',x,T(t_idx,:)');
end
A_fill = [A;T_est];

%num_A_val = sum(val_map);
%sum_A_col = sum(A);

%mean_A = sum_A_col./num_A_val;
%std_A = zeros(size(mean_A));

%for s_idx = 1:size(A,2)
%    A_col = A(val_map(:,s_idx));
%    std_A(s_idx) =  std(A_col);
%end

%A_mean_fill = A - mean_A.*val_map;
%A_stand = A_mean_fill./std_A;


%[U,S,V] = svd(A_stand);
%A_stand_fill = U*S*V';

%A_fill = (A_stand_fill.*std_A)+ mean_A(ones(size(A,1),1),:);


end