function approx_out = mat_fill_svd(svd_values, svd_directions,partial_vec)
% This function makes a Singular Value Decomposition of an incomplete
% matrix and returns the estimated full matrix value

%History            
%Version    Date        Who     Summary
%1          10/10/2017  JesseB  Initial Version
%2          12/03/2017  JesseB  Integrated into line search


% format singular values into square matrix
s_n = size(svd_values,2);
svd_values = svd_values(1:s_n,1:s_n);

loadings = svd_values*svd_directions';
loadings = loadings(1:3,:);

%gradient decent to tune partial_fill U to fit known results as much as possible
loading_n = size(loadings,1);
u_approx = zeros(1,loading_n);
known_idx = find(partial_vec);        
for er_idx = 1:250   
    for j_idx = 1:length(known_idx)
        err = partial_vec -  u_approx*loadings;
        k_idx = known_idx(j_idx);
        delta =  .00000000013*err(k_idx).*loadings(:,k_idx);
        u_approx = u_approx + delta'; 
    end
end
approx_out = u_approx*loadings;
part_fill_logic = (partial_vec~=0);
approx_out(part_fill_logic)= partial_vec(part_fill_logic);
end

