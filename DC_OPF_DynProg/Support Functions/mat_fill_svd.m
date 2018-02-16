function [approx_out, u_approx] = mat_fill_svd(svd_values, svd_directions, partial_vec)
% This function makes a Singular Value Decomposition of an incomplete
% matrix and returns the estimated full matrix value

%History            
%Version    Date        Who     Summary
%1          10/10/2017  JesseB  Initial Version
%2          12/03/2017  JesseB  Integrated into line search


% format singular values into square matrix
s_n = min(size(svd_values));
svd_values = svd_values(1:s_n,1:s_n);

loadings = svd_values*svd_directions(:,1:s_n)';
loadings = loadings(:,:);

%gradient decent to tune partial_fill U to fit known results as much as possible
loading_n = size(loadings,1);
u_approx = zeros(1,loading_n);
known_idx = find(partial_vec);        
for er_idx = 1:400   
    for j_idx = 1:length(known_idx)
        err = partial_vec -  u_approx*loadings;
        k_idx = known_idx(j_idx);
        delta =  .00000000006*err(k_idx).*loadings(:,k_idx);
        u_approx = u_approx + delta'; 
    end
    approx_out_log(er_idx,:) = u_approx*loadings;
end
%plot(approx_out_log)
approx_out = u_approx*loadings;
part_fill_logic = (partial_vec~=0);
approx_out(part_fill_logic)= partial_vec(part_fill_logic);

end

