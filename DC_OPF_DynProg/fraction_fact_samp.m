function samp_id = fraction_fact_samp(cand_n)
% This function returns the fractional factorial level IV design, so main
% effects and two way interactions are preserved

%History            
%Version    Date        Who     Summary
%1          10/07/2017  JesseB  Initial Version

% ToDo: 


one_idx = nchoosek(1:cand_n,2);
[row, col] = size(one_idx);
A = zeros(row,cand_n);
for r_idx = 1:row
    for c_idx = 1:col
        A(r_idx, one_idx(r_idx,c_idx)) = 1;
    end
end

samp = fracfact(fracfactgen([eye(cand_n);A],(1+floor(2*sqrt(cand_n)))));
samp_id = bi2de((samp+1)./2);
end