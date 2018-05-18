function samp = fraction_fact_samp(problem)
% This function returns the fractional factorial level IV design, so main
% effects and two way interactions are preserved

%History            
%Version    Date        Who     Summary
%1          10/07/2017  JesseB  Initial Version
%2          12/02/2017  JesseB  Reworked with walsh_index fractional factorial

cand_n = problem.params.cand.n(problem.z_idx);

if problem.z_idx > 1
    samp_size = problem.params.refine_samp_n;
else
    samp_size = problem.params.initial_samp_n;
end

hadamard_power = ceil(log(max(samp_size,cand_n)+1)/log(2));
hadamard_size = 2^hadamard_power;

walsh_index = [1 2 4 8 15 16 32 51 64 85 106 128 150 171 219 237 247 256 ...
    279 297 455 512 537 557 597 643 803 863 898 1024 1051 1070 1112 1169 ...
    1333 1345 1620 1866 2048 2076 2085 2185 2372 2456 2618 2800 2873 ...
    3127 3284 3483 3557 3763 4096 4125 4135 4174 4435 4459 4469 4497 ...
    4752 5255 5732 5804 5915 6100 6369 6907 7069 8192 8263 8351 8422 ...
    8459 8571 8750 8858 9124 9314 9500 10026 10455 10556 11778 11885 ...
    11984 13548 14007 14514 14965 15125 15554 16384 16457 16517 16609 ...
    16771 16853 17022 17453 17891 18073 18562 18980 19030 19934 20075 ...
    20745 21544 22633 23200 24167 25700 26360 26591 26776 28443 28905 ...
    29577 32705]+1;

walsh_avail = walsh_index(walsh_index < samp_size);
hadamard_samp = hadamard(hadamard_size, 'int8');
if (length(walsh_avail) > cand_n) 
    samp = hadamard_samp(1:(max(walsh_avail)-1),walsh_avail(1:cand_n));
else
    samp = hadamard_samp(1:samp_size,2:cand_n+1);
end
clear hadamard_samp
samp = ((samp +1)./2);
samp = uint8(samp);
samp = samp(:,randperm(cand_n));
samp = [eye(cand_n,'uint8');samp];
%{
if params.interaction 
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
else
    samp = fracfact(fracfactgen(eye(cand_n)));
    samp_id = bi2de((samp+1)./2); 
end
%}
end