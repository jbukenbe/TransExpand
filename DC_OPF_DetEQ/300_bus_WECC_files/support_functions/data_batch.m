function data_batch(filename, n)
% This is a batch function to quickly group the output data from a large
% array of runs on the HPC where the output is of the form: filename_n for
% n output files.

%History            
%Version    Date        Who     Summary
%1          04/04/2017  JesseB  Initial Version


data_list = [1:16];
output = [];
n = length(data_list);
for n_idx = 1:n
    data_idx = data_list(n_idx); 
    outfile_name = sprintf('%s_%d',filename,data_idx);
    load(outfile_name);
    output = [output;scen_op_cost(:,1:10)'];
end
out_name = sprintf('%ss_1_to_%d',filename,n);
m = matfile(out_name);
m.output = output;

end