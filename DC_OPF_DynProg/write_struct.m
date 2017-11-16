function write_struct(problem)
% This function takes the problem structure for the transmission expansion
% problem and saves the relavent output as a matfile. This way results from
% HPC batch runs can be saved and experimented on later

%History            
%Version    Date        Who     Summary
%1          11/02/2017  JesseB  Initial Version

filename = 'output.mat';
output = matfile(filename);
output.problem = problem;

end