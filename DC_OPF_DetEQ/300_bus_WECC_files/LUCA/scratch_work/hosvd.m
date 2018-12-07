function problem = hosvd()
% This function is scratch work to work out the high dimentional singular
% value decomposition for the TNEP problem with generation uncertainty

%History            
%Version    Date        Who     Summary
%1          11/12/2018  JesseB  scratch work

%%
n = 200;
n_hour = 10;
n_network = 5;
n_gen = 4;


A = reshape((rand(1,n)+(1:n)),n_hour, n./n_hour);
A = A - mean(A);
A = reshape(A,n_hour, n_network, n_gen);


[Uh,Sh,Vh] = svd(reshape(A,n_hour,n./n_hour));
[Ut,St,Vt] = svd(reshape(permute(A,[1,3,2]), n./n_network, n_network)');
[Ug,Sg,Vg] = svd(reshape(A, n./n_gen, n_gen)');

hr_dim = 3;
t_dim = 2;
g_dim = 2;

Uh = Uh(:, 1:hr_dim);
Ut = Ut(:, 1:t_dim);
Ug = Ug(:, 1:g_dim);

g = reshape(Uh'*reshape(A,n_hour,n./n_hour), hr_dim, n_network, n_gen);
ng = nnz(g);
g = reshape((Ug'*reshape(g, ng./n_gen, n_gen)')', hr_dim, n_network, g_dim);
nt = nnz(g);
g = permute(reshape((Ut'*reshape(permute(g,[1,3,2]), nt./n_network, n_network)')', hr_dim, g_dim, t_dim),[1 3 2]);

A_fit = reshape(Uh*reshape(g,n_hour,n./n_hour), n_hour, n_network, n_gen);
A_fit = reshape((Ug*reshape(A_fit, n./n_gen, n_gen)')', n_hour, n_network, n_gen);
A_fit = permute(reshape((Ut*reshape(permute(A_fit,[1,3,2]), n./n_network, n_network)')', n_hour, n_gen, n_network),[1 3 2]);

end




