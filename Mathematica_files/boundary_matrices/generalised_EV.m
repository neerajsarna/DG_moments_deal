clear all;


for i = 4 : 15
 Ntensors = get_ntensors(i);
filename = strcat("heat_conduction/AxOld",num2str(Ntensors),".txt");
fileID = fopen(filename,"r");
indices = dlmread(filename);

ia = indices(2:end,1);
ja = indices(2:end,2);
v = indices(2:end,3);

m = indices(1,1);
n = indices(1,2);

Ax = full(sparse(ia,ja,v,m,n));
P0 = -eye(m,n);
P0(1,1) = 0;

[V,D] = eig(P0,Ax);

[~,permutation]=sort(diag(D));

num_inf = isinf(diag(D));
num_inf = sum(num_inf(:)==1);

D=D(permutation,permutation);V=V(:,permutation);

% eigenvalues in pair
pair = [];
for num_pairs = 1 : (length(diag(D))-num_inf)/2
    pair = [num_pairs length(diag(D))-num_pairs-num_inf+1 pair];
end

% add the infinite eigenvalues
pair = [pair length(diag(D))-num_inf+1:length(diag(D))];


D=D(pair,pair);V=V(:,pair);

disp(norm(P0*V-Ax*V*D));

filename = strcat("heat_conduction/VOld",num2str(Ntensors),".txt");
dlmwrite(filename,V,'precision',16,'delimiter','\t');

filename = strcat("heat_conduction/DOld",num2str(Ntensors),".txt");
dlmwrite(filename,diag(D),'precision',16,'delimiter','\t');
end

function f = get_ntensors(id_theory)
f = sum(1:1:ceil(id_theory/2)) + sum(1:1:floor(id_theory/2));
end

function[sortV,sortD] = sort_ev(V,D)

% position of the nz entries in D
pos_nz = find(diag(D));



end
