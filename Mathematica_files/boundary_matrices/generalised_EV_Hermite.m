clear all;

nEqn = [10, 16, 28, 40, 60, 80, 110, 140, 182, 224, 280, 336, 408, 480, 570, ...
660, 770, 880];

for i = nEqn
filename = strcat("../../system_matrices/monotone_systems/A1_2D_",num2str(i),"Trun.txt");
fileID = fopen(filename,"r");
indices = dlmread(filename);

ia = indices(2:end,1)+1;
ja = indices(2:end,2)+1;
v = indices(2:end,3);

Ax = full(sparse(ia,ja,v,i,i));

filename = strcat("../../system_matrices/monotone_systems/BGK_P_2D_",num2str(i),"Trun.txt");
fileID = fopen(filename,"r");
indices = dlmread(filename);

ia = indices(2:end,1)+1;
ja = indices(2:end,2)+1;
v = indices(2:end,3);

P0 = full(sparse(ia,ja,v,i,i));

[V,D] = eig(P0,Ax);

[~,permutation]=sort(diag(D));

D=D(permutation,permutation);V=V(:,permutation);

% eigenvalues in pair
pair = [];
for num_pairs = 1 : (length(diag(D)))/2
    pair = [num_pairs length(diag(D))-num_pairs+1 pair];
end


D=D(pair,pair);V=V(:,pair);

disp(norm(P0*V-Ax*V*D));

filename = strcat("heat_conduction_hermite/V",num2str(i),".txt");
dlmwrite(filename,V,'precision',16,'delimiter','\t');

filename = strcat("heat_conduction_hermite/D",num2str(i),".txt");
dlmwrite(filename,diag(D),'precision',16,'delimiter','\t');
end


function[sortV,sortD] = sort_ev(V,D)

% position of the nz entries in D
pos_nz = find(diag(D));

end
