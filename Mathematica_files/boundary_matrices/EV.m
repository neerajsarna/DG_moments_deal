clear all;

nEqn = [10, 16, 28, 40, 60, 80, 110, 140, 182, 224, 280, 336, 408, 480, 570, ...
660, 770, 880];
 
for i = nEqn

disp("Equations: ");
disp(i);
filename = strcat("../../system_matrices/monotone_systems/A1_2D_",num2str(i),"Trun.txt");
fileID = fopen(filename,"r");
indices = dlmread(filename);
fclose(fileID);

ia = indices(2:end,1)+1;
ja = indices(2:end,2)+1;
v = indices(2:end,3);

Ax = full(sparse(ia,ja,v,880,880));

[V,D] = eig(Ax);

Ax_mod = sparse(V * abs(D) * V');

[row, col, v] = find(Ax_mod);

row = row-1;
col = col-1;

nnz = length(row);

filename = strcat("../../system_matrices/monotone_systems/Amod1_2D_",num2str(i),"Trun.txt");
dlmwrite(filename,nnz);
dlmwrite(filename,[row col v],'-append',...
'delimiter',' ');

end


