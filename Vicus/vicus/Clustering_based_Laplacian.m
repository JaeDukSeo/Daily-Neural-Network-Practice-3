function [group,best_group_index,best_eigengap] = Clustering_based_Laplacian(L,NUMC);

[n] = length(L);

%  M = (eye(n) - A)'*(eye(n) - A);




%[group,best_group_index,Quality] = LS_cluster_rotate_Laplacian_Bo(L,NUMC,0,1);

M = full(L);
[V,D] = eig(M);
[a,b]=sort(real(diag(D)),'ascend');



for JJ = 1 : length(NUMC)
    C = NUMC(JJ);
    UU = V(:,b(1:C));
    UU = UU./repmat(sqrt(sum(UU.^2+eps)),size(UU,1),1);
    [EigenvectorsDiscrete,EigenVectors]=discretisation(UU);
    [temp,group1] = max(EigenvectorsDiscrete,[],2);
    %group1 = litekmeans(UU,C);
    group{JJ} = group1;
end


eigengap = abs(diff(a(1:max(NUMC+1))));
eigengap = eigengap(NUMC);
[best_eigengap, best_group_index] = max(eigengap);

if length(NUMC)==1
    group = group{1};
end