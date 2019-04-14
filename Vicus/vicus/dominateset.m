function  [PNN_matrix, loc]  = dominateset(aff_matrix,NR_OF_KNN);

[res,loc] = maxk(aff_matrix, NR_OF_KNN, 2 );
inds = repmat((1:size(loc,1))',[1 size(loc,2)]);

PNN_matrix1 = zeros(size(aff_matrix));
PNN_matrix1(sub2ind(size(aff_matrix),inds(:),loc(:))) = res(:);
PNN_matrix = full(PNN_matrix1+PNN_matrix1')/2;

end

