function [LocalSpectrum] = LS_CallocalLaplacian_fast1(W,K)
tic
A = LS_localdiffuse_KNN(W,K);
toc
M = (A + A' - A'*A);
%M = (A+A')/2;
LocalSpectrum = (eye(length(M)) - (M));
end

function Coef2 = LS_localdiffuse_KNN(W,K)
[n] = size(W,1);
W = W - diag(diag(W));
W = W/max(max(W));%%Normalization
if sum(sum(W>0))>0
W = W + rand(size(W))*min(min(W(W>0)))/length(W)/length(W);
end

W = W+W';
%W = TransitionFields(dominateset(W,min(round(K),n-1)));
[NB,NBind] = dominateset(W,K);
NB = NB+NB';%%%This is to find the top K neighbors for each node using MaxminSelection Package
beta = 0.5;
%[V, D] = ncuts(NB, n-1);

%beta = mean(D)

Coef2 = zeros(size(W));

num=zeros(n,1);
%W = (W-min(min(W)))/(max(W(:))-min(W(:)));
for i =1: n
    indexi = (find(NB(i,:)>0));
    temp = NBind(indexi,:);
    indexi = union(i,temp(:));
    %indexi = indexi(indexi>=i);
    if length(indexi)==1
        temp = 1;
    else
        Wi = NB(indexi,indexi);
        WiD = (Wi+Wi')/2;
        
%         WiD = WiD - diag(diag(WiD));
%         WiD = NE_dn(WiD, 'ave');WiD = (WiD+WiD')/2;
%         
         WiD = WiD + (eye(length(WiD)) + diag(sum(WiD)));
%         [Vi,Di] = eig(WiD);
%         Di = diag(Di);[a,b] = sort(Di,'descend');
%         beta = b(2)/b(1)
        if sum(isnan(WiD(:)))+sum(isinf(WiD(:)))>0
            temp = zeros(size(WiD));
        else
            S = NE_dn(WiD, 'ave');
            S = S^2;
            IW = (1-beta)*inv((1+eps)*eye(length(WiD))-beta*S);
            temp = IW./repmat(1-diag(IW)+eps,1,length(IW));
            temp = temp - diag(diag(temp));
            if sum(sum(isnan(temp)))>0
                i
            end
        end
    end
    if sum(sum(isnan(temp)))==0
        Coef2(indexi,indexi) = (temp + Coef2(indexi,indexi));
        %Coef2(i,indexi) = (temp(indexi==i,:) + Coef2(i,indexi));
        num(indexi) = num(indexi) + 1;
    end
%%%%%%    
%     indexj = [];
%     for j = 1:length(indexi)
%         indexj = union(indexj,find(NB(indexi(j),:)>0));
%     end
%     indexj = union(indexj,indexi);
%     if length(indexj)==1
%         temp = 1;
%     else
%         Wi = NB(indexj,indexj)+W(indexj,indexj);
%         WiD = (Wi+Wi')/2;
%         WiD = WiD - diag(diag(WiD));
%         %WiD = WiD + eye(length(WiD)) + diag(max(WiD));
%         if sum(isnan(WiD(:)))+sum(isinf(WiD(:)))>0
%             temp = zeros(size(WiD));
%         else
%             IW = zeros(size(WiD));
%             S = dn(WiD, 'ave');
%             IW = IW + (1-beta)*inv((1+eps)*eye(length(WiD))-beta*S);
%             temp = IW./repmat(1-diag(IW)+eps,1,length(IW));
%             temp = temp - diag(diag(temp));
%         end
%     end
%     [C,IA,IB] = intersect(indexi,indexj);
%     Coef2(indexi,indexi) = (temp(IB,IB) + Coef2(indexi,indexi));
%     num(indexi) = num(indexi) + 1;
%     Coef2(i,indexj) = (temp(indexj==i,:) + Coef2(i,indexj));
%     num(indexj) = num(indexj) + 1;
    
%%%%    
    
end

Coef2 = (Coef2./repmat(num,1,n));

% %%%
%
% indexi = 1:n;
% if length(indexi)==1
%     temp = 1;
% else
%     Wi = NB(indexi,indexi);
%     WiD = (Wi+Wi')/2;
%     WiD = WiD - diag(diag(WiD));
%     WiD = WiD + (eye(length(WiD)) + diag(max(WiD)));
%     if sum(isnan(WiD(:)))+sum(isinf(WiD(:)))>0
%         temp = zeros(size(WiD));
%     else
%         IW = zeros(size(WiD));
%         S = dn(WiD, 'gph');
%         IW = IW + (1-beta)*inv((1+eps)*eye(length(WiD))-beta*S);
%         temp = IW./repmat(1-diag(IW)+eps,1,length(IW));
%         temp = temp - diag(diag(temp));
%     end
% end
% Coef2(indexi,indexi) = (temp.*(NB>0) + Coef2(indexi,indexi))/2;
%
%


end
