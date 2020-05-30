function [Uc]=convTPmodel(dim,U,type)
for i=1:dim
    Uc{i}=genhull(U{i},type);
    Ucp{i}=pinv(Uc{i});
end
%B=tprods(S,Ucp);