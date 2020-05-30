function [K P,kim,diagnostic] = lmi_solve(lmi)
%LMI_SOLVE

% solve
options = sdpsettings;
options.savesolveroutput=1;
%diagnostic=solvesdp(lmi.F);
diagnostic=solvesdp(lmi.F,[],options);
checkset(lmi.F)

P = inv(double(lmi.X));
R = size(lmi.A, 1);
K = zeros(R, lmi.m, lmi.n);
for r = 1:R
	K(r,:,:) = double(lmi.M{r}) * P;
    
end
kim=K;
K = reshape(K, [lmi.sizes lmi.m lmi.n]);
