function [K P pres dres] = lmi_solve(lmi)
%LMI_SOLVE

% solve
solvesdp(lmi.F);
[pres dres]=checkset(lmi.F)

%   pres : Primal constraint residuals
%   dres : Dual constraint residuals
%  
%   If no output argument is supplied, tabulated results are displayed
%  
%   Primal constraint residuals are calculated as:
%  
%    Semidefinite constraint F(x)>0 : min(eig(F))
%    Element-wise constraint F(x)>0 : min(min(F))
%    Equality constraint F==0       : -max(max(abs(F)))
%    Second order cone t>||x||      : t-||x||
%    Integrality constraint on x    : max(abs(x-round(x)))
%    Rank constraint rank(X) < r     : r-rank(X)
%    Sum-of-square constraint       : Minus value of largest (absolute value) coefficient 
%                                     in the polynomial p-v'*v
%  
%   Dual constraints are evaluated similarily.
% The constraint residuals are defined as smallest eigenvalue, smallest element, 
% negated largest absolute-value element and largest distance to an integer for 
% semidefinite inequality, element-wise inequality, equality and integrality constraints 
% respectively. Hence, a solution is feasible if all residuals related to inequalities are non-negative.
P = inv(double(lmi.X));
R = size(lmi.A, 1);
K = zeros(R, lmi.m, lmi.n);
for r = 1:R
	K(r,:,:) = double(lmi.M{r}) * P;
end
K = reshape(K, [lmi.sizes lmi.m lmi.n]);
