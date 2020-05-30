function u = towercraneux_CalU(x,K,U,domain,n,q)
%% PDC ¸ú×Ù¿ØÖÆÆ÷
% 2018-01-18
% u = tpcontroller(p, x, K, U, domain);
%TPCONTROLLER Outputs the control signal for a TP controller at a given state
%	u = TPCONTROLLER(p, x, K, U, domain)
%	
%	p      - parameter point
%	x      - state vector
%	K      - core tensor of the feedback TP model
%	U      - weighting function data
%	domain - parameter domain (intervals)
%
%	u      - resulting control signal

p = [x(3) x(4)];
W = queryw1(U,domain,p);
z = tprod(K(:,:,:,1:n),W);
z1 = tprod(K(:,:,:,n+1),W);
z = shiftdim(z);
[n1 n2] = size(z);
% TODO
if n1 > n2
	u = -z(1:n1)' * x(:)-z1(end)*q;% tilde K=[K_i,-L_i],need to revert to original
else
	u = -z(1:n1) * x(:)-z1(end)*q;
end