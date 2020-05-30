function ParaStruct=TPPDC180415Noise(dt,x_0,y_r)
%% TPPDC
%% LPV model definition for crane
if nargin<3
    dt=0.1;
    x_0=[0.3,0,pi/4,0]';
    y_r=0.4;    
end
    B_eqx = 14;
    B_palpha = 0.0015;
    eta_gx = 0.55;
    eta_mx = 0.65;
    g = 9.81;
    J_cm= 9.44E-5;
    K_gx = 76.64;
    K_tx= 0.032;
    K_m=0.03;
    R_mx = 25;
    M = 2.78;
    m = 0.32;
    r_mp = 0.0375;
    L=0.6;
% system matrix: lpv = [A(p) B(p)]
D_1=@(p)((M+m)*J_cm+m*L^2*(M+m*sin(p(1))^2));
% p(t): x(3),x(4)
lpvux = {...
    @(p)0         @(p)1    @(p)0    @(p)0               @(p)0; 
    @(p)0 ...
    @(p)-(J_cm+m*L^2*B_eqx+eta_gx*K_gx^2 *eta_mx*K_tx*K_m*R_mx*r_mp^2) /D_1(p)...
    @(p)m^2*L^2*g*sinc(p(1)/pi)*cos(p(1))/D_1(p)...
    @(p)((m^2*L^3+L*m*J_cm)*sin(p(1))* p(2)+m*L*cos(p(1))*B_palpha)/D_1(p)...
    @(p)(-(m*L^2+J_cm)^2*eta_gx*K_gx*eta_mx*K_tx)/(D_1(p) *R_mx*r_mp^2);
   @(p)0         @(p)0    @(p)0    @(p)1               @(p)0;
   @(p)0 ...
   @(p)(m*L*cos(p(1))*B_eqx+eta_gx*K_gx^2 *eta_mx*K_tx*K_m*R_mx*r_mp^2)/D_1(p)...
   @(p)-(M+m)*m *g*L*sinc(p(1))/D_1(p)...
   @(p)-(M+m)*B_palpha-m^2*L^2*sin(p(1))*cos(p(1))* p(2)/D_1(p)...
   @(p)(-(m*L*cos(p(1)))^2*eta_gx*K_gx* eta_mx*K_tx)/(D_1(p)* R_mx*r_mp);
};

% number of states (size of the A matrix)
n1 = 4;
% parameter dep1endencies:
% dep1(i,j,k) is 1 if Sp{i,j} dep1ends on p(k)
dep1 = zeros([size(lpvux) 2]);% x_3,p(1),x_4
dep1(2,2,:) = [1 0];
dep1(2,3,:) = [1 0]; 
dep1(2,4,:) = [1 1];
dep1(2,5,:) = [1 0];
dep1(4,2,:) = [1 0];
dep1(4,3,:) = [1 0];
dep1(4,4,:) = [1 1];
dep1(4,5,:) = [1 0];

% sampling intervals for each parameter
domain = [-27/180*pi 27/180*pi; -0.45 0.45];
% grid size: number of grid points for each parameter
gridsize= [101,101];

% time for compute ODE45
Tspan=6E-3;
% Cal Itration Times
SimTime=50;
VarSize=4;
IterationTimes=floor(SimTime/Tspan);
X=zeros(VarSize,IterationTimes);
X_PDC=zeros(VarSize,IterationTimes);
time=zeros(1,IterationTimes);
U_ct=zeros(1,IterationTimes);
U_pdc=zeros(1,IterationTimes);% zeros(2,IterationTimes)y有意思了
sigma_t=zeros(1,IterationTimes);
opt = odeset('RelTol',1e-4,'AbsTol',1e-4);
%% cal initial conditon
% sampling
lpvdata1 = sampling_lpv(lpvux, dep1, domain, gridsize);
% hosvd
svtol=1E-3;
keep=[3,2];   % reserved sigular value
hull = 'cno';
[S1 U1 sv tol] = hosvd_lpv(lpvdata1, dep1, gridsize,svtol, keep);
% generating tight polytopic representation
hull = 'cno';
U1 = genhull(U1, hull);
S1 = coretensor(U1, lpvdata1, dep1);

%% plot hull
for ii=1:length(U1)  
    plothull(U1(ii))
    pause(10)
    name=strcat('SPGTPDCshul',num2str(ii));
    epsname1=strcat(name,'.eps' );
    saveas(gcf,epsname1,'epsc2')
end
close all

% check model approximation error
[maxerr meanerr] = tperror(lpvux, S1, U1, domain, 1000);
disp('max and mean error:'); disp(maxerr); disp(meanerr);


% State feedback TP controller design
%lmi = lmistruct(S1, n1);

I = size(S1);
sizes = I(1:end-2);
R = prod(sizes);

S1 = reshape(S1, [R I(end-1) I(end)]);

m = I(end) - n1;
p = I(end-1) - n1;
A = S1(:, 1:n1, 1:n1);
B = S1(:, 1:n1, n1+1:n1+m);
C = S1(:, n1+1:n1+p, 1:n1);
D = S1(:, n1+1:n1+p, n1+1:n1+m);

X = sdpvar(n1+1, n1+1, 'symmetric');
M = cell(1, R);
for r = 1:R
	M{r} = sdpvar(m, n1+1, 'full');
end

% X > 0
lmi.F = set(X > 0, 'poz def');
%SET has been considered obsolete for many years, and the time has come...
%Update your code. http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Commands.set

%lmi = lmi_asym_decay(lmi1, 0.05);  % |x(t)| < c exp(-0.05 t)
R = size(A, 1);
alpha=0.06;
% X*Ar' + Ar*X - Br*Mr - Mr'*Br' + 2*alpha*X < 0
for r = 1:R
	Ar = [reshape(A(r,:,:), [n1 n1]),zeros(n1,1);-[1,0,0,0,0]];
	Br = [reshape(B(r,:,:), [n1 m]);0];
	lmi.F = lmi.F + set(X*Ar' + Ar*X - Br*M{r} - M{r}'*Br' + 2*alpha*X < 0, sprintf('type1 lmi %d', r));
end

% X*Ar' + Ar*X + X*As' + As*X - Br*Ms - Ms'*Br' - Bs*Mr - Mr'*Bs' + 4*alpha*X <= 0
for r = 1:R
	for s = r+1:R
		Ar = [reshape(A(r,:,:), [n1 n1]),zeros(n1,1);-[1,0,0,0,0]];
		As = [reshape(A(s,:,:), [n1 n1]),zeros(n1,1);-[1,0,0,0,0]];
		Br = [reshape(B(r,:,:), [n1 m]);0];
		Bs = [reshape(B(s,:,:), [n1 m]);0];
		lmi.F = lmi.F + set(X*Ar' + Ar*X + X*As' + As*X - Br*M{s} - M{s}'*Br'...
                  - Bs*M{r} - M{r}'*Bs' + 4*alpha*X <= 0, sprintf('type2 lmi %d', r));
	end
end
umax = 8;
phi = 1;
%lmi1 = lmi_input(lmi, umax, phi);  % if |x| < phi then |u| < umax
R = size(A, 1);
% constraints on the control value

% phi^2 I < X
lmi.F = lmi.F + set(phi^2 * eye(n1+1) < X, 'phi^2 I < X');

% [X, Mr'; Mr, mu^2 I] > 0
for r = 1:R
	lmi.F = lmi.F + set([X M{r}'; M{r} umax^2*eye(m)] > 0, sprintf('type3 lmi %d', r));
end
lmi.sizes = sizes;
lmi.n = n1+1;
lmi.m = m;
lmi.p = p;
lmi.A = A;
lmi.B = B;
lmi.C = C;
lmi.D = D;
lmi.X = X;
lmi.M = M;
[K1,P1] = lmi_solve(lmi);
% [K1,P1,Kim] = lmi_solve(lmi);
% Kim=squeeze(Kim);

x1_0=x_0;
opt = odeset('RelTol',1e-4,'AbsTol',1e-5);
IterationTimes=floor(SimTime/Tspan);
X=zeros(4,IterationTimes);
time=X(1,:);
tpdc_u=time;
htpdc_u=time;
%% reference signal y_r
q=(y_r-x1_0(1))*Tspan;% extended signal q
for i=1:IterationTimes
    x1_0
    TPPDC.time(i)=i*Tspan;
    TPPDC.X(:,i)=x1_0;
     % model of cranes
    TPPDC.u(i) = towercraneux_CalUq(x1_0,K1,U1,domain,n1,q);
    [t,x]=ode45(@SPG_plantnoise,[0,Tspan],x1_0,opt,TPPDC.u(i),lpvux,dt);
    q=q+(y_r-x1_0(1))*Tspan;
    x1_0=x(end,:)';
end
Hnamemat=strcat('TPPDC','.mat');
save(Hnamemat,'time','TPPDC')




