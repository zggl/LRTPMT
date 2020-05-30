clc
clear
close all
%% LPV model definition for 3D crane
B_eq=5.4 ;
B_p=0.0024 ;
eta_g= 1 ;eta_m=1;g=9.81;
I_p=0.0078838 ;J_m=3.9001e-7 ;
K_g=3.71;K_m=0.0076776;K_t=0.007683;
l_p= 0.3302;
M_c=1.0731;M_p=0.23;
R_m=2.6;r_mp=0.00635;
m = 0.32;
%% p(1)-x_3
SPG_LPV1={@(p)0,@(p)1,@(p)0,@(p)0,@(p)0;...
 @(p)0,...
             @(p)-(I_p+M_p*l_p^2)*( B_eq+(eta_g *K_g^2*eta_m* K_t *K_m/(R_m *r_mp^2))) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2),...
             @(p)M_p^2*l_p^2 *g *cos(p(1))*sinc(p(1)/pi) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2),...
             @(p)((M_p^2*l_p^3+M_p*l_p^2)*sin(p(1))*p(2)+M_p* l_p* B_p*cos(p(1))) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2),...
             @(p)-(I_p*M_p *l_p)^2*(eta_g *K_g*eta_m* K_t)/(R_m* r_mp) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2);...
  @(p)0,@(p)0,@(p)0,@(p)1,@(p)0;...
  @(p)0,...
             @(p)M_p* l_p*cos(p(1))*( B_eq-(eta_g* K_g^2*eta_m* K_t* K_m) / ((R_m* r_mp^2))*((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2)),...
             @(p)-(M_p+M_c)*M_p* l_p*sinc(p(1)/pi) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2),...
             @(p)((-(M_p+M_c)*B_p-M_p^2*l_p^2*cos(p(1))*sin(p(1))*p(2)) ) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2),...
             @(p)(-M_p* l_p*cos(p(1))*(eta_g^2* K_g*K_t)/(R_m* r_mp)) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2);
} 
% Parameter relation tensor
dep1 = zeros([size(SPG_LPV1) 2]); % x_3,p(1),x_4

%---------------------------
dep1(2,2,:) = [1 0];
dep1(2,3,:) = [1 0];
dep1(2,4,:) = [1 1];
dep1(2,5,:) = [1 0];

dep1(4,2,:) = [1 0 ];
dep1(4,3,:) = [1 0];
dep1(4,4,:) = [1 1];
dep1(4,5,:) = [1 0];
n = 4;  
% sampling intervals for each parameter
domain = [-27*pi/180,27*pi/180; -0.8,0.8];
gridsize = [137,137];
%%  [Y, X] = minandmax(F, 'local')
%  % 2;1 x_3 
lp1=domain(1,1);rp1=domain(1,2);
lp2=domain(2,1);rp1=domain(2,2);
chebf{2,2,1} = chebfun(@(x) -(I_p+M_p*l_p^2)*( B_eq+(eta_g *K_g^2*eta_m* K_t *K_m/(R_m *r_mp^2))) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(x)^2),[lp1,rp1]);
[ignored,chebfextrema{2,2,1}] = minandmax(chebf{2,2,1},'local');

% 2;4 x_3 
chebf{2,3,1} = chebfun(@(x)M_p^2*l_p^2 *g *cos(x)*sinc(x/pi) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(x)^2),[lp1,rp1]);
[ignored,chebfextrema{2,3,1}] = minandmax(chebf{2,3,1},'local');

% 3;4 x_3 
f = @(x,y)((M_p^2*l_p^3+M_p*l_p^2)*sin(x)*y+M_p* l_p* B_p*cos(x)) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(x)^2);
%  Create a chebfun2
F = chebfun2(f, [domain(1,1),domain(1,2),domain(2,1),domain(2,2)],'vectorize');
[minf,minx] = min2(F);
[maxf,maxx] = max2(F);
chebf{2,4,1}=minx';
chebf{2,4,2}=maxx';
% 2;5 x_3 x_4-----
 chebf{2,5,1} = chebfun(@(x)-(I_p*M_p *l_p)^2*(eta_g *K_g*eta_m* K_t)/(R_m* r_mp) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(x)^2),[lp1,rp1]);
 [ignored,chebfextrema{2,5,1}] = minandmax(chebf{2,5,1},'local');
% 4;2
chebf{4,2,1} = chebfun(@(x)M_p* l_p*cos(x)*( B_eq-(eta_g* K_g^2*eta_m* K_t* K_m) / ((R_m* r_mp^2))*((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(x)^2)),[lp1,rp1]);
[ignored,chebfextrema{4,2,1}] = minandmax(chebf{4,2,1},'local');
% % 4;3 x_3 
chebf{4,3,1} = chebfun(@(x)-(M_p+M_c)*M_p* l_p*sinc(x/pi) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(x)^2),[lp1,rp1]);
[ignored,chebfextrema{4,3,1}] = minandmax(chebf{4,3,1},'local');
% % 4;4 x_3 x_4
f = @(x,y)((-(M_p+M_c)*B_p-M_p^2*l_p^2*cos(x)*sin(x)*y) ) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(x)^2);
%  Create a chebfun2
F = chebfun2(f, [domain(1,1),domain(1,2),domain(2,1),domain(2,2)],'vectorize');
[minf,minx] = min2(F);
[maxf,maxx] = max2(F);
chebf{4,4,1}=minx';
chebf{4,4,2}=maxx';
% % 4;5 x_3 
chebf{4,5,1} = chebfun(@(x)(-M_p* l_p*cos(x)*(eta_g^2* K_g*K_t)/(R_m* r_mp)) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(x)^2),[lp1,rp1]);
[ignored,chebfextrema{4,5,1}] = minandmax(chebf{4,5,1},'local');
%% TP transformation, same as:
%   [S U] = tptrans(lpv, dep, domain, gridsize, 'close');
TP.lpv=SPG_LPV1;
TP.dep=dep1;
TP.domain =domain;
TP.gridsize=gridsize;
TP.siz=prod(gridsize);
%% HS parameters
p = primes(410);    %  prime base
pbaseindx=3;        %  prime base vector index
p=p(pbaseindx:end); % p_1=5
Datanum=TP.siz;
Xij= DIST_HAMMERSLEYDiy(domain, Datanum,p);
TP.X_scaled=Xij;
lpvdata = sampling_lpvud(TP);

% hosvd
[S U sv tol] = hosvd_lpv(lpvdata, dep1, gridsize, 0.001,[3,2]);

% generating tight polytopic representation
hull = 'cno';
U = genhull(U, hull);
S = coretensor(U, lpvdata, dep1);

%% MVS
ii=1
[UsizeM,UsizeN]=size(U{ii});
U1=[U{ii},ones(UsizeM,1)];

[V,Sigm,Bt]=svd(U1,'econ'); %A=USV' econ
T0=Sigm*Bt'*[eye(UsizeN);zeros(1,UsizeN)];% Sigma B'E
J=size(V,2);
[Ae, indice, Rp ]= VCA(V,'Endmembers',J,'SNR',0,'verbose',1);%VCA to guess R0

R=V(indice,:);
figure(1)
for ij=1:size(R,2)-1
    subplot(2,2,ij)
    plot(V(:,ij),V(:,ij+1),'o')
    hold on 
    plot(R(:,ij),R(:,ij+1),'*')
end
R0=V(indice,:)+randi(J,1)*(V(indice,:)-sum(V(indice,:))/J);
pause(2)
close all
[M1,Up,Q,my,sing_values] = mvsa(V',UsizeN+1,'M0',0,'spherize','yes'); %  Minimum Volume Simplex Analysis (MVSA),
%Up'*Up=I_N, U{ii}'*U{ii}=I_N£¬M=Up/Q=>M*Q=Up
% norm(M-Up/Q,'fro')
S1=inv(M1)*V';% V'=MS->V=S'M'
mean(sum(S1),2)% equiv 1,mean the right weight
norm(S1'*M1'-V,1)% W=S',Rinv=M'.
W1=S1';
%sum(W1,2)
T1=M1'*T0; %60x3
norm(W1*T1-U{ii},'fro')% U=WT

ii=2
[UsizeM,UsizeN]=size(U{ii});
U1=[U{ii},ones(UsizeM,1)];
[V,Sigm,Bt]=svd(U1,'econ'); %A=USV' econ
J=size(V,2);
T0=Sigm*Bt'*[eye(UsizeN);zeros(1,UsizeN)];
[Ae, indice, Rp ]= VCA(V,'Endmembers',J,'SNR',0,'verbose',1);%VCA to guess R0
R=V(indice,:);
R0=V(indice,:)+randi(J,1)*(V(indice,:)-sum(V(indice,:))/J);
[M2,Up,Q,my,sing_values] = mvsa(V',UsizeN+1,'M0',0,'spherize','yes'); %  Minimum Volume Simplex Analysis (MVSA),
%Up'*Up=I_N, U{ii}'*U{ii}=I_N£¬M=Up/Q=>M*Q=Up
% norm(M-Up/Q,'fro')

S2=inv(M2)*V';% V'=MS->V=S'M'
mean(sum(S2),2)% equiv 1,mean the right weight
norm(S2'*M2'-V,1)% W=S',Rinv=M'.
W2=S2';
%sum(W2,2)
T2=M2'*T0; %60x3
norm(W2*T2-U{ii},'fro')% U=WT

% ii=2 could not make M feasible
S_mvs = tproddimi(S, T1, 1); % S \times _1 T
S_mvs = tproddimi(S_mvs, T2, 2); % S \times _1 T
U{1}=W1;
U{2}=W2;

hull = 'cno';
U1 = genhull(U, hull);
lpvux=SPG_LPV1;
[maxerr meanerr] = tperror(lpvux, S_mvs, U1, domain, 100);
alpha=0.06;
umax=8;
phi=1;
n1=n;
[K1,P1,Kim,diagnostic]=LMIsolve(S_mvs,U1,n1,m,alpha,umax,phi);
pause(3)

Kim=squeeze(Kim);
x_0=[0.3,0,pi/4,0]';
opt = odeset('RelTol',1e-4,'AbsTol',1e-5);
Tspan=6E-3;
SimTime=50;
IterationTimes=floor(SimTime/Tspan);
MVSLRHSTPPDC.X=zeros(4,IterationTimes);
MVSLRHSTPPDC.time=MVSLRHSTPPDC.X(1,:);
%%y_r
y_r=0.4;
%% cal initial conditon
x_0=[0.3,0,pi/4,0]';
x1_0=x_0;
dt=0.1;
q=(y_r-x_0(1))*Tspan;% extended signal q
for i=1:IterationTimes
    x_0
    MVSLRHSTPPDC.time(i)=i*Tspan;
    MVSLRHSTPPDC.X(:,i)=x_0;
     % model of cranes
    MVSLRHSTPPDC.u(i) = towercraneux_CalUq(x_0,K1,U1,domain,n1,q);
    [t,x]=ode45(@SPG_plantnoise,[0,Tspan],x_0,opt,MVSLRHSTPPDC.u(i),lpvux,dt);
    q=q+(y_r-x_0(1))*Tspan;
    x_0=x(end,:)';
end
save('MVSLRHSTPPDC.mat','MVSLRHSTPPDC')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K1,P1,Kim,diagnostic]=LMIsolve(S1,U1,n1,m,alpha,umax,phi)
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

    %lmi = lmi_asym_decay(lmi1, 0.05);  % |x(t)| < c exp(-0.05 t)
    R = size(A, 1);
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

    [K1,P1,Kim,diagnostic] = lmi_solveplus(lmi);
    % best value of t: 
    %diagnostic.solveroutput.copt
end

