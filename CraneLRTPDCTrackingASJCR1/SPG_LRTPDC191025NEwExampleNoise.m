% function ToraSim_NoSimulink()
% LPV model

% state vector:
%	x1: position of the cart
%	x2: velocity of the cart
%	x3: angular position of the proof body
%	x4: angular velocity of the proof body

% parameters of the LPV model:
%	p1: x3 (angular position)
%	p2: x4 (angular velocity)
clc
clear
%% LPV model definition for quadrotor
B_eq=5.4 ;
B_p=0.0024 
eta_g= 1 ;eta_m=1;g=9.81
I_p=0.0078838 ;J_m=3.9001e-7 
K_g=3.71;K_m=0.0076776;K_t=0.007683
l_p= 0.3302
M_c=1.0731;M_p=0.23
R_m=2.6;r_mp=0.00635
m = 0.32;
%%
%   a_1=-(I_p+M_p*l_p^2)*( B_eq+(eta_g *K_g^2*eta_m* K_t *K_m/(R_m *r_mp^2))),
%   a_2=M_p^2*l_p^2 *g *cos(p(1))*sin(p(1))/p(1),
%   a_3=((M_p^2*l_p^3+M_p*l_p^2)*sin(p(1))*p(2)+M_p* l_p* B_p*cos(p(1)))
%   a_4=M_p* l_p*cos(p(1))*( B_eq-(eta_g* K_g^2*eta_m* K_t* K_m)/(R_m* r_mp^2))
%   a_5=-(M_p+M_c)*M_p* l_p*sin(p(1))/p(1)
%   a_6=(-(M_p+M_c)*B_p-M_p^2*l_p^2*cos(p(1))*sin(p(1))*p(2))
%   a_x=((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2)
%   b_1=-(I_p*M_p *l_p)^2*(eta_g *K_g*eta_m* K_t)/(R_m* r_mp)
%   b_2=(-M_p* l_p*cos(p(1))*(eta_g^2* K_g*K_t)/(R_m* r_mp))
% system matrix: lpv = [A(p) B(p)]
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
gridsize = [136,136];
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

TP.lpv=SPG_LPV1;
TP.dep=dep1;
TP.domain =domain;
TP.gridsize=gridsize;
TP.siz=prod(gridsize);
TP.chebfextrema=chebfextrema;

X1=linspace(domain(1,1),domain(1,2),gridsize(1));
X2=linspace(domain(2,1),domain(2,2),gridsize(2));
[X2,X1] = ndgrid(X2,X1);
combins=[X1(:) X2(:)]; 
TP.X_scaled=combins;

% hosvd
lpvdata = sampling_lpvud(TP);
% hosvd
[S1 U1 sv tol] = hosvd_lpv(lpvdata, dep1, gridsize, 0.001,[3,2]);
% generating tight polytopic representation
hull = 'cno';
U1 = genhull(U1, hull);
S1 = coretensor(U1, lpvdata, dep1);

%% plot hull
for ii=1:length(U1)  
    plothull(U1(ii))
    pause(10)
    name=strcat('SPGLRTPDCshul',num2str(ii));
    epsname1=strcat(name,'.eps' );
    saveas(gcf,epsname1,'epsc2')
end
close all


lpvux=SPG_LPV1;
[maxerr meanerr] = tperror(lpvux, S1, U1, domain, 100);
alpha=0.06;
umax=8;
phi=1;
n1=n;
[K1,P1,Kim,diagnostic]=LMIsolve(S1,U1,n1,m,alpha,umax,phi);
pause(5)
Kim=squeeze(Kim);
x_0=[0.3,0,pi/4,0]';
opt = odeset('RelTol',1e-4,'AbsTol',1e-5);
Tspan=6E-3;
SimTime=50;
IterationTimes=floor(SimTime/Tspan);
LRHSTPPDC.X=zeros(4,IterationTimes);
LRHSTPPDC.time=LRHSTPPDC.X(1,:);
%%y_r
y_r=0.4;
%% cal initial conditon
x_0=[0.3,0,pi/4,0]';
x1_0=x_0;
dt=0.1;
q=(y_r-x_0(1))*Tspan;% extended signal q
for i=1:IterationTimes
    x_0
    LRTPPDC.time(i)=i*Tspan;
    LRTPPDC.X(:,i)=x_0;
     % model of cranes
    LRTPPDC.u(i) = towercraneux_CalUq(x_0,K1,U1,domain,n1,q);
    [t,x]=ode45(@SPG_plantnoise,[0,Tspan],x_0,opt,LRTPPDC.u(i),lpvux,dt);
    q=q+(y_r-x_0(1))*Tspan;
    x_0=x(end,:)';
end
save('LRTPPDC136pointscheb.mat','LRTPPDC')
%% plot
% figure 1


% figure(1)
% %pos_screen=get(0,'screensize');
% %fig_pos=linspace(20,pos_screen(3),4);
% %figure('position',[fig_pos(1),10,600,860])
% subplot(2,1,1)
% plot(time,X(1,:),'--');
% xlabel('time (sec)');ylabel('x_1 (m)');
% % legend({'$q_1$','$q_{d1}$','$\hat{q}_1$'},'Interpreter','latex');
% subplot(2,1,2)
% plot(time,X(2,:),'--')
% xlabel('time (sec)');ylabel('x_2 (m)');
% % legend({'$q_2$','$q_{d2}$','$\hat{q}_2$'},'Interpreter','latex');
% name='SPGHSTPDCx1x2';
% epsname1=strcat(name,'.eps' );
% saveas(gcf,epsname1,'epsc2')
% 
% figure(2)
% subplot(2,1,1)
% plot(time,X(3,:),'--');
% xlabel('time (sec)');ylabel('x_3 (m)');
% % legend({'$q_1$','$q_{d1}$','$\hat{q}_1$'},'Interpreter','latex');
% 
% subplot(2,1,2)
% plot(time,X(4,:),'--')
% xlabel('time (sec)');ylabel('x_4 (m)');
% % legend({'$q_2$','$q_{d2}$','$\hat{q}_2$'},'Interpreter','latex');
% name='SPGHSTPDCx3x4';
% epsname1=strcat(name,'.eps' );
% saveas(gcf,epsname1,'epsc2')
% figure(3)
% plot(time,u,'--')
% xlabel('time (sec)');ylabel('u(t)');
% name='SPGHSTPDCu';
% epsname1=strcat(name,'.eps' );
% saveas(gcf,epsname1,'epsc2')
% close all

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