clc
clear
close all
%% LPV model definition for 3D crane
%% LPV model definition for quadrotor
B_eq=5.4 ;
B_p=0.0024 ;
eta_g= 1 ;eta_m=1;g=9.81;
I_p=0.0078838 ;J_m=3.9001e-7 ;
K_g=3.71;K_m=0.0076776;K_t=0.007683;
l_p= 0.3302;
M_c=1.0731;M_p=0.23;
R_m=2.6;r_mp=0.00635;
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
filename=strcat('UDData',num2str(prod(gridsize)),'.mat');
load(filename);
TP.X_scaled=Xij;
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
%% UD parameters

TP.coli=1;
TP.s=1;
TP.UDflag=1;
if TP.s==2
    min_ranges_p=domain(:,1); 
    max_ranges_p=domain(:,2);
elseif TP.s==3
    min_ranges_p=domain(:,1); 
    max_ranges_p=domain(:,2); 
end
lpvdata = sampling_lpvud(TP);
% hosvd
[S U sv tol] = hosvd_lpv(lpvdata, dep1, gridsize, 0.001,[2,2]);
sv
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

% generating tight polytopic representation      
lmi = lmistruct(S_mvs, n);
lmi = lmi_asym_decay(lmi, 0.05);  % |x(t)| < c exp(-0.05 t)
umax = 8;
phi = 1;
lmi = lmi_input(lmi, umax, phi);  % if |x| < phi then |u| < umax
K = lmi_solve(lmi);
% time for compute ODE45
Tspan=6E-3;
% Cal Itration Times
SimTime=50;
VarSize=4;
ItrationTimes=floor(SimTime/Tspan);
UDTPPDC.X=zeros(VarSize,ItrationTimes);
UDTPPDC.time=zeros(1,ItrationTimes);
UDTPPDC.U_pdc=zeros(2,ItrationTimes);
sigma_t=zeros(2,ItrationTimes);
opt = odeset('RelTol',1e-4,'AbsTol',1e-4);
% error is 2 demension, $[e,\dot{e}]$

%% cal initial conditon
x_0=[0.3,0,pi/4,0]';
x1_0=x_0;
index1=[1,3];
index2=[2,4];
flag=3;
Psi1=1E-1;
alpha1=1;
%% simu 
for i=1:ItrationTimes
    x1_0
    p =[x1_0(3), x1_0(4)];
    refx=0.185;
    p_pdc = [x1_0(3), x1_0(4)];
    W = queryw1(U,domain,p);
    z = shiftdim(tprod(K,W));
    u_pdc = tpcontroller(p_pdc, x1_0, K, U, domain); 
    [t1,x1]=ode45(@SPG_plant,[0,Tspan],x1_0,opt,u_pdc,SPG_LPV1); 
    %fun1(t),a new type help ode45
    MVSLRUDTPPDC.time(i) = Tspan * i;
    MVSLRUDTPPDC.X(:,i) = x1(end,:)';
    x1_0 = x1(end,:)';
    MVSLRUDTPPDC.u(i)=u_pdc;
end
save('MVSLRUDTPPDC.mat','MVSLRUDTPPDC')
% add the start time [0;0]
% figure(1) 
% plot(time,X([index1],:),'-',time,X_PDC([index1],:),'-.','LineWidth',1.2);
% legend('TPVUD:x_1(t)','TPVUD:x_3(t)','TPPDC:x_1(t)','TPPDC:x_3(t)')
% ylabel('x_1(t),x_2(t)')
% xlabel('Time (s)');
% name='SGCUDTPDCx1x3';
% epsname1=strcat(name,'.eps' );
% saveas(gcf,epsname1,'epsc2')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% plot(time,X([index2],:),'-',time,X_PDC([index2],:),'-.','LineWidth',1.2);
% legend('TPVUD:x_2(t)','TPVUD:x_4(t)','TPPDC:x_2(t)','TPPDC:x_4(t)')
% ylabel('x_2(t),x_4(t)')
% xlabel('Time (s)');
% name='SGCUDTPDCx2x4';
% epsname1=strcat(name,'.eps' );
% saveas(gcf,epsname1,'epsc2')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% plot(time,U_ct,'-',time,U_pdc,'-.','LineWidth',1.2);
% legend('TPVUD:u(t)','TPPDC:u(t)')
% xlabel('Time (s)');
% ylabel('u')
% name='SGCUDTPDCu';
% epsname1=strcat(name,'.eps' );
% saveas(gcf,epsname1,'epsc2')
% %set(gca,'XTick',0:SimTime);
% pause(5)
% close all



