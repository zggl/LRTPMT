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
gridsize = [137,137];
lpvdata = sampling_lpv(SPG_LPV1, dep1, domain, gridsize)
% hosvd
[S U sv tol] = hosvd_lpv(lpvdata, dep1, gridsize, 0.001,[2,2]);
sv
% generating tight polytopic representation
hull = 'cno';
U = genhull(U, hull);
S = coretensor(U, lpvdata, dep1);

%plot the results
%plothull(U, domain);
%% plot hull
for ii=1:length(U)  
    plothull(U(ii))
    pause(10)
    name=strcat('SPGTPDCshul',num2str(ii));
    epsname1=strcat(name,'.eps' );
    saveas(gcf,epsname1,'epsc2')
end
close all

% check model approximation error
[maxerr meanerr] = tperror(SPG_LPV1, S, U, domain, 1000);
disp('max and mean error:'); disp(maxerr); disp(meanerr);
% save('tora_data', 'S', 'U', 'n', 'domain', 'gridsize');

% State feedback TP controller design
lmi = lmistruct(S, n);
lmi = lmi_asym_decay(lmi, 0.05);  % |x(t)| < c exp(-0.05 t)
umax = 8;
phi = 1;
lmi = lmi_input(lmi, umax, phi);  % if |x| < phi then |u| < umax
K = lmi_solve(lmi);
x_0=[0.3,0,pi/4,0]';
opt = odeset('RelTol',1e-4,'AbsTol',1e-5);
Tspan=0.01;
SimTime=50;
IterationTimes=floor(SimTime/Tspan);
TPPDC.X=zeros(4,IterationTimes);
TPPDC.time=TPPDC.X(1,:);
for i=1:IterationTimes
    TPPDC.time(i)=i*Tspan;
    TPPDC.X(:,i)=x_0;
     % model of TORA
    u(i) =SPG_CalU(x_0,K,U,domain);
    [t,x]=ode45(@SPG_plant,[0,Tspan],x_0,opt,u(i),SPG_LPV1);
    x_0=x(end,:)';
    TPPDC.u(i)=u(i);
end
%% plot
% figure 1
save('TPPDC.mat','TPPDC')

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

