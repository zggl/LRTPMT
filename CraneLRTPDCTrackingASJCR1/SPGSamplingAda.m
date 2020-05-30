%PDIP LRTPPDC 20191011
% TP model transformation based controller design for the parallel-type double inverted pendulum
% Szabolcs Nagy, Zoltan Petres, Peter Baranyi
% FUZZ-IEEE 2008 p1374-1380
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
SPG_LPV1=@(p)[0,1,0,0,0;...
 0,-(I_p+M_p*l_p^2)*( B_eq+(eta_g *K_g^2*eta_m* K_t *K_m/(R_m *r_mp^2))) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2),...
      M_p^2*l_p^2 *g *cos(p(1))*sinc(p(1)/pi) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2),...
      ((M_p^2*l_p^3+M_p*l_p^2)*sin(p(1))*p(2)+M_p* l_p* B_p*cos(p(1))) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2),...
       -(I_p*M_p *l_p)^2*(eta_g *K_g*eta_m* K_t)/(R_m* r_mp) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2);...
  0,0,0,1,0;...
  0,...
    M_p* l_p*cos(p(1))*( B_eq-(eta_g* K_g^2*eta_m* K_t* K_m) / ((R_m* r_mp^2))*((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2)),...
    -(M_p+M_c)*M_p* l_p*sinc(p(1)/pi) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2),...
   ((-(M_p+M_c)*B_p-M_p^2*l_p^2*cos(p(1))*sin(p(1))*p(2)) ) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2),...
    (-M_p* l_p*cos(p(1))*(eta_g^2* K_g*K_t)/(R_m* r_mp)) / ((M_p+M_c)*I_p+M_c*M_p*l_p^2+ M_p^2*l_p^2*sin(p(1))^2);
]
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

first2dimp=gridsize(1);
second2dimp=gridsize(2);

gp3=sort([0,linspace(domain(1,1),domain(1,2),first2dimp)]);    % x3-adaptive
gp4=sort([linspace(domain(2,1),domain(2,2),second2dimp),0]);    % x4-adaptive

x3size=length(gp3);
x4size=length(gp4);
tic

for m3=1:x3size
    for m4=1:x4size
        x3=gp3(m3);x4=gp4(m4);
        S(m3,m4,:,:)=SPG_LPV1([x3,x4]);
    end
end

toc

%% 
tic
valueflag=[1,1,1,1];
[S,U,sv]=hosvd(S,valueflag,1e-9);

dim=2

[Uc]=convTPmodel(dim,U,'cno');
toc

gi={};
gi{1}=gp3;
gi{2}=gp4;
for i=1:dim
    figure(i);
    plot(gi{i},Uc{i},'LineWidth',3);
    grid on;
    %
    epsname = strcat('TPMTVInputNumModel','Fig',num2str(i), '.eps' );
    saveas(gcf,epsname,'epsc2')
end
pause(5)
close all

tic
hull = 'cno';
U1 = genhull(Uc, hull);
S = coretensor(U1, S, dep1);
lmi = lmistruct(S, n);
lmi = lmi_asym_decay(lmi, 0.05);  % |x(t)| < c exp(-0.05 t)
umax = 8;
phi = 1;
lmi = lmi_input(lmi, umax, phi);  % if |x| < phi then |u| < umax
[K P pres dres] = lmi_solvediag(lmi);
save('SPGLRTPPDC_data', 'S', 'U', 'n', 'domain', 'gridsize','K','P','pres','dres');
K
toc

(toc-tic)/3600
