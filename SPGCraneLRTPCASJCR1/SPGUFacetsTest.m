%% LPV model definition for SPG
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
lpv={@(p)0,@(p)1,@(p)0,@(p)0,@(p)0;...
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
dep1 = zeros([size(lpv) 2]); % x_3,p(1),x_4

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
%gridsize = [137,137];
% grid size: number of grid points for each parameter
TP.UDn=[60];
iind=0;
for ii=TP.UDn
    iind=iind+1;
    gridsize = [ii,ii];
    %% TP transformation, same as:
    %   [S U] = tptrans(lpv, dep, domain, gridsize, 'close');
    TP.lpv=lpv;
    TP.dep=dep1;
    TP.domain =domain;
    TP.gridsize=gridsize;
    TP.siz=prod(gridsize);
    %% UD parameters

    TP.coli=1;
    TP.s=2;
    TP.UDflag=1;
    if TP.s==2
        min_ranges_p=domain(:,1); 
        max_ranges_p=domain(:,2);
    elseif TP.s==3
        min_ranges_p=domain(:,1); 
        max_ranges_p=domain(:,2); 
    end
    %[TP.X_scaled,Xij]=UniformDesignWithScale(TP.siz,TP.s,TP.coli,min_ranges_p,max_ranges_p);
    datafile=strcat('UDData',num2str(TP.siz),'.mat')
    load(datafile)
    % sampling via AUDTP
    % lpvdata = sampling_lpvud(TP);
    % sampling via TP
    TP.X_scaled=Xij;
    lpvdata = sampling_lpvud(TP);

    % hosvd
    [S U sv tol] = hosvd_lpv(lpvdata, TP.dep, gridsize, 0.01,[3,2]);
end
[m,n]=size(U{1});
nck=nchoosek(m,n)
[matind,facets]=UijFacets(U{1});



