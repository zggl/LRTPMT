%% TORA system plant
function [dx] = SPG_plant(t,x,u,lpv,dt)
% p = [x(3) x(4)];
% 
% AB = querylpv(lpv, p);
% 
% % dx = A x + B u
% dx = AB * [x(:); u];
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
    D_1=@(p)((M+m)*J_cm+m*L^2*(M+m*sin(p(1))^2));
     p = [x(3) x(4)];
    %%%%%%%%%%%%%%%%%%%%
   dx= [x(2);
        -(J_cm+m*L^2*B_eqx+eta_gx*K_gx^2 *eta_mx*K_tx*K_m*R_mx*r_mp^2) *x(2)/D_1(p)+...
          m^2*L^2*g*sin(p(1))*cos(p(1))/D_1(p)+...
          ((m^2*L^3+L*m*J_cm)*sin(p(1))* p(2)+m*L*cos(p(1))*B_palpha)*x(4)/D_1(p)+...
          (-(m*L^2+J_cm)^2*eta_gx*K_gx*eta_mx*K_tx)/(D_1(p) *R_mx*r_mp^2)*u;
          x(4)+dt;
          (m*L*cos(p(1))*B_eqx+eta_gx*K_gx^2 *eta_mx*K_tx*K_m*R_mx*r_mp^2)*x(2)/D_1(p)+...
          -(M+m)*m *g*L*sin(p(1))/D_1(p)+...
          -(M+m)*B_palpha-m^2*L^2*sin(p(1))*cos(p(1))* p(2)*x(4)/D_1(p)+...
          (-(m*L*cos(p(1)))^2*eta_gx*K_gx* eta_mx*K_tx)/(D_1(p)* R_mx*r_mp)*u;
   ];