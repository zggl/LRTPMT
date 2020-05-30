clc
clear all
coli=1;
n=100;
s=2;
if s==2
    min_ranges_p=[-1,-1]  
    max_ranges_p=[1, 1]  
elseif s==3
    min_ranges_p=[-1,-1,-1]  
    max_ranges_p=[1, 1,1]  
end
udflag=0;
[X_scaled,Xij]=UniformDesignWithScale(n,s,coli,min_ranges_p,max_ranges_p,udflag);
plot(X_scaled(:,1),X_scaled(:,2),'r*')
hold on
n=200;
s=2;
[X_scaled]=UniformDesignWithScale(n,s,coli,min_ranges_p,max_ranges_p,udflag);
plot(X_scaled(:,1),X_scaled(:,2),'bo')
legend({'$p=60$','$p=120$'}, 'Interpreter','latex' )
xlabel('$X_1$','Interpreter', 'latex');
ylabel('$X_2$','Interpreter', 'latex');
epsname = strcat('UDdistributePoints', '.eps' );
saveas(gcf,epsname,'epsc2')
pause(3)
epsname = strcat('UDdistributePoints.png');
saveas(gcf,epsname)
close all