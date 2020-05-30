clc
clear
filenam={'TPPDC136points.mat','LRTPPDC136pointscheb.mat'} %8
filen=length(filenam);
for filei=1:filen
    load(filenam{filei})
end
linsty={'--',':','-.','-'};
markesty={'*','o','v','d','x','s','^','<','>','p','h'};
indexi=1:4;
indxstep=500;
for figi=indexi
    figure(figi) 
    index1=figi;
    plot(TPPDC.time,TPPDC.X(index1,:),'-',...
                LRTPPDC.time,LRTPPDC.X(index1,:),'-.',...
                'LineWidth',1.2);
    if figi==1
        hold on
        plot(TPPDC.time,0.4*ones(size(TPPDC.X(index1,:))),'-');
    end
    if figi==3
       legend('TPDC','LRTPDC')        
    else
       legend('TPDC','LRTPDC','Location','southeast')     
    end
    ylabel(strcat('x_',num2str(figi),'(t)'))
    xlabel('Time (s)');
    epsname1=strcat('SPGLRTPDCs136Trackingx',num2str(figi),'.eps');
    saveas(gcf,epsname1,'epsc2')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figi+1)
plot(TPPDC.time,TPPDC.u,'-',...
        LRTPPDC.time,LRTPPDC.u,'--','LineWidth',1.2);
legend('TPDC','LRTPDC','Location','southeast')
xlabel('Time (s)');
ylabel('u')
name='SPGLRTPDCs136Trackingu';
epsname1=strcat(name,'.eps' );
saveas(gcf,epsname1,'epsc2')

 pause(60)
 close all