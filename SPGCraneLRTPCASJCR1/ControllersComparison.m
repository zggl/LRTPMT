
clc
clear
filenam={'TPPDC.mat','LRTPPDC.mat','MVSLRHSTPPDC.mat','MVSLRUDTPPDC.mat'} %9
filen=length(filenam);
for filei=1:filen
    load(filenam{filei})
end
linsty={'--',':','-.','-'};
markesty={'*','o','v','d','x','s','^','<','>','p','h'};
indexi=1:4;
indxstep=600;
for figi=indexi
    figure(figi) 
    index1=figi;
    plot(TPPDC.time,TPPDC.X(index1,:),'-*',...
                LRTPPDC.time,LRTPPDC.X(index1,:),'--o',...
                MVSLRHSTPPDC.time,MVSLRHSTPPDC.X(index1,:),'-x',...
                MVSLRUDTPPDC.time,MVSLRUDTPPDC.X(index1,:),'-.s',...
               'MarkerIndices',1:indxstep:size(TPPDC.time,2),'LineWidth',1.2);
    if ((figi==1) || (figi==2))
    legend('TPDC','LRTPDC','MLHTPPDC','MLUTPDC')        
    else
       legend('TPDC','LRTPDC','MLHTPDC','MLUTPDC','Location','southeast')     
    end
    ylabel(strcat('x_',num2str(figi),'(t)'))
    xlabel('Time (s)');
    epsname1=strcat('SPGLRTPDCsx',num2str(figi),'.eps');
    saveas(gcf,epsname1,'epsc2')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figi+1)
plot(TPPDC.time,TPPDC.u,'-*',...
        LRTPPDC.time,LRTPPDC.u,'--o',...
        MVSLRHSTPPDC.time,MVSLRHSTPPDC.u,'-h',...
        MVSLRUDTPPDC.time,MVSLRUDTPPDC.u,'-.s',...
        'MarkerIndices',1:indxstep:size(TPPDC.time,2),'LineWidth',1.2);
legend('TPDC','LRTPDC','MLHTPDC','MLUTPDC')
xlabel('Time (s)');
ylabel('u')
name='SPGLRTPDCsu';
epsname1=strcat(name,'.eps' );
saveas(gcf,epsname1,'epsc2')

pause(20)
close all
