close all
H2O2 = [0.005,0.1,0.25,0.5,1,2];
OPD = 0.5;
ForFit = [0.7 1];
for k = 1:6
    %figure
    file_name = ['DiffBoth_Result_H2O2',num2str(H2O2(k)*1e3),'uM_OPD_Conc_',num2str(OPD*1e3),'Complicated.mat'];
    load(file_name);
    eval(['DAP_H2O2_',num2str(H2O2(k)*1e3),'uMLimit1e3 = DAP_Coa;']);
    DAP_row_trace = squeeze(DAP_Coa(15,:,:))./2;
    DAP_col_trace = squeeze(DAP_Coa(:,7,:))./2;
    subplot(2,3,k)
    plot(time./60,DAP_col_trace)
    %legend('RD = 0.22','RD = 0.11','RD = 0','RD = -0.11','RD = -0.22')
    title(['H2O2 = ',num2str(H2O2(k)*1e3),'uM Col'])
    %xlim([0 3000]);
    %ylim([0 0.2])
    %subplot(1,2,i-2)
    %figure;
    %imagesc(time./60,15:1:1,DAP_col_trace);
    %caxis([0 0.6]);
    %ylabel('Relative positions')
    %xlabel('Time')
    %colorbar
    %legend('RD = 0.44','RD = 0.33','RD = 0.22','RD = 0.11','RD = 0')
    %legend('RD = 0.22','RD = 0.11','RD = 0','RD = -0.11','RD = -0.22')
    %xlim([5000 10000]);
    %ylim([0 0.3]);
    %title(['H2O2=',num2str(H2O2(i)*1e3),'uM DAP col'])
    %subplot(1,2,i-1)
    %figure;
    %imagesc(time./60,15:1:1,DAP_row_trace);
    %caxis([0 0.6]);
    %ylabel('Relative positions')
    %xlabel('Time')
    %colorbar
   % xlim([5000 10000]);
    %ylim([0 0.3]);
    %legend('RD = 0.44','RD = 0.33','RD = 0.22','RD = 0.11','RD = 0')
    %legend('RD = 0.22','RD = 0.11','RD = 0','RD = -0.11','RD = -0.22')
    %title(['H2O2=',num2str(H2O2(i)*1e3),'uM DAP row'])
    
    for j = 1:length(DAP_col_trace(:,1))
        Time_Point1e3(k,j) = interp1(DAP_col_trace(j,100:end),time(100:end),0.1,'linear');
        %[M,I_1] = min(abs(DAP_col_trace(j,:)-0.2)); %[M,I_2] = min(abs(DAP_col_trace(j,:)-0.4));
        time_ForFit = time(100):0.1:time(end);
        DAP_Col_ForFit = interp1(time(100:end),DAP_col_trace(j,100:end),time_ForFit,'linear');
        DAP_Row_ForFit = interp1(time(100:end),DAP_row_trace(j,100:end),time_ForFit,'linear');
        [M,I_1] = min(abs(DAP_Col_ForFit-DAP_Col_ForFit(end)*ForFit(1))); [M,I_2] = min(abs(DAP_Col_ForFit-DAP_Col_ForFit(end)*ForFit(2)));
        p1(j,:,k) = polyfit(time_ForFit(I_1:I_2),DAP_Col_ForFit(I_1:I_2),1);
        [M,I_1] = min(abs(DAP_Row_ForFit-DAP_Row_ForFit(end)*ForFit(1))); [M,I_2] = min(abs(DAP_Row_ForFit-DAP_Row_ForFit(end)*ForFit(2)));
        p2(j,:,k) = polyfit(time_ForFit(I_1:I_2),DAP_Row_ForFit(I_1:I_2),1);
        Onset_Time_col(j,k) = -p1(j,2,k)/p1(j,1,k);
        Onset_Time_row(j,k) = -p2(j,2,k)/p2(j,1,k);
    end
    
end

%%
Diff_Time = [Time_Point1e3(:,1:end-1) - Time_Point1e3(:,2:end)]';
Diff_Time_row = Onset_Time_row(1:end-1,:) - Onset_Time_row(2:end,:);
Diff_Time2 = Onset_Time_col(1:end-1,:) - Onset_Time_col(2:end,:);
Diff_Time2_2 = Diff_Time2(1:end-1,:) - Diff_Time2(2:end,:);
figure;
for k = 1:14
    subplot(3,6,k)
    plot([0.005 0.125 0.25 0.5 1 2], Diff_Time2(k,[1 2 3 4 5 6]),'o-');
end
%legend('1','2','3','4','5');
hold off;
figure;
plot(Diff_Time2(:,[1 2 3 4 5 6]),'o-');
%%
figure;
PlotData = flip(Conc_record_1_col(1:end-1,:));
time_interval = [1 226];
Fac = 4;
subplot(2,2,1)
h = imagesc(time(time_interval(1):time_interval(2))./60,1:15,PlotData(:,time_interval(1):time_interval(2)));
set(gca,'YTick',[0:2:15]);
set(gca,'YTickLabels',{'0','110','220','330','440','550','660','770'});
%xlim([0 451])
title('H_2O_2 concentration evolution in column')
xlabel('time (min)'); 
ylabel('Relative position (\mum)');
colorbar 
colormap jet
caxis([0 max(max(PlotData(:,time_interval(1):time_interval(2))))/Fac])

PlotData = flip(Conc_record_1_row(1:end-1,:));
subplot(2,2,2)
h = imagesc(time(time_interval(1):time_interval(2))./60,1:15,PlotData(:,time_interval(1):time_interval(2)));
set(gca,'YTick',[0:2:15]);
set(gca,'YTickLabels',{'0','110','220','330','440','550','660','770'});
%xlim([0 451])
title('H_2O_2 concentration evolution in row')
xlabel('time (min)'); 
ylabel('Relative position (\mum)');
colorbar 
colormap jet
caxis([0 max(max(PlotData(:,time_interval(1):time_interval(2))))/Fac])

PlotData = flip(Conc_record_2_col(1:end-1,:));
subplot(2,2,3)
h = imagesc(time(time_interval(1):time_interval(2))./60,1:15,PlotData(:,time_interval(1):time_interval(2)));
set(gca,'YTick',[0:2:15]);
set(gca,'YTickLabels',{'0','110','220','330','440','550','660','770'});
%xlim([0 451])
title('o-PD concentration evolution in column')
xlabel('time (min)'); ylabel('Relative position (\mum)');
colorbar 
colormap jet
caxis([0 max(max(PlotData(:,time_interval(1):time_interval(2))))/Fac])

PlotData = flip(Conc_record_2_row(1:end-1,:));
subplot(2,2,4)
h = imagesc(time(time_interval(1):time_interval(2))./60,1:15,PlotData(:,time_interval(1):time_interval(2)));
set(gca,'YTick',[0:2:15]);
set(gca,'YTickLabels',{'0','110','220','330','440','550','660','770'});
%xlim([0 451])
title('o-PD concentration evolution in column')
xlabel('time (min)'); ylabel('Relative position (\mum)');
colorbar 
colormap jet
caxis([0 max(max(PlotData(:,time_interval(1):time_interval(2))))/Fac])
%%
% Diff_Desc_array = [1,0;0,1;0,0];
% offset = 300;
% for k = 1:3
%     H2O2_Conc = 0.5; %mM
%     OPD_Conc = 0.5; %mM
%     Diff_Desc = Diff_Desc_array(k,:);
%     file_name = ['DiffDec',num2str(Diff_Desc(1)),num2str(Diff_Desc(2)),'_Result_H2O2',num2str(H2O2_Conc*1e3),'uM_OPD_Conc_',num2str(OPD_Conc*1e3),'Complicated'];
%     load(file_name)
%     eval(['DAP_H2O2_',num2str(H2O2_Conc*1e3),'uMLimit1e3 = DAP_Coa;']);
%     DAP_row_trace = squeeze(DAP_Coa(7,:,:))./2;
%     DAP_col_trace = squeeze(DAP_Coa(:,7,:))./2;
%     for p = 1:5
%         DAP_col_Acc(:,p,k) = sum(DAP_col_trace([15 14 13 12 11 10 9],1+113*(p-1)+offset:113*p+offset).*8./60,2);
%     end
%     subplot(1,3,k)
%     bar(DAP_col_Acc(:,:,k),'stack')
%     %legend('RD = 0.22','RD = 0.11','RD = 0','RD = -0.11','RD = -0.22')
%     title(['DiffDec',num2str(Diff_Desc(1)),num2str(Diff_Desc(2))])
% end
%%
