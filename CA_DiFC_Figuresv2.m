%Script for creating Figures for CA DiFC paper
%Josh Pace 20250322


%Place data file 'Contrast Agent and DiFC Data.mat' (download from Pennsieve repository  DOI: 10.26275/btt6-puqw)

load('/Users/jspace/Library/CloudStorage/OneDrive-SharedLibraries-NortheasternUniversity/Niedre, Mark - Niedre_Lab/Josh/Papers and Conferences/CA and DiFC/Figures/Data/Contrast Agent and DiFC Data.mat');

%%
[noise5, background5] = getBackgroundAndNoise(Nude_No_CA.data_bs,5);

figure(1);
tiledlayout(2,1)
nexttile(1)
plot(prelabeled_leg.time,prelabeled_leg.data_bs(:,2));
nexttile(2)
plot(prelabeled_leg.time,prelabeled_leg.data_bs(:,2));
hold on
plot(Nude_No_CA.time,background5(:,2));

%%
%Plotting no injection, 24Hrs post 10nmol OTL38 and VGT-309 background
xlimMax = 1800;
figure(2)
tiledlayout(3,2)
nexttile(1,[3 1])
plot(Nude_No_CA.time,Nude_No_CA.data(:,2))
hold on
plot(Nude_24Hrs_10nmol_OTL38.time,-Nude_24Hrs_10nmol_OTL38.data(:,2))
plot(Nude_24Hrs_10nmol_VGT.time,Nude_24Hrs_10nmol_VGT.data(:,2))
hold off
ylabel('NIR-DiFC Signal (nA)')
ylim([0 50000])
yticks(0:5000:50000)
xlim([0 xlimMax])
%xticks(0:900:3600)
%xticklabels({'0','','30','','60'})
xticks(0:300:3600)
%xticklabels({'0','','30','','60'})
xticklabels({'0','','10','','20','','30'})
xlabel('Time (min)')
% ylim([ 0 14])
% yticks(0:1:14)
% yticklabels({'0','','2','','4','','6','','8','','10','','12','','14'})
% title('24 Hrs Post VGT-309 Injection')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
% nexttile(2)
% plot(Nude_24Hrs_10nmol_OTL38.time,Nude_24Hrs_10nmol_OTL38.data_bs(:,2))
% hold on
% plot(Nude_24Hrs_10nmol_OTL38.time, Nude_24Hrs_10nmol_OTL38.thresh_curve(:,2), '-k','LineWidth',1)
% hold off
% legend({'','Detection Threshold'},'Location','northoutside','Orientation','horizontal')
% 
% %ylabel('NIR-DiFC Signal (nA)')
% ylim([-350 1050])
% yticks(-350:350:1050)
% xlim([0 xlimMax])
% xticks(0:300:3600)
% xticklabels({'0','','10','','20','','30'})
% %xlabel('Time (min)')
% axis = gca;
% axis.FontSize =30;
% axis.FontWeight = 'normal';
% axis.FontName = 'Arial';
% axis.LineWidth = 1;
% nexttile(4)
% plot(Nude_24Hrs_10nmol_VGT.time,Nude_24Hrs_10nmol_VGT.data_bs(:,2))
% hold on
% plot(Nude_24Hrs_10nmol_VGT.time, Nude_24Hrs_10nmol_VGT.thresh_curve(:,2), '-k','LineWidth',1)
% hold off
% ylabel('NIR-DiFC Signal (nA)')
% xlim([0 xlimMax])
% xticks(0:300:3600)
% ylim([-350 1050])
% yticks(-350:350:1050)
% xticklabels({'0','','10','','20','','30'})
% %xlabel('Time(min)')
% axis = gca;
% axis.FontSize =30;
% axis.FontWeight = 'normal';
% axis.FontName = 'Arial';
% axis.LineWidth = 1;
% nexttile(6)
% plot(Nude_No_CA.time,Nude_No_CA.data_bs(:,2))
% hold on
% plot(Nude_No_CA.time, Nude_No_CA.thresh_curve(:,2), '-k','LineWidth',1)
% hold off
% %ylabel('NIR-DiFC Signal (nA)')
% ylim([-350 1050])
% yticks(-350:350:1050)
% xlim([0 xlimMax])
% xticks(0:300:3600)
% xticklabels({'0','','10','','20','','30'})
% xlabel('Time (min)')
% axis = gca;
% axis.FontSize =30;
% axis.FontWeight = 'normal';
% axis.FontName = 'Arial';
% axis.LineWidth = 1;

%%
%Plotting no injection, 24Hrs post 10nmol OTL38 and VGT-309 background, with peak examples
samples_per_sec = 2000;
start_time = 1800;
end_time = 3600;
xlimMax = 1800;
MarkerSize = 15;
figure(3)
tiledlayout(3,2)
% nexttile(1,[3 1])
% plot(Nude_No_CA.time,Nude_No_CA.data(:,2))
% hold on
% plot(Nude_24Hrs_10nmol_OTL38.time,-Nude_24Hrs_10nmol_OTL38.data(:,2))
% plot(Nude_24Hrs_10nmol_VGT.time,Nude_24Hrs_10nmol_VGT.data(:,2))
% hold off
% ylabel('NIR-DiFC Signal (nA)')
% ylim([0 50000])
% yticks(0:5000:50000)
% xlim([0 xlimMax])
% %xticks(0:900:3600)
% %xticklabels({'0','','30','','60'})
% xticks(0:300:3600)
% %xticklabels({'0','','30','','60'})
% xticklabels({'0','','10','','20','','30'})
% xlabel('Time (min)')
% % ylim([ 0 14])
% % yticks(0:1:14)
% % yticklabels({'0','','2','','4','','6','','8','','10','','12','','14'})
% % title('24 Hrs Post VGT-309 Injection')
% axis = gca;
% axis.FontSize =30;
% axis.FontWeight = 'normal';
% axis.FontName = 'Arial';
% axis.LineWidth = 1;
nexttile(2)
plot(prelabeled_leg.time(1:samples_per_sec*1800),prelabeled_leg.data_bs(1:samples_per_sec*1800,2))
hold on
plot(prelabeled_leg.time(prelabeled_leg.fwd_peaks(2).locs), prelabeled_leg.fwd_peaks(2).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(prelabeled_leg.time(prelabeled_leg.rev_peaks(2).locs), prelabeled_leg.rev_peaks(2).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(Nude_24Hrs_10nmol_OTL38.time(1:samples_per_sec*1800),Nude_24Hrs_10nmol_OTL38.data_bs(1:samples_per_sec*1800,2))

plot(Nude_24Hrs_10nmol_OTL38.time(1:samples_per_sec*1800), Nude_24Hrs_10nmol_OTL38.thresh_curve(1:samples_per_sec*1800,2), '-k','LineWidth',1)
hold off
legend({'','2 Probe Matched Cell','','Detection Threshold'},'Location','northoutside','Orientation','horizontal')
%ylabel('NIR-DiFC Signal (nA)')
ylim([-350 1050])
yticks(-350:350:1050)
xlim([0 xlimMax])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30'})
%xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(4)
plot(prelabeled_leg.time(1:samples_per_sec*1800),prelabeled_leg.data_bs(1:samples_per_sec*1800,2))
hold on
plot(prelabeled_leg.time(prelabeled_leg.fwd_peaks(2).locs), prelabeled_leg.fwd_peaks(2).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(prelabeled_leg.time(prelabeled_leg.rev_peaks(2).locs), prelabeled_leg.rev_peaks(2).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(Nude_24Hrs_10nmol_VGT.time(1:samples_per_sec*1800),Nude_24Hrs_10nmol_VGT.data_bs(1:samples_per_sec*1800,2))

plot(Nude_24Hrs_10nmol_VGT.time(1:samples_per_sec*1800), Nude_24Hrs_10nmol_VGT.thresh_curve(1:samples_per_sec*1800,2), '-k','LineWidth',1)
hold off
ylabel('NIR-DiFC Signal (nA)')
ylim([-350 1050])
yticks(-350:350:1050)
xlim([0 xlimMax])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30'})
%xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(6)
plot(prelabeled_leg.time(1:samples_per_sec*1800),prelabeled_leg.data_bs(1:samples_per_sec*1800,2))
hold on
plot(prelabeled_leg.time(prelabeled_leg.fwd_peaks(2).locs), prelabeled_leg.fwd_peaks(2).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(prelabeled_leg.time(prelabeled_leg.rev_peaks(2).locs), prelabeled_leg.rev_peaks(2).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(Nude_No_CA.time(1:samples_per_sec*1800),Nude_No_CA.data_bs(1:samples_per_sec*1800,2))

plot(Nude_No_CA.time(1:samples_per_sec*1800), Nude_No_CA.thresh_curve(1:samples_per_sec*1800,2), '-k','LineWidth',1)
hold off
%ylabel('NIR-DiFC Signal (nA)')
ylim([-350 1050])
yticks(-350:350:1050)
xlim([0 xlimMax])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
%%
%Plotting no injection, 24Hrs post 10nmol OTL38 and VGT-309 background, with peak examples (redo for the new labeling strategy fig).
samples_per_sec = 2000;
start_time = 1800;
end_time = 3600;
xlimMax = 1800;
MarkerSize = 15;
figure(4)
tiledlayout(1,2)
nexttile(1)
plot(prelabeled_leg.time(1:samples_per_sec*1800),prelabeled_leg.data_bs(1:samples_per_sec*1800,2))
hold on
plot(prelabeled_leg.time(prelabeled_leg.fwd_peaks(2).locs), prelabeled_leg.fwd_peaks(2).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(prelabeled_leg.time(prelabeled_leg.rev_peaks(2).locs), prelabeled_leg.rev_peaks(2).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(Nude_No_CA.time(1:samples_per_sec*1800),Nude_No_CA.data_bs(1:samples_per_sec*1800,2))

plot(Nude_No_CA.time(1:samples_per_sec*1800), Nude_No_CA.thresh_curve(1:samples_per_sec*1800,2), '-k','LineWidth',1)
hold off
legend({'','2 Probe Matched Cell','','Detection Threshold'},'Location','northoutside','Orientation','horizontal')
ylabel('NIR-DiFC Signal (nA)')
ylim([-350 1050])
yticks(-350:350:1050)
xlim([0 xlimMax])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(2)
plot(prelabeled_leg.time(1:samples_per_sec*1800),prelabeled_leg.data_bs(1:samples_per_sec*1800,2))
hold on
plot(prelabeled_leg.time(prelabeled_leg.fwd_peaks(2).locs), prelabeled_leg.fwd_peaks(2).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(prelabeled_leg.time(prelabeled_leg.rev_peaks(2).locs), prelabeled_leg.rev_peaks(2).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(Nude_24Hrs_10nmol_OTL38.time(1:samples_per_sec*1800),Nude_24Hrs_10nmol_OTL38.data_bs(1:samples_per_sec*1800,2))

plot(Nude_24Hrs_10nmol_OTL38.time(1:samples_per_sec*1800), Nude_24Hrs_10nmol_OTL38.thresh_curve(1:samples_per_sec*1800,2), '-k','LineWidth',1)
hold off

%ylabel('NIR-DiFC Signal (nA)')
ylim([-350 1050])
yticks(-350:350:1050)
xlim([0 xlimMax])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;

%% Creating a new version of the fluoresence clearance plots
lineWidth = 2;
BarWidth = 0.4;
figure(5)
tiledlayout(3,2)
nexttile(1, [1 2])
plot(VGT_Nude_10nmol_24Hr_1.time, VGT_Nude_10nmol_24Hr_1.data(:,2),'Color','k','LineWidth',lineWidth)
hold on
plot(Nude_NIR_control.time,Nude_NIR_control.data(:,2),'Color','b','LineWidth',lineWidth)
hold off
legend({'CA','Baseline'})
ylim([0 12000])
yticks(0:2000:12000)
yticklabels({'0','','4','','8','','12'})
xlim([0 3600])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30','','40','','50','','60'})
ylabel({'Fluorescence'; 'Background (nA)'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(3)
bar(1,mean(OTL38_Nude_first10min_bg_ratio),BarWidth,'m');
hold on
error_bar = errorbar(1,mean(OTL38_Nude_first10min_bg_ratio),[],std(OTL38_Nude_first10min_bg_ratio),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
bar(2,mean(OTL38_Nude_last10min_bg_ratio),BarWidth,'c');
error_bar = errorbar(2,mean(OTL38_Nude_last10min_bg_ratio),[],std(OTL38_Nude_last10min_bg_ratio),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
hold off
xticks(1:2)
xticklabels({'First 10', 'Last 10'});
xtickangle(0)
ylim([0 21])
yticks(0:1.5:21)
yticklabels({'0','','3','','6','','9','','12','','15','','18','','21'})
ylabel({'Fluorescence Background'; '(Normalized to Baseline)'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('OTL38 - 24Hrs')
nexttile(4)
bar(1,mean(VGT_Nude_first10min_bg_ratio),BarWidth,'m');
hold on
error_bar = errorbar(1,mean(VGT_Nude_first10min_bg_ratio),[],std(VGT_Nude_first10min_bg_ratio),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
bar(2,mean(VGT_Nude_last10min_bg_ratio),BarWidth,'c');
error_bar = errorbar(2,mean(VGT_Nude_last10min_bg_ratio),[],std(VGT_Nude_last10min_bg_ratio),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
hold off
xticks(1:2)
xticklabels({'First 10', 'Last 10'});
xtickangle(0)
ylim([0 21])
yticks(0:1.5:21)
yticklabels({'0','','3','','6','','9','','12','','15','','18','','21'})
%ylabel({'Fluorescence Background'; '(Normalized to Baseline)'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('VGT-309 - 24Hrs')
nexttile(5)
bar(1,PSMA_2_Nude_first10min_bg_ratio,BarWidth,'m');
hold on
bar(2,PSMA_2_Nude_last10min_bg_ratio,BarWidth,'c');
hold off
xticks(1:2)
xticklabels({'First 10', 'Last 10'});
xtickangle(0)
ylim([0 21])
yticks(0:1.5:21)
yticklabels({'0','','3','','6','','9','','12','','15','','18','','21'})
ylabel({'Fluorescence Background'; '(Normalized to Baseline)'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('Cy5-PSMA-02 - 24Hrs')
nexttile(6)
bar(1,PSMA_4_Nude_first10min_bg_ratio,BarWidth,'m');
hold on
bar(2,PSMA_4_Nude_last10min_bg_ratio,BarWidth,'c');
hold off
xticks(1:2)
xticklabels({'First 10', 'Last 10'});
xtickangle(0)
ylim([0 21])
yticks(0:1.5:21)
yticklabels({'0','','3','','6','','9','','12','','15','','18','','21'})
%ylabel({'Fluorescence Background'; '(Normalized to Baseline)'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('Cy5-PSMA-04 - 4Hrs')

%% Creating a new version of the fluoresence clearance plots (landscape for presentations)
lineWidth = 2;
BarWidth = 0.4;
figure(6)
tiledlayout(2,4)
nexttile(1, [2 2])
plot(VGT_Nude_10nmol_24Hr_1.time, VGT_Nude_10nmol_24Hr_1.data(:,2),'Color','k','LineWidth',lineWidth)
hold on
plot(Nude_NIR_control.time,Nude_NIR_control.data(:,2),'Color','b','LineWidth',lineWidth)
hold off
legend({'CA','Baseline'})
ylim([0 12000])
yticks(0:2000:12000)
yticklabels({'0','','4','','8','','12'})
xlim([0 3600])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30','','40','','50','','60'})
ylabel({'Fluorescence'; 'Background (nA)'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(3)
bar(1,mean(OTL38_Nude_first10min_bg_ratio),BarWidth,'m');
hold on
error_bar = errorbar(1,mean(OTL38_Nude_first10min_bg_ratio),[],std(OTL38_Nude_first10min_bg_ratio),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
bar(2,mean(OTL38_Nude_last10min_bg_ratio),BarWidth,'c');
error_bar = errorbar(2,mean(OTL38_Nude_last10min_bg_ratio),[],std(OTL38_Nude_last10min_bg_ratio),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
hold off
xticks(1:2)
xticklabels({'First 10', 'Last 10'});
xtickangle(0)
ylim([0 21])
yticks(0:1.5:21)
yticklabels({'0','','3','','6','','9','','12','','15','','18','','21'})
ylabel({'Fluorescence Background'; '(Normalized to Baseline)'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('OTL38 - 24Hrs')
nexttile(7)
bar(1,mean(VGT_Nude_first10min_bg_ratio),BarWidth,'m');
hold on
error_bar = errorbar(1,mean(VGT_Nude_first10min_bg_ratio),[],std(VGT_Nude_first10min_bg_ratio),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
bar(2,mean(VGT_Nude_last10min_bg_ratio),BarWidth,'c');
error_bar = errorbar(2,mean(VGT_Nude_last10min_bg_ratio),[],std(VGT_Nude_last10min_bg_ratio),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
hold off
xticks(1:2)
xticklabels({'First 10', 'Last 10'});
xtickangle(0)
ylim([0 21])
yticks(0:1.5:21)
yticklabels({'0','','3','','6','','9','','12','','15','','18','','21'})
%ylabel({'Fluorescence Background'; '(Normalized to Baseline)'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('VGT-309 - 24Hrs')
nexttile(4)
bar(1,PSMA_2_Nude_first10min_bg_ratio,BarWidth,'m');
hold on
bar(2,PSMA_2_Nude_last10min_bg_ratio,BarWidth,'c');
hold off
xticks(1:2)
xticklabels({'First 10', 'Last 10'});
xtickangle(0)
ylim([0 21])
yticks(0:1.5:21)
yticklabels({'0','','3','','6','','9','','12','','15','','18','','21'})
ylabel({'Fluorescence Background'; '(Normalized to Baseline)'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('Cy5-PSMA-02 - 24Hrs')
nexttile(8)
bar(1,PSMA_4_Nude_first10min_bg_ratio,BarWidth,'m');
hold on
bar(2,PSMA_4_Nude_last10min_bg_ratio,BarWidth,'c');
hold off
xticks(1:2)
xticklabels({'First 10', 'Last 10'});
xtickangle(0)
ylim([0 21])
yticks(0:1.5:21)
yticklabels({'0','','3','','6','','9','','12','','15','','18','','21'})
%ylabel({'Fluorescence Background'; '(Normalized to Baseline)'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('Cy5-PSMA-04 - 4Hrs')
%% False postive in vivo examples of OTL38 and VGT

samples_per_sec = 2000;
endTime = 2700;
yMax = 500;
yMin = -300;
MarkerSize = 12;
figure(7)
tiledlayout(3,1)
nexttile(1)
plot(Nude_NIR_control.time(1:samples_per_sec*endTime),Nude_NIR_control.data_bs(1:samples_per_sec*endTime,1))
hold on
plot(Nude_NIR_control.time(Nude_NIR_control.unmatched_peaks(1).locs), Nude_NIR_control.unmatched_peaks(1).pks, '*', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)

plot(Nude_NIR_control.time(Nude_NIR_control.fwd_peaks(1).locs), Nude_NIR_control.fwd_peaks(1).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(Nude_NIR_control.time(Nude_NIR_control.rev_peaks(1).locs), Nude_NIR_control.rev_peaks(1).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(Nude_NIR_control.time(1:samples_per_sec*endTime), Nude_NIR_control.thresh_curve(1:samples_per_sec*endTime,1), '-k','LineWidth',1)
hold off
legend({'','Detection Threshold','',''},'Location','north','Orientation','horizontal')
ylabel('NIR-DiFC Signal (nA)')
ylim([yMin yMax])
yticks(yMin:100:yMax)
yticklabels({'-300','','-100','0','100','','300','','500'})
xlim([0 endTime])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30','','40'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(2)
plot(OTL38_Nude_24Hr_Peaks.time(1:samples_per_sec*endTime),OTL38_Nude_24Hr_Peaks.data_bs(1:samples_per_sec*endTime,1))
hold on
plot(OTL38_Nude_24Hr_Peaks.time(OTL38_Nude_24Hr_Peaks.unmatched_peaks(1).locs), OTL38_Nude_24Hr_Peaks.unmatched_peaks(1).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',MarkerSize)

plot(OTL38_Nude_24Hr_Peaks.time(OTL38_Nude_24Hr_Peaks.fwd_peaks(1).locs), OTL38_Nude_24Hr_Peaks.fwd_peaks(1).pks, '^', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(OTL38_Nude_24Hr_Peaks.time(OTL38_Nude_24Hr_Peaks.rev_peaks(1).locs), OTL38_Nude_24Hr_Peaks.rev_peaks(1).pks, '^', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(OTL38_Nude_24Hr_Peaks.time(1:samples_per_sec*endTime), OTL38_Nude_24Hr_Peaks.thresh_curve(1:samples_per_sec*endTime,1), '-k','LineWidth',1)

hold off
legend({'','1 Probe Unmatched Detection','2 Probe Matched Detection','Detection Threshold'},'Location','northoutside','Orientation','horizontal')
ylabel('NIR-DiFC Signal (nA)')
ylim([yMin yMax])
yticks(yMin:100:yMax)
yticklabels({'-300','','-100','0','100','','300','','500'})
xlim([0 endTime])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30','','40'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
startTime = 3000;
nexttile(3)
plot(VGT_Nude_24Hr_Peaks.time(samples_per_sec*startTime:samples_per_sec*(endTime+startTime)),VGT_Nude_24Hr_Peaks.data_bs(samples_per_sec*startTime:samples_per_sec*(endTime+startTime),1))
hold on
plot(VGT_Nude_24Hr_Peaks.time(VGT_Nude_24Hr_Peaks.unmatched_peaks(1).locs), VGT_Nude_24Hr_Peaks.unmatched_peaks(1).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',MarkerSize)

plot(VGT_Nude_24Hr_Peaks.time(VGT_Nude_24Hr_Peaks.fwd_peaks(1).locs), VGT_Nude_24Hr_Peaks.fwd_peaks(1).pks, '^', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(VGT_Nude_24Hr_Peaks.time(VGT_Nude_24Hr_Peaks.rev_peaks(1).locs), VGT_Nude_24Hr_Peaks.rev_peaks(1).pks, '^', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(VGT_Nude_24Hr_Peaks.time(samples_per_sec*startTime:samples_per_sec*(endTime+startTime)), VGT_Nude_24Hr_Peaks.thresh_curve(samples_per_sec*startTime:samples_per_sec*(endTime+startTime),1), '-k','LineWidth',1)

hold off
legend({'','1 Probe Unmatched Detection','2 Probe Matched Detection','Detection Threshold'},'Location','northoutside','Orientation','horizontal')
ylabel('NIR-DiFC Signal (nA)')
ylim([yMin yMax])
yticks(yMin:100:yMax)
yticklabels({'-300','','-100','0','100','','300','','500'})
xlim([startTime startTime+endTime])
xticks(startTime:300:startTime+endTime)
xticklabels({'0','','10','','20','','30','','40'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
%% False postive in vivo examples of OTL38 and VGT and Count rate plot
BarWidth = 0.4;
samples_per_sec = 2000;
endTime = 1800;
yMax = 500;
yMin = -300;
MarkerSize = 12;
figure(8)
tiledlayout(4,2)
nexttile(2)
bar(1,0,BarWidth)
hold on
bar(2,mean(OTL38_10nmol_24Hr_CR),BarWidth);
error_bar = errorbar(2,mean(OTL38_10nmol_24Hr_CR),[],std(OTL38_10nmol_24Hr_CR),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
bar(3,mean(VGT_10nmol_24Hr_CR),BarWidth,'c');
error_bar = errorbar(3,mean(VGT_10nmol_24Hr_CR),[],std(VGT_10nmol_24Hr_CR),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
hold off
xticks(1:3)
xticklabels({'No Injection', 'OTL38','VGT-309'});
xtickangle(0)
ylim([0 8])
yticks(0:1:8)
yticklabels({'0','','2','','4','','6','','8'})
ylabel('Cell Count Rate (/Hr)')
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(3, [1 2])
plot(Nude_NIR_control.time(1:samples_per_sec*endTime),Nude_NIR_control.data_bs(1:samples_per_sec*endTime,1))
hold on
plot(Nude_NIR_control.time(Nude_NIR_control.unmatched_peaks(1).locs), Nude_NIR_control.unmatched_peaks(1).pks, '*', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)

plot(Nude_NIR_control.time(Nude_NIR_control.fwd_peaks(1).locs), Nude_NIR_control.fwd_peaks(1).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(Nude_NIR_control.time(Nude_NIR_control.rev_peaks(1).locs), Nude_NIR_control.rev_peaks(1).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(Nude_NIR_control.time(1:samples_per_sec*endTime), Nude_NIR_control.thresh_curve(1:samples_per_sec*endTime,1), '-k','LineWidth',1)
hold off
legend({'','Detection Threshold','',''},'Location','north','Orientation','horizontal')
ylabel('NIR-DiFC Signal (nA)')
ylim([yMin yMax])
yticks(yMin:100:yMax)
yticklabels({'-300','','-100','0','100','','300','','500'})
xlim([0 endTime])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30','','40'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(5, [1 2])
plot(OTL38_Nude_24Hr_Peaks.time(1:samples_per_sec*endTime),OTL38_Nude_24Hr_Peaks.data_bs(1:samples_per_sec*endTime,1))
hold on
plot(OTL38_Nude_24Hr_Peaks.time(OTL38_Nude_24Hr_Peaks.unmatched_peaks(1).locs), OTL38_Nude_24Hr_Peaks.unmatched_peaks(1).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',MarkerSize)

plot(OTL38_Nude_24Hr_Peaks.time(OTL38_Nude_24Hr_Peaks.fwd_peaks(1).locs), OTL38_Nude_24Hr_Peaks.fwd_peaks(1).pks, '^', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(OTL38_Nude_24Hr_Peaks.time(OTL38_Nude_24Hr_Peaks.rev_peaks(1).locs), OTL38_Nude_24Hr_Peaks.rev_peaks(1).pks, '^', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(OTL38_Nude_24Hr_Peaks.time(1:samples_per_sec*endTime), OTL38_Nude_24Hr_Peaks.thresh_curve(1:samples_per_sec*endTime,1), '-k','LineWidth',1)

hold off
legend({'','1 Probe Unmatched Detection','2 Probe Matched Detection','Detection Threshold'},'Location','northoutside','Orientation','horizontal')
ylabel('NIR-DiFC Signal (nA)')
ylim([yMin yMax])
yticks(yMin:100:yMax)
yticklabels({'-300','','-100','0','100','','300','','500'})
xlim([0 endTime])
xticks(0:300:3600)
xticklabels({'0','','10','','20','','30','','40'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
startTime = 2400;
nexttile(7, [1 2])
plot(VGT_Nude_24Hr_Peaks.time(samples_per_sec*startTime:samples_per_sec*(endTime+startTime)),VGT_Nude_24Hr_Peaks.data_bs(samples_per_sec*startTime:samples_per_sec*(endTime+startTime),1))
hold on
plot(VGT_Nude_24Hr_Peaks.time(VGT_Nude_24Hr_Peaks.unmatched_peaks(1).locs), VGT_Nude_24Hr_Peaks.unmatched_peaks(1).pks, 'o', 'Color', '#0072B2','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',MarkerSize)

plot(VGT_Nude_24Hr_Peaks.time(VGT_Nude_24Hr_Peaks.fwd_peaks(1).locs), VGT_Nude_24Hr_Peaks.fwd_peaks(1).pks, '^', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(VGT_Nude_24Hr_Peaks.time(VGT_Nude_24Hr_Peaks.rev_peaks(1).locs), VGT_Nude_24Hr_Peaks.rev_peaks(1).pks, '^', 'Color', '#0072B2','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',MarkerSize)
plot(VGT_Nude_24Hr_Peaks.time(samples_per_sec*startTime:samples_per_sec*(endTime+startTime)), VGT_Nude_24Hr_Peaks.thresh_curve(samples_per_sec*startTime:samples_per_sec*(endTime+startTime),1), '-k','LineWidth',1)

hold off
legend({'','1 Probe Unmatched Detection','2 Probe Matched Detection','Detection Threshold'},'Location','northoutside','Orientation','horizontal')
ylabel('NIR-DiFC Signal (nA)')
ylim([yMin yMax])
yticks(yMin:100:yMax)
yticklabels({'-300','','-100','0','100','','300','','500'})
xlim([startTime startTime+endTime])
xticks(startTime:300:startTime+endTime)
xticklabels({'0','','10','','20','','30','','40'})
xlabel('Time (min)')
axis = gca;
axis.FontSize =30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;

%% Contrast Agent BALB/C Cardiac Puncture Phantom Results (without ABY-029 and condensed)
width = 0.4;
figure(9);
tiledlayout(2,2);
nexttile(1)
bar(1,cells_per_ml_blood_No_CA_BALB_C_NIR,width);
hold on
bar(2,cells_per_ml_blood_3Hrs_OTL38_BALB_C,width);

bar(3,cells_per_ml_blood_24Hrs_OTL38_BALB_C,width);
hold off
xticks(1:3)
xticklabels({'Blood Only','3Hrs', '24Hrs'});
xtickangle(0)
ylim([0 2000])
yticks(0:250:2000)
yticklabels({'0','','500','','1000','','1500','','2000'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('OTL38')
nexttile(2)
bar(1,cells_per_ml_blood_No_CA_BALB_C_NIR,width);
hold on
bar(2,cells_per_ml_blood_48Hrs_VGT_BALB_C,width);
hold off
xticks(1:2)
xticklabels({'Blood Only','48Hrs'});
xtickangle(0)
ylim([0 2000])
yticks(0:250:2000)
yticklabels({'0','','500','','1000','','1500','','2000'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('VGT-309')
nexttile(3)
bar(1,cells_per_ml_blood_No_CA_BALB_C_NIR,width);
hold on
bar(2,cells_per_ml_blood_4Hrs_Cy5_PSMA_02_BALB_C,width);
bar(3,cells_per_ml_blood_24Hrs_Cy5_PSMA_02_BALB_C,width);
hold off
xticks(1:3)
xticklabels({'Blood Only','4Hrs', '24Hrs'});
xtickangle(0)
ylim([0 2000])
yticks(0:250:2000)
yticklabels({'0','','500','','1000','','1500','','2000'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('Red-PSMA-02')
nexttile(4)
bar(1,cells_per_ml_blood_No_CA_BALB_C_NIR,width);
hold on
bar(2,cells_per_ml_blood_4Hrs_Cy5_PSMA_04_BALB_C,width);
bar(3,cells_per_ml_blood_24Hrs_Cy5_PSMA_04_BALB_C,width);
hold off
xticks(1:3)
xticklabels({'Blood Only','4Hrs', '24Hrs'});
xtickangle(0)
ylim([0 2000])
yticks(0:250:2000)
yticklabels({'0','','500','','1000','','1500','','2000'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('Red-PSMA-04')
fig = gcf;
tl= fig.Children;
ylabel(tl, 'Non-Cancerous Cell DiFC Detections (/mL Blood)', 'FontSize',40,'FontWeight','normal','FontName','Arial');
%xlabel(tl,'Peak Amplitude (mV)', 'FontSize', 40,'FontWeight','bold','FontName','Arial');

%%
%New background figure ( shoter and longer timepoints for OTL38, VGT-309, PSMA-02 & 04)
yMax = 12;
width = 0.4;
figure(10);
tiledlayout(2,2);
nexttile(1)
bar(1,0,width);
hold on
bar(2,mean(OTL38_24hr_background_ratio_2),width);
error_bar = errorbar(2,mean(OTL38_24hr_background_ratio_2),[],std(OTL38_24hr_background_ratio_2),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
hold off
xticks(1:2)
xticklabels({'3Hrs', '24Hrs'});
xtickangle(0)
ylim([0 yMax])
yticks(0:1:yMax)
yticklabels({'0','','2','','4','','6','','8','','10','','12'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(2)
bar(1,mean(VGT_2hr_background_ratio_1),width);
hold on
error_bar = errorbar(1,mean(VGT_2hr_background_ratio_1),[],std(VGT_2hr_background_ratio_1),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];

bar(2,mean(VGT_24hr_background_ratio_1),width);
error_bar = errorbar(2,mean(VGT_24hr_background_ratio_1),[],std(VGT_24hr_background_ratio_1),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
hold off
xticks(1:2)
xticklabels({'2Hrs', '24Hrs'});
xtickangle(0)
ylim([0 yMax])
yticks(0:1:yMax)
yticklabels({'0','','2','','4','','6','','8','','10','','12'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(3)
bar(1,0,width);
hold on
bar(2,mean(Cy5_PSMA_2_Nude_10nmol_24Hr_ratio_1),width);
error_bar = errorbar(2,mean(Cy5_PSMA_2_Nude_10nmol_24Hr_ratio_1),[],std(Cy5_PSMA_2_Nude_10nmol_24Hr_ratio_1),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
hold off
xticks(1:2)
xticklabels({'4Hrs', '24Hrs'});
xtickangle(0)
ylim([0 yMax])
yticks(0:1:yMax)
yticklabels({'0','','2','','4','','6','','8','','10','','12'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
nexttile(4)
bar(1,mean(Cy5_PSMA_4_Nude_10nmol_4Hr_ratio_1),width);
hold on
error_bar = errorbar(1,mean(Cy5_PSMA_4_Nude_10nmol_4Hr_ratio_1),[],std(Cy5_PSMA_4_Nude_10nmol_4Hr_ratio_1),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];

bar(2,mean(Cy5_PSMA_2_Nude_10nmol_24Hr_ratio_1),width);
error_bar = errorbar(2,mean(Cy5_PSMA_2_Nude_10nmol_24Hr_ratio_1),[],std(Cy5_PSMA_2_Nude_10nmol_24Hr_ratio_1),'LineWidth',1,'CapSize',15);
error_bar.Color = [0 0 0];
hold off
xticks(1:2)
xticklabels({'4Hrs', '24Hrs'});
xtickangle(0)
ylim([0 yMax])
yticks(0:1:yMax)
yticklabels({'0','','2','','4','','6','','8','','10','','12'})

axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
fig = gcf;
tl= fig.Children;
ylabel(tl, 'Fluorescence Background (Normalized to Baseline)', 'FontSize',40,'FontWeight','normal','FontName','Arial');

%% Red PSMA-02 04 CD1 Repeats
%Plotting CD-1 PSMA-02 04 and PBS


% swarmchart(x1_2Hrs , LLC_VGT_CR_2Hrs,100,'filled','k','XJitterWidth',0.2)
% txt = num2str(cells_per_ml_blood_4Hrs_PBS_CD_1_1);
% text(0.8,20,txt,'FontSize',fSize)
x1= ones(1,3);
x2= ones(1,3)*2;
x3= ones(1,3)*3;
x4= ones(1,3)*4;
width = 0.4;
fSize  = 30;
height = 130;
figure(11);
tiledlayout(1,2);
%4Hrs
nexttile(1)
bar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_mean,width);% Male CD-1 4Hrs post PBS injection in phantom 
hold on
error_bar = errorbar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_mean,[],cells_per_ml_blood_24Hrs_PBS_CD_1_mean,'k','LineWidth',1,'CapSize',15);
%swarmchart(x1 , cells_per_ml_blood_4Hrs_PBS_CD_1,100,'filled','k','XJitterWidth',0.2)

bar(2,cells_per_ml_blood_4Hrs_Cy5_PSMA_02_CD_1_mean,width);% Male CD-1 4Hrs Post 10nmol Cy5 PSMA-02 in phantom 
error_bar = errorbar(2,cells_per_ml_blood_4Hrs_Cy5_PSMA_02_CD_1_mean,[],cells_per_ml_blood_4Hrs_Cy5_PSMA_02_CD_1_std,'k','LineWidth',1,'CapSize',15);
%swarmchart(x3 , cells_per_ml_blood_4Hrs_Cy5_PSMA_02_CD_1,100,'filled','k','XJitterWidth',0.2)

bar(3,cells_per_ml_blood_4Hrs_Cy5_PSMA_04_CD_1_mean,width);% Male CD-1 4Hrs Post 10nmol Cy5 PSMA-04 in phantom 
error_bar = errorbar(3,cells_per_ml_blood_4Hrs_Cy5_PSMA_04_CD_1_mean,[],cells_per_ml_blood_4Hrs_Cy5_PSMA_04_CD_1_std,'k','LineWidth',1,'CapSize',15);
%swarmchart(x2 , cells_per_ml_blood_4Hrs_Cy5_PSMA_04_CD_1,100,'filled','k','XJitterWidth',0.2)


hold off
xticks(1:3)
xticklabels({'PBS','PSMA-02','PSMA-04'});
ylabel('Male CD-1 Immune Cell DiFC Detections (/mL Blood)','FontSize',40,'FontWeight','normal','FontName','Arial');
title('4Hrs Post Injection','FontSize',40,'FontWeight','normal','FontName','Arial');
ylim([0 250])
yticks(0:25:250)
yticklabels({'0','','50','','100','','150','','200','','250'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
%24Hrs
nexttile(2)

% bar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_high_mean,width);% Male CD-1 24Hrs post PBS injection in phantom 
% hold on
% error_bar = errorbar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_high_mean,[],cells_per_ml_blood_24Hrs_PBS_CD_1_high_std,'k','LineWidth',1,'CapSize',15);
% %swarmchart(x1 , cells_per_ml_blood_24Hrs_PBS_CD_1_high,100,'filled','k','XJitterWidth',0.2)

bar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_mean,width);% Male CD-1 24Hrs post PBS injection in phantom 
hold on
error_bar = errorbar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_mean,[],cells_per_ml_blood_24Hrs_PBS_CD_1_std,'k','LineWidth',1,'CapSize',15);
%swarmchart(x2, cells_per_ml_blood_24Hrs_PBS_CD_1,100,'filled','k','XJitterWidth',0.2)

bar(2,cells_per_ml_blood_24Hrs_Cy5_PSMA_02_CD_1_mean,width);% Male CD-1 24Hrs Post 10nmol Cy5 PSMA-02 in phantom 
error_bar = errorbar(2,cells_per_ml_blood_24Hrs_Cy5_PSMA_02_CD_1_mean,[],cells_per_ml_blood_24Hrs_Cy5_PSMA_02_CD_1_std,'k','LineWidth',1,'CapSize',15);
%swarmchart(x4 , cells_per_ml_blood_24Hrs_Cy5_PSMA_02_CD_1,100,'filled','k','XJitterWidth',0.2)


bar(3,cells_per_ml_blood_24Hrs_Cy5_PSMA_04_CD_1_mean,width);% Male CD-1 24Hrs Post 10nmol Cy5 PSMA-04 in phantom 
error_bar = errorbar(3,cells_per_ml_blood_24Hrs_Cy5_PSMA_04_CD_1_mean,[],cells_per_ml_blood_24Hrs_Cy5_PSMA_04_CD_1_std,'k','LineWidth',1,'CapSize',15);
%swarmchart(x3 , cells_per_ml_blood_24Hrs_Cy5_PSMA_04_CD_1,100,'filled','k','XJitterWidth',0.2)

hold off
xticks(1:3)
xticklabels({'PBS','PSMA-02','PSMA-04'});
%ylabel('','FontSize',40,'FontWeight','bold','FontName','Arial');
title('24Hrs Post Injection','FontSize',40,'FontWeight','normal','FontName','Arial');
ylim([0 250])
yticks(0:25:250)
yticklabels({'0','','50','','100','','150','','200','','250'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;

%% Combining BALBC phatom with the PSMA CD-1 repeats

width = 0.4;
figure(12);
tiledlayout(3,2);
nexttile(1)
bar(1,cells_per_ml_blood_No_CA_BALB_C_NIR,width);
hold on
bar(2,cells_per_ml_blood_3Hrs_OTL38_BALB_C,width);

bar(3,cells_per_ml_blood_24Hrs_OTL38_BALB_C,width);
hold off
xticks(1:3)
xticklabels({'Blood Only','3Hrs', '24Hrs'});
xtickangle(0)
ylim([0 2000])
yticks(0:250:2000)
yticklabels({'0','','500','','1000','','1500','','2000'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('OTL38')
nexttile(2)
bar(1,cells_per_ml_blood_No_CA_BALB_C_NIR,width);
hold on
bar(2,cells_per_ml_blood_48Hrs_VGT_BALB_C,width);
hold off
xticks(1:2)
xticklabels({'Blood Only','48Hrs'});
xtickangle(0)
ylim([0 2000])
yticks(0:250:2000)
yticklabels({'0','','500','','1000','','1500','','2000'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('VGT-309')
nexttile(3)
bar(1,cells_per_ml_blood_No_CA_BALB_C_NIR,width);
hold on
bar(2,cells_per_ml_blood_4Hrs_Cy5_PSMA_02_BALB_C,width);
bar(3,cells_per_ml_blood_24Hrs_Cy5_PSMA_02_BALB_C,width);
hold off
xticks(1:3)
xticklabels({'Blood Only','4Hrs', '24Hrs'});
xtickangle(0)
ylim([0 2000])
yticks(0:250:2000)
yticklabels({'0','','500','','1000','','1500','','2000'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('Red-PSMA-02')
nexttile(4)
bar(1,cells_per_ml_blood_No_CA_BALB_C_NIR,width);
hold on
bar(2,cells_per_ml_blood_4Hrs_Cy5_PSMA_04_BALB_C,width);
bar(3,cells_per_ml_blood_24Hrs_Cy5_PSMA_04_BALB_C,width);
hold off
xticks(1:3)
xticklabels({'Blood Only','4Hrs', '24Hrs'});
xtickangle(0)
ylim([0 2000])
yticks(0:250:2000)
yticklabels({'0','','500','','1000','','1500','','2000'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
title('Red-PSMA-04')

nexttile(5)
bar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_mean,width);% Male CD-1 4Hrs post PBS injection in phantom 
hold on
error_bar = errorbar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_mean,[],cells_per_ml_blood_24Hrs_PBS_CD_1_mean,'k','LineWidth',1,'CapSize',15);

bar(2,cells_per_ml_blood_4Hrs_Cy5_PSMA_02_CD_1_mean,width);% Male CD-1 4Hrs Post 10nmol Cy5 PSMA-02 in phantom 
error_bar = errorbar(2,cells_per_ml_blood_4Hrs_Cy5_PSMA_02_CD_1_mean,[],cells_per_ml_blood_4Hrs_Cy5_PSMA_02_CD_1_std,'k','LineWidth',1,'CapSize',15);

bar(3,cells_per_ml_blood_4Hrs_Cy5_PSMA_04_CD_1_mean,width);% Male CD-1 4Hrs Post 10nmol Cy5 PSMA-04 in phantom 
error_bar = errorbar(3,cells_per_ml_blood_4Hrs_Cy5_PSMA_04_CD_1_mean,[],cells_per_ml_blood_4Hrs_Cy5_PSMA_04_CD_1_std,'k','LineWidth',1,'CapSize',15);


hold off
xticks(1:3)
xticklabels({'PBS','PSMA-02','PSMA-04'});
%ylabel('Male CD-1 Immune Cell DiFC Detections (/mL Blood)','FontSize',40,'FontWeight','normal','FontName','Arial');
title('4Hrs Post Injection','FontSize',40,'FontWeight','normal','FontName','Arial');
ylim([0 250])
yticks(0:25:250)
yticklabels({'0','','50','','100','','150','','200','','250'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;
%24Hrs
nexttile(6)

% bar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_high_mean,width);% Male CD-1 24Hrs post PBS injection in phantom 
% hold on
% error_bar = errorbar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_high_mean,[],cells_per_ml_blood_24Hrs_PBS_CD_1_high_std,'k','LineWidth',1,'CapSize',15);
% %swarmchart(x1 , cells_per_ml_blood_24Hrs_PBS_CD_1_high,100,'filled','k','XJitterWidth',0.2)

bar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_mean,width);% Male CD-1 24Hrs post PBS injection in phantom 
hold on
error_bar = errorbar(1,cells_per_ml_blood_24Hrs_PBS_CD_1_mean,[],cells_per_ml_blood_24Hrs_PBS_CD_1_std,'k','LineWidth',1,'CapSize',15);
%swarmchart(x2, cells_per_ml_blood_24Hrs_PBS_CD_1,100,'filled','k','XJitterWidth',0.2)

bar(2,cells_per_ml_blood_24Hrs_Cy5_PSMA_02_CD_1_mean,width);% Male CD-1 24Hrs Post 10nmol Cy5 PSMA-02 in phantom 
error_bar = errorbar(2,cells_per_ml_blood_24Hrs_Cy5_PSMA_02_CD_1_mean,[],cells_per_ml_blood_24Hrs_Cy5_PSMA_02_CD_1_std,'k','LineWidth',1,'CapSize',15);
%swarmchart(x4 , cells_per_ml_blood_24Hrs_Cy5_PSMA_02_CD_1,100,'filled','k','XJitterWidth',0.2)


bar(3,cells_per_ml_blood_24Hrs_Cy5_PSMA_04_CD_1_mean,width);% Male CD-1 24Hrs Post 10nmol Cy5 PSMA-04 in phantom 
error_bar = errorbar(3,cells_per_ml_blood_24Hrs_Cy5_PSMA_04_CD_1_mean,[],cells_per_ml_blood_24Hrs_Cy5_PSMA_04_CD_1_std,'k','LineWidth',1,'CapSize',15);
%swarmchart(x3 , cells_per_ml_blood_24Hrs_Cy5_PSMA_04_CD_1,100,'filled','k','XJitterWidth',0.2)

hold off
xticks(1:3)
xticklabels({'PBS','PSMA-02','PSMA-04'});
%ylabel('','FontSize',40,'FontWeight','bold','FontName','Arial');
title('24Hrs Post Injection','FontSize',40,'FontWeight','normal','FontName','Arial');
ylim([0 250])
yticks(0:25:250)
yticklabels({'0','','50','','100','','150','','200','','250'})
axis = gca;
axis.FontSize = 30;
axis.FontWeight = 'normal';
axis.FontName = 'Arial';
axis.LineWidth = 1;

fig = gcf;
tl= fig.Children;
ylabel(tl, 'Non-Cancerous Cell DiFC Detections (/mL Blood)', 'FontSize',40,'FontWeight','normal','FontName','Arial');
%xlabel(tl,'Peak Amplitude (mV)', 'FontSize', 40,'FontWeight','bold','FontName','Arial');