%% plot the figures for the truncation paper
%% housekeeping
clear all; close all; clc
%% load the data
load("simulation_results_no_fade.mat")
% load("simulation_results_fade_in.mat")
fs=48000;
time=0:1/fs:length(r_cov_wcr)/fs -1/fs;
SNR = 5:5:30;
%% colors for plotting
numPlots = 6;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

set(groot,'defaultAxesTickLabelInterpreter','latex'); 

%% plot unimodal regression
norm_factor = max(r_cov_cor,[],1);
ff = figure(2); clf; hold on
n = 6

sp1 = subplot(1,2,1)
plot(time, 10^-16+abs(r_c)/norm_factor(:,n), 'k', 'LineWidth',2), hold on
plot(time, 10^-16+abs(r_cov_unc(:,n))/norm_factor(:,n), 'r', 'LineWidth',2)
set(gca, 'yscale', 'log', 'YTick', 10.^[-6, -4, -2, 0])
xlim([0.2 2])
ylim([10^-5, 10])



ylabel('Covariance', 'Interpreter','latex', 'FontSize',12)
xlabel('Time (s)', 'Interpreter','latex', 'FontSize',12)
sp1.Position(1)=0.08
sp1.Position(3)=0.41
sp1.Position(4)=0.75
sp1.Position(2)=0.17





% plot the rectangle showing the range of the inline figure
rectangle('Position', [0.65 10^-1 0.08 0.1 ],  'LineWidth',1)

% plot the dashed lines from rectangle in IR to the inlay axes
plot([0.65 1.05], 10^-1.*[1 0.03], 'k--',  'LineWidth',1)
plot([0.73 1.9], 10^-1.*[1 0.06], 'k--',  'LineWidth',1)
plot([0.65 1.05], 10^-1.*[2 7], 'k--',  'LineWidth',1)
plot([0.73 1.9], 10^-1.*[2 7], 'k--',  'LineWidth',1)


lgd1 = legend({'Raw covariance', 'With unimodal regression'}, "Location","north", 'Interpreter','latex', 'numcolumns', 2)
lgd1.Position(2) = 0.83
lgd1.Position(1) = 0.09

ax = axes('position',[.27 .48 .2 .30]);
plot(time, r_c, 'k', 'LineWidth',2), hold on
plot(time, r_cov_unc(:,n), 'r', 'LineWidth',2)
xlim([0.7 0.705])
set(gca, 'YTickLabel', [], 'XTickLabel', [])
% ylim([-0.002 0.006])


sp2 = subplot(1,2,2)
plot(time, r_cov_cor(:,n)/norm_factor(:,n), 'color', cMap(:, 1), 'LineWidth',2); hold on
plot(time, r_cov_wcr(:,n)/norm_factor(:,n), 'color', cMap(:, 2), 'LineWidth',2)
plot(time, r_cov_unc(:,n)/norm_factor(:,n), 'color', cMap(:, 3), 'LineWidth',2)
plot(time, r_cov_ant(:,n)/norm_factor(:,n), 'color', cMap(:, 4), 'LineWidth',2)

xline(op_ground_truth, 'r-.', LineWidth=2)
xline(tp_ground_truth(n), 'r', LineWidth=2)

yline(v_n_cor(n)/norm_factor(:,n), 'k', LineWidth=2)


set(gca, 'yscale', 'log', 'YTick', 10.^[-6, -4, -2, 0])
xlim([0.2 2])
ylim([10^-5, 10])

xlabel('Time (s)', 'Interpreter','latex', 'FontSize',12)
sp2.Position(1)=0.57
sp2.Position(3)=0.41
sp2.Position(4)=0.75
sp2.Position(2)=0.17

ff.Position(4)=250
ff.Position(3)=800

%% print figure
set(ff,'Units','Inches');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[ff.Position(3), ff.Position(4)])
print(ff,'Covariance_unimodal','-dpdf','-r0')


%% table 1 results
% truncation
100*(tp_unc-tp_ground_truth)./tp_ground_truth
100*(tp_cor-tp_ground_truth)./tp_ground_truth
100*(tp_wcr-tp_ground_truth)./tp_ground_truth
100*(tp_ant-tp_ground_truth)./tp_ground_truth

% onset
100*(op_unc-op_ground_truth)./op_ground_truth
100*(op_cor-op_ground_truth)./op_ground_truth
100*(op_wcr-op_ground_truth)./op_ground_truth
100*(op_ant-op_ground_truth)./op_ground_truth

%% table 2 results
% noise variance difference
100*(v_noise_len-v_n_cor)./v_n_cor

% truncation time difference
mean(100*(tp_unc_len-tp_ground_truth)./tp_ground_truth, 1)
std(100*(tp_unc_len-tp_ground_truth)./tp_ground_truth,0, 1)
% onset time difference
mean(100*(op_unc_len-op_ground_truth)./op_ground_truth, 1)
std(100*(op_unc_len-op_ground_truth)./op_ground_truth, 0, 1)
