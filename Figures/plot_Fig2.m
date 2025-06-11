%% plot the figures for the truncation paper
%% housekeeping
clear all; close all; clc
%% load the data
load("Arni_results.mat")
t_n(177,:) = 0;
fs=44100;
time=0:1/fs:length(filtered_signal)/fs -1/fs;
%% colors for plotting
numPlots = 7;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

set(groot,'defaultAxesTickLabelInterpreter','latex'); 

freq = [250, 500, 1000, 2000, 4000, 8000];
%% 
num_rir = 599;
f = figure(1); clf; hold on
sp1 = subplot(2,3,1); hold on
for n = 1:length(freq)
scatter(rand(1, num_rir)/10+n-0.15, t_n(1:num_rir, n),100,'k', 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.15, 'HandleVisibility','off')
scatter(rand(1, num_rir)/10+n+0.15, inter_time(1:num_rir, n),100, 'r', 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.15, 'HandleVisibility','off')

scatter(n-0.10, median(t_n(1:num_rir, n)),35, 'k', 'filled', 'o', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',1, 'HandleVisibility','off', 'LineWidth',2)
scatter(n+0.20, median(inter_time(1:num_rir, n)), 35, 'r', 'filled', 'o', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',1, 'HandleVisibility','off', 'LineWidth',2)

end 

xlim([0.5 6.5])
ylim([0.2 1])
set(gca, 'Xtick', 1:6,  'XTickLabel', {'250', '500', '1k', '2k', '4k', '8k'}, 'fontsize', 11)



scatter(150,150, 'k', 'filled', 'o', 'MarkerEdgeColor','none')
scatter(150,150, 'r', 'filled', 'o', 'MarkerEdgeColor','none')

lgd = legend('Proposed', 'Lundeby', 'location', 'northeast', 'interpreter', 'latex', 'numcolumns',2);
% lgd.Position(2) = 0.83; 
% lgd.Position(1) = 0.14;


box on
xlabel({'Frequency (Hz)', '(a)'}, 'Interpreter','latex', 'fontsize', 11)
ylabel('$\hat{t}_\textrm{T}$ (s)', 'Interpreter','latex', 'fontsize', 11)


nir = 250
nf = 2
sp2 = subplot(2,3,4); hold on
plot(time, db(filtered_signal(:, nir, nf)./max(abs(filtered_signal(:, nir, nf)))), 'Color', [1 1 1]*0.75, 'HandleVisibility','off')
xline(t_n( nir, nf), 'k', 'LineWidth',2)%, 'label', 'Proposed', 'Interpreter','latex')
xline(inter_time(nir, nf), 'r', 'LineWidth',2)%, 'label', 'Lundeby', 'Interpreter','latex')

xlim([0 0.8])
ylim([-90 0])

xlabel({'Time (s)', '(c)'}, 'Interpreter','latex', 'fontsize', 11)
ylabel('Signal value (dB)', 'Interpreter','latex', 'fontsize', 11)
set(gca, 'fontsize', 11)

box on

lgd = legend('Proposed', 'Lundeby', 'location', 'southeast', 'interpreter', 'latex');

%%
load("MRTD_results_2.mat")
num_rir_mrtd = 50;
fs_=48000
time_2 = 0:1/fs_:length(filtered_signal_mrtd)/fs_ -1/fs_;
inter_time_mrtd(inter_time_mrtd>2)=0;
%% trucation point proposed vs Lubndeby

sp3 = subplot(2,3,2);cla(sp3); hold on
for n = 1:length(freq)+1

scatter(rand(1, num_rir_mrtd)/10+n-0.12, t_n_is(:,  n),80, 'k', 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.15, 'HandleVisibility','off')
scatter(rand(1, num_rir_mrtd)/10+n+0.12, inter_time_mrtd(:, n),80, 'r', 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.15, 'HandleVisibility','off')


end 

xlim([0.5 7.5])
ylim([-0.1 2.5])
set(gca, 'Xtick', 1:7,  'XTickLabel', {'250', '500', '1k', '2k', '4k', '8k', '16k'}, 'fontsize', 11)

scatter(150,150, 'k', 'filled', 'o', 'MarkerEdgeColor','none')
scatter(150,150, 'r', 'filled', 'o', 'MarkerEdgeColor','none')


box on
xlabel({'Frequency (Hz)', '(b)'}, 'Interpreter','latex', 'fontsize', 11)
ylabel('$\hat{t}_\textrm{T}$ (s)', 'Interpreter','latex', 'fontsize', 11)

lgd = legend('Proposed', 'Lundeby', 'location', 'northeast', 'interpreter', 'latex', 'numcolumns', 2);

%%
nf = 1
nir = 10
sp4 = subplot(2,3,5);cla(sp4); hold on
plot(time_2, db(filtered_signal_mrtd(:, 1, nir, nf)./max(abs(filtered_signal_mrtd(:, 1, nir, nf)))), 'Color', [1 1 1]*0.75, 'HandleVisibility','off')
xline(t_n_is(nir, nf), 'k', 'LineWidth',2)%, 'label', 'Proposed', 'Interpreter','latex')
xline(inter_time_mrtd(nir, nf), 'r', 'LineWidth',2)%, 'label', 'Lundeby', 'Interpreter','latex')

xlim([0 2])
ylim([-60 0])

xlabel({'Time (s)', '(d)'}, 'Interpreter','latex', 'fontsize', 12)
% ylabel('Signal value (dB)', 'Interpreter','latex')
set(gca, 'fontsize', 11)
box on

lgd = legend('Proposed', 'Lundeby', 'location', 'southeast', 'interpreter', 'latex');

%%
nf = 7
nir = 11
sp5 = subplot(2,3,6); cla(sp5); hold on
plot(time_2, db(filtered_signal_mrtd(:, 1, nir, nf)./max(abs(filtered_signal_mrtd(:, 1, nir, nf)))), 'Color', [1 1 1]*0.75, 'HandleVisibility','off')
xline(t_n_is(nir, nf), 'k', 'LineWidth',2)%, 'label', 'Proposed', 'Interpreter','latex')
xline(inter_time_mrtd(nir, nf), 'r', 'LineWidth',2)%, 'label', 'Lundeby', 'Interpreter','latex')

xlim([0 1.5])
ylim([-60 0])
set(gca, 'fontsize', 11)
xlabel({'Time (s)', '(e)'}, 'Interpreter','latex', 'fontsize', 12)
% ylabel('Signal value (dB)', 'Interpreter','latex')

box on

lgd = legend('Proposed', 'Lundeby', 'location', 'southeast', 'interpreter', 'latex');

%%
sp1.Position(3) = 0.355
sp1.Position(4) = 0.35
sp1.Position(2) = 0.6

sp3.Position(1) = 0.55%sp4.Position(1)
sp3.Position(3) = 0.355
sp3.Position(4) = 0.35
sp3.Position(2) = 0.6


%%
f.Position(3) = 1200
f.Position(4) = 500
%% print figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'Truncation_evaluation','-dpdf','-r0')

%% additional colors
cy = [255, 190, 11]./255;
co = [106, 153, 78]./255

%% onset point proposed vs several methods

ff = figure(2); clf, hold on

spp1 = subplot(1,3,1); hold on
for n = 1:length(freq)+1
scatter(rand(1, num_rir_mrtd)/10+n-0.3, t_s_is(:,n), 'k', 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.25, 'HandleVisibility','off')
scatter(rand(1, num_rir_mrtd)/10+n-0.15, M(:,  n), 20, cMap(:, 1).', 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.25, 'HandleVisibility','off')
scatter(rand(1, num_rir_mrtd)/10+n, M_5(:,  n),20,  cMap(:, 6).', 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.25, 'HandleVisibility','off')
scatter(rand(1, num_rir_mrtd)/10+n+0.15, D_E(:,  n),20,  cy, 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.25, 'HandleVisibility','off')
scatter(rand(1, num_rir_mrtd)/10+n+0.3, E(:,  n),20,  co, 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.25, 'HandleVisibility','off')

scatter(150,150, 20, 'k', 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',1)
scatter(150,150,  20, cMap(:, 1).', 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',1)
scatter(150,150, 20,  cMap(:, 6).', 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',1)
scatter(150,150, 20,  cy, 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',1)
scatter(150,150, 20,  co, 'filled', 'o', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',1)

end 
box on



xlim([0.5 6.5])
ylim([0.0 0.4])
set(gca, 'Xtick', 1:6,  'XTickLabel', {'250', '500', '1k', '2k', '4k', '8k'}, 'fontsize', 11)

xlabel({'Frequency (Hz)', '(a)'}, 'Interpreter','latex', 'fontsize', 11)
ylabel('$\hat{t}_\textrm{O}$ (s)', 'Interpreter','latex', 'fontsize', 11)

lgd = legend('Proposed', '$M$', '$M_5$', '$D_E$','$E$','location', 'northeast', 'interpreter', 'latex', 'numcolumns', 2, 'fontsize', 8);

%% 
spp2 = subplot(1,3,2);cla(spp2);hold on

nf = 7
nir = 36

plot(time_2, db(filtered_signal_mrtd(:, 1, nir, nf)./max(abs(filtered_signal_mrtd(:, 1, nir, nf)))), 'Color', [1 1 1]*0.75, 'HandleVisibility','off')
plot(t_s_is( nir, nf), -10, 'ko-', 'LineWidth',1, 'MarkerFaceColor','k', 'MarkerSize',4)
xline(t_s_is( nir, nf), 'k', 'HandleVisibility','off')
plot(M( nir, nf), -10,'o', 'color', cMap(:, 1).', 'LineWidth',2, 'MarkerFaceColor',cMap(:, 1).', 'MarkerSize',3)
plot(M_5( nir, nf), -20,'o', 'color', cMap(:, 6).', 'LineWidth',2, 'MarkerFaceColor',cMap(:, 6).', 'MarkerSize',3)
plot(D_E( nir, nf), -20,'o', 'color', cy.', 'LineWidth',2, 'MarkerFaceColor',cy.', 'MarkerSize',3)
plot(E( nir, nf), -15,'o', 'color', co.', 'LineWidth',2, 'MarkerFaceColor',co.', 'MarkerSize',3)

xlim([-0.01 0.7])
ylim([-60 2])

xlabel({'Time (s)', '(b)'}, 'Interpreter','latex', 'fontsize', 11)
ylabel('Signal value (dB)', 'Interpreter','latex', 'fontsize', 11)
set(gca, 'fontsize', 11)

box on

lgd = legend('Proposed', '$M$', '$M_5$', '$D_E$','$E$','location', 'southeast', 'interpreter', 'latex', 'numcolumns', 2, 'fontsize', 8);
%% 
spp3 = subplot(1,3,3);cla(spp3);hold on

nf = 2
nir = 41

plot(time_2, db(filtered_signal_mrtd(:, 1, nir, nf)./max(abs(filtered_signal_mrtd(:, 1, nir, nf)))), 'Color', [1 1 1]*0.75, 'HandleVisibility','off')
plot(t_s_is( nir, nf), -10, 'ko-', 'LineWidth',1, 'MarkerFaceColor','k', 'MarkerSize',4)
xline(t_s_is( nir, nf), 'k', 'HandleVisibility','off')
plot(M( nir, nf), -10,'o', 'color', cMap(:, 1).', 'LineWidth',2, 'MarkerFaceColor',cMap(:, 1).', 'MarkerSize',3)
plot(M_5( nir, nf), -20,'o', 'color', cMap(:, 6).', 'LineWidth',2, 'MarkerFaceColor',cMap(:, 6).', 'MarkerSize',3)
plot(D_E( nir, nf), -20,'o', 'color', cy.', 'LineWidth',2, 'MarkerFaceColor',cy.', 'MarkerSize',3)
plot(E( nir, nf), -15,'o', 'color', co.', 'LineWidth',2, 'MarkerFaceColor',co.', 'MarkerSize',3)

xlim([-0.01 0.7])
ylim([-60 2])

xlabel({'Time (s)', '(e)'}, 'Interpreter','latex')
% ylabel('Signal value (dB)', 'Interpreter','latex')

box on


xlabel({'Time (s)', '(c)'}, 'Interpreter','latex', 'fontsize', 11)
ylabel('Signal value (dB)', 'Interpreter','latex', 'fontsize', 11)
set(gca, 'fontsize', 11)
lgd = legend('Proposed', '$M$', '$M_5$', '$D_E$','$E$','location', 'southeast', 'interpreter', 'latex', 'numcolumns', 2, 'fontsize', 8);
%%
ff.Position(3) = 1200
ff.Position(4) = 250
%% print figure
set(ff,'Units','Inches');
set(ff,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[ff.Position(3), ff.Position(4)])
print(ff,'Onset_evaluation','-dpdf','-r0')

%% Lundeby failure rate
100 * length(inter_time_mrtd(inter_time_mrtd==0))/(size(inter_time_mrtd,1)*size(inter_time_mrtd,2))