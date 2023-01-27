%% Single-cell dose response simulator
clear; clc; close all; %set_plot_defaults_export();


pop_ec50 = 1e-8;
pop_e0 = 0.9;
pop_emax = 0.2;
pop_hs = 2^0.5;

dose_conc = logspace(-11, -5, 12);
dose_conc_fine = logspace(-11, -5, length(dose_conc)*100);
pop_dr = pop_emax + (pop_e0-pop_emax)./(1+(dose_conc./pop_ec50).^pop_hs)+ random('Normal', 0, 0.01, length(dose_conc), 1)';

figure()
semilogx(dose_conc, pop_dr, 'k')
ylim([0 1])
ylabel('Drug Effect')
xlabel('Dose (nM)')

%% Simulate single cells
N_cells = 1000;

%pm stands for plus or minus
pm_ec50 = [0.01:0.2:3.4];
pm_ec50 = -3:0.05:3;
pm_hs = [0.01:0.05:3.5];
pm_hs = -3:0.2:3.1;
pm_e0 = [0.05];
pm_emax = 0.05;
%Create new uniform distributions centered on population parameters
% dist_ec50 = makedist('Uniform', log10(pop_ec50)-pm_ec50, log10(pop_ec50)+pm_ec50);
% dist_e0 = makedist('Uniform', pop_e0 - pm_e0, min(pop_e0+pm_e0, 1));
% dist_emax = makedist('Uniform', max(pop_emax-pm_emax,0), min(pop_emax+pm_emax, 1));
% dist_hs = makedist('Uniform', log2(pop_hs)-pm_hs, log2(pop_hs)+pm_hs);
for ii = 1:length(pm_ec50)
    for jj = 1:length(pm_hs)
        %         sc_dist_ec50 = makedist('Uniform', log10(pop_ec50)-pm_ec50(ii), log10(pop_ec50)+pm_ec50(ii));
        %         sc_dist_hs = makedist('Uniform', log2(pop_hs)-pm_hs(jj), log2(pop_hs)+pm_hs(jj));
        
        sc_dist_ec50 = makedist('Normal', log10(pop_ec50) + pm_ec50(ii), 0.5);
        sc_dist_hs = makedist('Normal', log2(pop_hs) + pm_hs(jj), 0.5);
        
        sc_ec50 = 10.^random(sc_dist_ec50, N_cells,1);
        sc_e0 = random('Uniform', pop_e0 - pm_e0, pop_e0+pm_e0);
        %         sc_e0 = pop_e0;
        sc_emax = random('Uniform', pop_emax - pm_emax, pop_emax+pm_emax);
        %         sc_emax = pop_emax;
        sc_hs = 2.^random(sc_dist_hs, N_cells,1);
        
        sc_dr(ii,jj,:,:) = transpose(sc_emax + ((sc_e0-sc_emax)./(1+(dose_conc./sc_ec50).^sc_hs)))+ random('Normal', 0, 0.01, length(dose_conc), 1);
    end
end

sc_dr_flat = reshape(sc_dr,  [length(pm_ec50)*length(pm_hs), length(dose_conc), N_cells]);
sim_pop_dr = squeeze(mean(sc_dr, 4));
sim_pop_dr_flat = reshape(sim_pop_dr, [length(pm_ec50)*length(pm_hs), length(dose_conc)]);
figure()
semilogx(dose_conc, sim_pop_dr_flat, 'LineWidth', 1)
hold on
semilogx(dose_conc, pop_dr, 'k')
ylim([0 1])
ylabel('Drug Effect')
xlabel('Dose (nM)')

%% plot single-cell responses
figure()
to_plot = squeeze(sc_dr(1,1,:,:))';
auc_sort = trapz(to_plot, 2);
sorted = sortrows([auc_sort, to_plot], 1, 'descend');
imagesc(sorted(:,2:end))
ylabel('cells')
xlabel('Increasing Dose ->')
colorbar()

for ii = 1:length(pm_ec50)
    for jj = 1:length(pm_hs)
        resid(ii,jj) = sum(abs(squeeze(sim_pop_dr(ii,jj,:)) - pop_dr'));
    end
end
figure()
[X,Y] = meshgrid(pm_hs,pm_ec50);
contourf(X,Y,resid)
xlabel('HS mean')
ylabel('EC50 mean')
cb = colorbar;
ylabel(cb, 'Residual with population dose response')
hold on
scatter(1,1,'r', 'filled')

%% Isolate specific points and plot
pts = [[-2 2.6]; [1,0];[2,-2];[0, 0]];
colors = {'r', 'g', 'm', 'c'};
figure()
subplot(2,1,1)
[X,Y] = meshgrid(pm_hs,pm_ec50);
contourf(X,Y,resid)
xlabel('HS mean')
ylabel('EC50 mean')
cb = colorbar;
ylabel(cb, 'Residual with population dose response')
subplot(2,1,2)
semilogx(dose_conc, pop_dr, 'k')
ylim([0 1])
ylabel('Drug Effect')
xlabel('Dose (nM)')
for ii = 1:length(pts)
    row_idx = find(pm_ec50 == pts(ii,2));
    col_idx = find(pm_hs == pts(ii,1));
    subplot(2,1,1)
    hold on
    scatter(pts(ii,1), pts(ii,2),250, colors{ii}, 'filled')
    subplot(2,1,2)
    hold on
    semilogx(dose_conc, squeeze(sim_pop_dr(row_idx, col_idx,:)), 'Color', colors{ii})
    
end

figure()
to_plot = squeeze(sc_dr(row_idx,col_idx,:,:))';
auc_sort = trapz(to_plot, 2);
sorted = sortrows([auc_sort, to_plot], 1, 'descend');
imagesc(sorted(:,2:end))
ylabel('cells')
xlabel('Increasing Dose ->')
colorbar()

%% Try out uniform EC50 distribution w/ sharp switch SC behavior
N_cells = 1000;
sc_dist_ec50 = makedist('Uniform', log10(pop_ec50)-1.5, log10(pop_ec50) + 1.5);
sc_dist_hs = makedist('Normal', log2(pop_hs)+1, 1);

sc_ec50 = 10.^random(sc_dist_ec50, N_cells,1);
sc_e0 = random('Uniform', pop_e0 - pm_e0, pop_e0+pm_e0);
%         sc_e0 = pop_e0;
sc_emax = random('Uniform', pop_emax - pm_emax, pop_emax+pm_emax);
%         sc_emax = pop_emax;
sc_hs = 2.^random(sc_dist_hs, N_cells,1);

sc_dr = transpose(sc_emax + ((sc_e0-sc_emax)./(1+(dose_conc_fine./sc_ec50).^sc_hs)));
to_plot = sc_dr';
auc_sort = trapz(to_plot, 2);
sorted = sortrows([auc_sort, to_plot], 1, 'descend');
imagesc(sorted(:,2:end))
ylabel('cells')
xlabel('Increasing Dose ->')
colorbar()

figure()
semilogx(dose_conc_fine,mean(sc_dr'), 'k')
ylim([0 1])
ylabel('Drug Effect')
xlabel('Dose (nM)')
hold on
semilogx(dose_conc_fine, sc_dr(:,1:4), 'k--')
legend('Population dose response', 'Single cells')

x0 = [pop_emax pop_e0, log10(pop_ec50), log2(pop_hs)];
F = @(x, xdata)x(1) + (x(2)-x(1))./(1 + (xdata./10^x(3)).^(2^x(4)));
lb = [0 0 -14 -4];
ub  = [1 1 -3 4];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,dose_conc_fine,mean(sc_dr'), lb, ub);
figure()
subplot(1,2,1)
semilogx(dose_conc_fine,mean(sc_dr'), 'k')
ylim([0 1])
ylabel('Drug Effect')
xlabel('Dose (nM)')
hold on
semilogx(dose_conc_fine, sc_dr(:,1:8), 'k--', 'LineWidth', 0.5)
xlim([min(dose_conc_fine) max(dose_conc_fine)])
legend('Population dose response', 'Single cells')
subplot(1,2,2)
X = [log10(sc_ec50), log2(sc_hs)];
hist3(X,'CdataMode','auto')
xlabel('Log_{10}(EC_{50})')
ylabel('Log_{2}(HS)')
cb = colorbar()
ylabel(cb, 'Number of cells')
view(2)
hold on
scatter3(x(3), x(4), 100, 100, 'r', 'filled')
xlim([min(log10(sc_ec50)) max(log10(sc_ec50))]);
ylim([min(log2(sc_hs)) max(log2(sc_hs))]);


fprintf('Size first switch compare fig')
% pause
% saveas(gcf, ['figuresv2/switchCompareAndFit.png'])
% exportgraphics(gcf, ['figuresv2/switchCompareAndFit.pdf'], 'ContentType','vector') 


figure()
subplot(1,3,1)
semilogx(dose_conc_fine, mean(sc_dr'), 'b')
ylim([0 1])
ylabel('Drug Effect')
xlabel('Dose (nM)')
hold on
semilogx(dose_conc_fine, sc_dr(:,1:10), 'r', 'LineWidth', 1)
legend('Population dose response', 'Single cells', 'Location', 'southoutside')
xlim([min(dose_conc_fine) max(dose_conc_fine)])
subplot(1,3,2)
[f, xi] = ksdensity(log10(sc_ec50));
meanEC50 = x(3);
plot(xi, f, 'r')
xline(meanEC50, 'b', 'LineWidth', 1)
xlim([min(log10(sc_ec50)) max(log10(sc_ec50))]);
xlim([log10(pop_ec50)-1.5 log10(pop_ec50) + 1.5]);
xlabel('Log_{10}(EC_{50})')
legend('Single-cell Distribution', 'Population Fit', 'Location', 'southoutside')

subplot(1,3,3)
[f, xi] = ksdensity(log2(sc_hs));
meanHS = x(4);
plot(xi, f, 'r')
xline(meanHS, 'b', 'LineWidth', 1)
xlim([min(log2(sc_hs)) max(log2(sc_hs))]);
legend('Single-cell Distribution', 'Population Fit', 'Location', 'southoutside')
xlabel('Log_{2}(HS)')

fprintf('Size first switch compare fig')
% pause
% saveas(gcf, ['figuresv2/switchCompareAndFit2.png'])
% exportgraphics(gcf, ['figuresv2/switchCompareAndFit2.pdf'], 'ContentType','vector') 
%% Compare 2 distributions with different emax

N_cells = 1000;

%pm stands for plus or minus
pm_ec50 = [0.01:0.2:3.4];
pm_ec50 = -3:0.05:3;
pm_hs = [0.01:0.05:3.5];
pm_hs = -3:0.2:3.1;
pm_e0 = [0.05];
pm_emax = 0.05;

%Create new uniform distributions centered on population parameters
% dist_ec50 = makedist('Uniform', log10(pop_ec50)-pm_ec50, log10(pop_ec50)+pm_ec50);
% dist_e0 = makedist('Uniform', pop_e0 - pm_e0, min(pop_e0+pm_e0, 1));
% dist_emax = makedist('Uniform', max(pop_emax-pm_emax,0), min(pop_emax+pm_emax, 1));
% dist_hs = makedist('Uniform', log2(pop_hs)-pm_hs, log2(pop_hs)+pm_hs);
sc_dr = zeros(length(pm_ec50), length(pm_hs), length(dose_conc), N_cells);
for ii = 1:length(pm_ec50)
    for jj = 1:length(pm_hs)
        %         sc_dist_ec50 = makedist('Uniform', log10(pop_ec50)-pm_ec50(ii), log10(pop_ec50)+pm_ec50(ii));
        %         sc_dist_hs = makedist('Uniform', log2(pop_hs)-pm_hs(jj), log2(pop_hs)+pm_hs(jj));
        
        sc_dist_ec50 = makedist('Normal', log10(pop_ec50) + pm_ec50(ii), 0.5);
        sc_dist_hs = makedist('Normal', log2(pop_hs) + pm_hs(jj), 0.5);
        
        sc_ec50 = 10.^random(sc_dist_ec50, N_cells,1);
        sc_e0 = random('Uniform', pop_e0 - pm_e0, pop_e0+pm_e0);
        %         sc_e0 = pop_e0;
        sc_emax = random('Uniform', pop_emax - pm_emax, pop_emax+pm_emax);
        %         sc_emax = pop_emax;
        sc_hs = 2.^random(sc_dist_hs, N_cells,1);
        
        sc_dr(ii,jj,:,:) = transpose(sc_emax + ((sc_e0-sc_emax)./(1+(dose_conc./sc_ec50).^sc_hs)))+ random('Normal', 0, 0.01, length(dose_conc), 1);
    end
end

sc_dr_flat = reshape(sc_dr,  [length(pm_ec50)*length(pm_hs), length(dose_conc), N_cells]);
sim_pop_dr = squeeze(mean(sc_dr, 4));
sim_pop_dr_flat = reshape(sim_pop_dr, [length(pm_ec50)*length(pm_hs), length(dose_conc)]);
figure()
semilogx(dose_conc, sim_pop_dr_flat, 'LineWidth', 1)
hold on
semilogx(dose_conc, pop_dr, 'k')
ylim([0 1])
ylabel('Drug Effect')
xlabel('Dose (nM)')


