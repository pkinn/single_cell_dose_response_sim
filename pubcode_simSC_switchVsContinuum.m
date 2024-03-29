%% Single-cell dose response simulator
clear; clc; close all; 

%Input dose response parameters and generate population dose response

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



%% Simulate populations with variable EC50 and sharp switch behavior
N_cells = 1000;

%Create distributions of ec50, hill slope, e0, and emax parameters
pm_e0 = 0.05;
pm_emax = 0.05;

sc_dist_ec50 = makedist('Uniform', log10(pop_ec50)-1.5, log10(pop_ec50) + 1.5);
sc_dist_hs = makedist('Normal', log2(pop_hs)+1, 1);

sc_ec50 = 10.^random(sc_dist_ec50, N_cells,1);
sc_e0 = random('Uniform', pop_e0 - pm_e0, pop_e0+pm_e0);
sc_emax = random('Uniform', pop_emax - pm_emax, pop_emax+pm_emax);
sc_hs = 2.^random(sc_dist_hs, N_cells,1);
%Create single-cell dose response matrix with dimensions [N_doses x
%N_cells];
sc_dr = transpose(sc_emax + ((sc_e0-sc_emax)./(1+(dose_conc_fine./sc_ec50).^sc_hs)));

%Fit dose response curve to population average dose response
x0 = [pop_emax pop_e0, log10(pop_ec50), log2(pop_hs)];
F = @(x, xdata)x(1) + (x(2)-x(1))./(1 + (xdata./10^x(3)).^(2^x(4)));
lb = [0 0 -14 -4];
ub  = [1 1 -3 4];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,dose_conc_fine,mean(sc_dr'), lb, ub);


% Plot 3 panels:
%Panel 1: Mean and single-cell dose responses
figure()
subplot(1,3,1)
semilogx(dose_conc_fine, mean(sc_dr'), 'b')
ylim([0 1])
ylabel('Drug Effect')
xlabel('Dose (nM)')
hold on
%Plot random single cell dose responses (There's no order to the sc_dr
%matrix so first 10 dose responses are random). 
semilogx(dose_conc_fine, sc_dr(:,1:10), 'r', 'LineWidth', 1)
legend('Population dose response', 'Single cells', 'Location', 'southoutside')
xlim([min(dose_conc_fine) max(dose_conc_fine)])
%Panel 2: Plot distribution of EC50 values and mean from population fit
subplot(1,3,2)
[f, xi] = ksdensity(log10(sc_ec50));
meanEC50 = x(3);
plot(xi, f, 'r')
xline(meanEC50, 'b', 'LineWidth', 1)
xlim([min(log10(sc_ec50)) max(log10(sc_ec50))]);
xlim([log10(pop_ec50)-1.5 log10(pop_ec50) + 1.5]);
xlabel('Log_{10}(EC_{50})')
legend('Single-cell Distribution', 'Population Fit', 'Location', 'southoutside')
%Panel 3: Plot distribution of HS values and mean from population fit. 
subplot(1,3,3)
[f, xi] = ksdensity(log2(sc_hs));
meanHS = x(4);
plot(xi, f, 'r')
xline(meanHS, 'b', 'LineWidth', 1)
xlim([min(log2(sc_hs)) max(log2(sc_hs))]);
legend('Single-cell Distribution', 'Population Fit', 'Location', 'southoutside')
xlabel('Log_{2}(HS)')

%Save figures (optional)
% saveas(gcf, ['figures/switchCompareAndFit2.png'])
% exportgraphics(gcf, ['figures/switchCompareAndFit2.pdf'], 'ContentType','vector') 
