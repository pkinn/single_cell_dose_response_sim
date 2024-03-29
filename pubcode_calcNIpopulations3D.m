%% Single-cell dose response simulator - comparing populations and generating threshold inhibition plots
clear; clc; close all; set_plot_defaults_export();
% set_plot_defaults
omiPopEC50 = 10^-7.7246;
omiPopE0 = 1;
omiPopEmax = 0.1585;
omiPopHS = 2^-1.024;

doseConc = logspace(-11, -5, 12);
omiPopDR = omiPopEmax + (omiPopE0-omiPopEmax)./(1+(doseConc./omiPopEC50).^omiPopHS);%+ random('Normal', 0, 0.01, length(doseConc), 1)';

figure()
semilogx(doseConc, omiPopDR, 'k')
ylim([0 1])
ylabel('Drug Effect')
xlabel('Dose (nM)')

%% Simulate four cases that all converge to the same population DR:
% low heterogeneity, high heterogeneity, bimodal/two populations, and
% resistant subpopulation

%For each population, make distributions and generate single-cell dose
%responses and plot.

%Population 1: Low heterogeneity
pdfEC50 = makedist('Normal', 'mu', log10(omiPopEC50), 'sigma', 0.1);
pdfEmax = makedist('Normal', 'mu', omiPopEmax, 'sigma', 0.1);
pdfHS = makedist('Normal', 'mu', log2(omiPopHS), 'sigma', 0.1);
pdfE0 = makedist('Normal', 'mu', 1, 'sigma', 0.1);
nCells = 1000;

params1 = [10.^random(pdfEC50, nCells,1) random(pdfE0,nCells,1) random(pdfEmax, nCells,1) 2.^random(pdfHS, nCells,1)];


scDR1 = zeros(nCells, length(doseConc));
for ii = 1:nCells
    scDR1(ii,:) = returnDR(params1(ii,:),doseConc);
end

scDR1mean = mean(scDR1, 1);

figure()
semilogx(doseConc, scDR1, 'LineWidth', 0.1)
hold on
semilogx(doseConc, omiPopDR, 'k')
semilogx(doseConc, scDR1mean, 'k--')

%Population 2: high heterogeneity
pdfEC50 = makedist('Normal', 'mu', log10(omiPopEC50), 'sigma', 0.5);
pdfEmax = makedist('Normal', 'mu', omiPopEmax, 'sigma', 0.1);
pdfHS = makedist('Normal', 'mu', log2(omiPopHS), 'sigma', 0.5);
pdfE0 = makedist('Normal', 'mu', 1, 'sigma', 0.1);
nCells = 1000;

params2 = [10.^random(pdfEC50, nCells,1) random(pdfE0,nCells,1) random(pdfEmax, nCells,1) 2.^random(pdfHS, nCells,1)];


scDR2 = zeros(nCells, length(doseConc));
for ii = 1:nCells
    scDR2(ii,:) = returnDR(params2(ii,:),doseConc);
end

scDR2mean = mean(scDR2, 1);

figure()
semilogx(doseConc, scDR2, 'LineWidth', 0.1)
hold on
semilogx(doseConc, omiPopDR, 'k')
semilogx(doseConc, scDR2mean, 'k--')


%Population 3: Bimodal 
pdfEC50 = makedist('Normal', 'mu', log10(omiPopEC50), 'sigma', 0.5);
pdfEC50 = gmdistribution([log10(omiPopEC50)-1 log10(omiPopEC50)+1]', 0.1);
pdfEmax = makedist('Normal', 'mu', omiPopEmax, 'sigma', 0.1);
pdfHS = makedist('Normal', 'mu', log2(omiPopHS), 'sigma', 0.5);
pdfHS = gmdistribution([log2(omiPopHS)-0.4 log2(omiPopHS)+0.4]', 0.1);
pdfE0 = makedist('Normal', 'mu', 1, 'sigma', 0.1);
nCells = 1000;

params3 = [10.^random(pdfEC50, nCells) random(pdfE0,nCells,1) random(pdfEmax, nCells,1) 2.^random(pdfHS, nCells)];


scDR3 = zeros(nCells, length(doseConc));
for ii = 1:nCells
    scDR3(ii,:) = returnDR(params3(ii,:),doseConc);
end

scDR3mean = mean(scDR3, 1);

figure()
semilogx(doseConc, scDR3, 'LineWidth', 0.1)
hold on
semilogx(doseConc, omiPopDR, 'k')
semilogx(doseConc, scDR3mean, 'k--')


%Population 4: Resistant subpopulation
pdfEC50 = makedist('Normal', 'mu', log10(omiPopEC50), 'sigma', 0.1);
pdfEmax = makedist('Normal', 'mu', omiPopEmax, 'sigma', 0.1);
pdfHS = makedist('Normal', 'mu', log2(omiPopHS), 'sigma', 0.1);
pdfE0 = makedist('Normal', 'mu', 1, 'sigma', 0.1);
pdfEC50Res = makedist('Normal', 'mu', log10(omiPopEC50)+3, 'sigma', 0.1);
pdfEmaxRes = makedist('Normal', 'mu', omiPopEmax+0.3, 'sigma', 0.1);


nCellsNormal = floor(0.95*nCells);
nCellsRes = floor(0.05*nCells);
params4Normal = [10.^random(pdfEC50, nCellsNormal,1) random(pdfE0,nCellsNormal,1) random(pdfEmax, nCellsNormal,1) 2.^random(pdfHS, nCellsNormal,1)];
params4Res = [10.^random(pdfEC50Res, nCellsRes,1) random(pdfE0,nCellsRes,1) random(pdfEmaxRes, nCellsRes,1) 2.^random(pdfHS, nCellsRes,1)];
params4 = [params4Normal; params4Res];

scDR4 = zeros(nCells, length(doseConc));
for ii = 1:nCells
    scDR4(ii,:) = returnDR(params4(ii,:),doseConc);
end

scDR4mean = mean(scDR4, 1);

figure()
semilogx(doseConc, scDR4, 'LineWidth', 0.1)
hold on
semilogx(doseConc, omiPopDR, 'k')
semilogx(doseConc, scDR4mean, 'k--')

scDRMeanAll = [scDR1mean; scDR2mean; scDR3mean; scDR4mean];
paramsAll = cat(3, params1, params2, params3, params4);
scDRAll = cat(3, scDR1, scDR2, scDR3, scDR4);

%% Calculate nInh
%For a specific threshold (0.6) calculate the number of cells inhibited
%below that threshold at each dose.
doseVecFine = logspace(min(log10(doseConc)), max(log10(doseConc)), 1000);

for ii = 1:size(paramsAll,3)
    nInh1D(ii,:) = calculatenInhThresh(paramsAll(:,:,ii), 0.6, doseVecFine);
end
%% Make plots for threshold inhibition curves
figure('Position',[100 100 850 450])
subplot(1,2,2)
semilogx(doseVecFine, nInh1D./nCells)
legText = { 'Low Variance', 'High Variance', 'Bimodal', 'Resistant Subpop.'};

xlabel('Dose (nM)')
ylabel('% of cells < 60% basal activity')
% legend(legText)
subplot(1,2,1)
semilogx(doseConc, scDRMeanAll)
ylabel('Mean dose response')
xlabel('Dose (nM)')
legend(legText, 'Location', 'SouthWest')
% exportgraphics(gcf, [ 'hypotheticalDRvsNIcompare.pdf'], 'ContentType','vector')
% exportgraphics(gcf, [ 'hypotheticalDRvsNIcompare.png'], 'Resolution', 1000)

%% Calculate nInh for different thresholds for threshold inhibition surfaces
threshRange = 0.2:0.001:0.9;
for ii = 1:size(paramsAll,3)
    for jj = 1:length(threshRange)
        nInh(ii,:,jj) = calculatenInhThresh(paramsAll(:,:,ii), threshRange(jj), doseVecFine);
        if mod(jj,100) == 0
            fprintf('jj = %i/%i', jj, length(threshRange))
        end
    end
end


%% Plot threshold inhibition surfaces
logDoseVecFine = log10(doseVecFine);
[dim1Grid, dim2Grid] = meshgrid(logDoseVecFine, 1-threshRange);
figure('Position',[100 100 900 900])
tl = tiledlayout(2,2);
for ii = 1:size(paramsAll,3)
    h(ii) = nexttile(tl);
    contourf(dim1Grid, (dim2Grid), squeeze(nInh(ii,:,:))'./nCells)
    caxis([0 1])
    title(legText{ii})
end
xlabel(tl,'Log_{10}(Dose)', 'FontSize', 20)
ylabel(tl,'Inhibition threshold', 'FontSize', 20)
cbh = colorbar(h(end));
% To position the colorbar as a global colorbar representing
% all tiles,
cbh.Layout.Tile = 'east';
ylabel(cbh, 'Percent below threshold')

% exportgraphics(gcf, [ 'hypotheticalDRsurfacesHighRes.pdf'], 'ContentType','vector')
% exportgraphics(gcf, [ 'hypotheticalDRsurfacesHighRes.png'], 'Resolution', 1000)



