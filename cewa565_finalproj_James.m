clc; clear; close all
load('final_project_data.mat')

bin = 63;
%% kill major spikes
r = find(SigTKE(:,bin) > 10e-6);
SigTKE(r,bin) = nan;
r = find(Sig_sediment_conc(:,bin) > 43);
Sig_sediment_conc(r,bin) = nan;
SigAmp(r,bin) = nan;

%% 3 month timeseries
figure()
subplot(5,1,1)
plot(Time, Sig_sediment_conc(:,bin))
title('Sediment Concentration')
ylabel('Normalized, %')

subplot(5,1,2)
plot(Time, SigAmp(:,bin))
title('ADCP Corrected Amplitude')
ylabel('dB')

% w velocity
subplot(5,1,3)
plot(Time, VecAmp(:,3))
title('ADV Amplitude')
ylabel('dB')

subplot(5,1,4)
plot(Time, SigTKE(:,bin))
title('ADCP Turbulence')
ylabel('Vertical TKE, m^2/s^2')

subplot(5,1,5)
plot(Time, VecTKE)
title('ADV Turbulence')
ylabel('TKE, m^2/s^2')

%% Regressions
fprintf("Sediment Conc vs ADCP TKE regression\n")
figure()
x = SigTKE(:,bin)*1e5;
y = Sig_sediment_conc(:,bin);
scatter(x, y)
xlabel('Signature1000 vertical TKE (m^2/s^2 x10^-5)')
ylabel('Signature1000 Sediment Conc (%)')
hold on
%regression(x, y, 1.96, Time);
figure()
%heatmaps
mpt = linspace(min(x),max(x),100);
npt = linspace(min(y),max(y),100);
N = histcounts2(y(:),x(:),npt,mpt);
imagesc(mpt,npt,N)
myColorMap = jet(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar
set(gca, 'XLim', mpt([1 end]), 'YLim', npt([1 end]), 'YDir', 'normal');

fprintf("Sediment Conc vs ADV TKE regression\n")
figure()
x = VecTKE*1e5;
y = Sig_sediment_conc(:,bin);
scatter(x, y)
xlabel('Vector TKE (m^2/s^2 x10^-5)')
ylabel('Signature1000 Sediment Conc (%)')
hold on
%regression(x, y, 1.96, Time);
figure()
%heatmaps
mpt = linspace(min(x),max(x),100);
npt = linspace(min(y),max(y),100);
N = histcounts2(y(:),x(:),npt,mpt);
imagesc(mpt,npt,N)
myColorMap = jet(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar
set(gca, 'XLim', mpt([1 end]), 'YLim', npt([1 end]), 'YDir', 'normal');

fprintf("ADCP vs ADV TKE regression\n")
figure()
x = VecTKE;
y = SigTKE(:,bin);
scatter(x, y)
xlabel('Vector TKE (m^2/s^2)')
ylabel('Signature TKE (m^2/s^2)')
hold on
%regression(x, y, 1.96, Time);
figure()
%heatmaps
mpt = linspace(min(x),max(x),100);
npt = linspace(min(y),max(y),100);
N = histcounts2(y(:),x(:),npt,mpt);
imagesc(mpt,npt,N)
myColorMap = jet(256);
myColorMap(1,:) = 1;
colormap(myColorMap);
colorbar
set(gca, 'XLim', mpt([1 end]), 'YLim', npt([1 end]), 'YDir', 'normal');

%% Autocorrelation
n = 100;
[SigTKE_acf,~,SigTKE_bounds] = autocorr(double(SigTKE(:,bin)), 'NumLags', n);
[VecTKE_acf,~,VecTKE_bounds] = autocorr(VecTKE, 'NumLags', n);

figure()
plot(SigTKE_acf); hold on
plot([0 n], [SigTKE_bounds(1) SigTKE_bounds(1)], 'r', [0 n], [SigTKE_bounds(2) SigTKE_bounds(2)], 'r')
ylabel('Correlation Coeff')
title('ADCP TKE Autocorrelation')
legend('autocor Rs','95% bounds');

figure()
plot(VecTKE_acf); hold on
plot([0 n], [VecTKE_bounds(1) VecTKE_bounds(1)], 'r', [0 n], [VecTKE_bounds(2) VecTKE_bounds(2)], 'r')
ylabel('Correlation Coeff')
title('ADV TKE Autocorrelation')
legend('autocor Rs','95% bounds');

%% EOFs and PCs, Aug18-Nov18
[U,S,V] = svd(Sig_sediment_conc,0); % timeseries with time on row axis, so U is the PC and V is EOF

n_bins = size(Sig_sediment_conc,2);
var_exp = NaN(n_bins,1);
for i = 1:n_bins
	var_exp(i) = S(i,i)/sum(S(:));
end
figure()
plot(1:n_bins,var_exp)
xlabel('EOF number')
ylabel('variance explained')
print('65% of the variance is explained in first pattern');

% first pattern
figure(); subplot(2,1,1)
plot(1:n_bins, -V(1,:));
xlabel('depth bin')
ylabel('spatial weight')
title('EOF 1')

subplot(2,1,2)
plot(Time, -U(:,1))
xlabel('time')
datetick('x')
ylabel('temporal weight')
title('PC 1')

% second pattern
figure(); subplot(2,1,1)
plot(1:n_bins, V(2,:));
xlabel('depth bin')
ylabel('spatial weight')
title('EOF 2')

subplot(2,1,2)
plot(Time, U(:,2))
xlabel('time')
datetick('x')
ylabel('temporal weight')
title('PC 2')
