clc; clear; close all
data = load('Vector_LakeK_8.16.18_11.5.18_ALL.mat');
idx = isnan(data.u);
u = data.u(~idx);
idy = isnan(data.v);
v = data.v(~idy);
idz = isnan(data.w);
w = data.w(~idz);

time = data.time;
flowTime = datetime(time,'ConvertFrom','datenum');
% time doesn't match up at all

% despiking
[uc, vc, wc] = func_despike_phasespace3d_3var(u, v, w, 2);

fs = 16; % ADV sounded at 16Hz
figure
plot(uc);
xlabel('time');
ylabel('velocity (m/s)');
title('u raw data');
figure
plot(vc);
xlabel('time');
ylabel('velocity (m/s)');
title('v raw data');
figure
plot(wc);
xlabel('time');
ylabel('velocity (m/s)');
title('w raw data');

%%
% Split data into ensemble windows (need to find variances)
window = 3*8192; % datapoints/hour, ensemble time windows

n = floor(length(u)/window);
ensemble_u = zeros(n,window); ensemble_v = zeros(n,window); ensemble_w = zeros(n,window);
for i=1:n
    ensemble_u(i,:) = uc(1+window*(i-1):window+window*(i-1));
    ensemble_v(i,:) = vc(1+window*(i-1):window+window*(i-1));
    ensemble_w(i,:) = wc(1+window*(i-1):window+window*(i-1));
end

% plotting ensemble means
u_mean = mean(ensemble_u,2);
v_mean = mean(ensemble_v,2);
w_mean = mean(ensemble_w,2);
figure
plot(u_mean);
xlabel('timestep');
ylabel('velocity (m/s)');
title('u ensemble means');
figure
plot(v_mean);
xlabel('timestep');
ylabel('velocity (m/s)');
title('v ensemble means');
figure
plot(w_mean);
xlabel('timestep');
ylabel('velocity (m/s)');
title('w ensemble means');

% ensemble variances
u_var = var(ensemble_u,0,2);
v_var = var(ensemble_v,0,2);
w_var = var(ensemble_w,0,2);

% take means out of ensembles (detrending for spectral density stuff)
ensemble_u = ensemble_u - u_mean;
ensemble_v = ensemble_v - v_mean;
ensemble_w = ensemble_w - w_mean;

%% Using pwelch to fft data into a power spectral density
for i=1:n
    [psd_u(:,i), f(:,i)] = pwelch(ensemble_u(i,:),[],[],[],fs);
    [psd_v(:,i)] = pwelch(ensemble_v(i,:),[],[],[],fs);
    [psd_w(:,i)] = pwelch(ensemble_w(i,:),[],[],[],fs);
end
figure
loglog(f,psd_u);
xlabel('Frequency (Hz)');
ylabel('PSD');
title('u dir');
figure
loglog(f,psd_v);
xlabel('Frequency (Hz)');
ylabel('PSD');
title('v dir');
figure
loglog(f,psd_w);
xlabel('Frequency (Hz)');
ylabel('PSD');
title('w dir');

%% Integrating the spectra
% % subtract out the white noise bits and integrating the rest should be turbulent variance
% % so area_um = var_u
% % find where f = 10Hz and subtract out the integration of that spectral density
% r = find(f(:,1)==10);
% area_u = trapz(f(:,1),psd_u) - trapz(f(r:end,1),psd_u(r:end,:));
% r = find(f(:,1)==12);
% area_v = trapz(f(:,1),psd_v) - trapz(f(r:end,1),psd_v(r:end,:));
% r = find(f(:,1)==15);
% area_w = trapz(f(:,1),psd_w) - trapz(f(r:end,1),psd_w(r:end,:));

%% TKE u'w'bar
ReStress = zeros(3,3,n);
for i=1:n
    u_prime = ensemble_u(i,:);
    v_prime = ensemble_v(i,:);
    w_prime = ensemble_w(i,:);
    u_prime_bar = mean(u_prime.^2); % u'^2bar
    v_prime_bar = mean(v_prime.^2); % v'^bar
    w_prime_bar = mean(w_prime.^2); % w'^2bar
    uv_prime_bar = mean(u_prime.*v_prime); % u'v'bar
    uw_prime_bar(i) = mean(u_prime.*w_prime); % u'w'bar (m^2/s^2) *10000; % cm^2/s^2
    vw_prime_bar = mean(v_prime.*w_prime); % v'w'bar
    
    ReStress(:,:,i) = [u_prime_bar  uv_prime_bar uw_prime_bar(i);...
                       uv_prime_bar v_prime_bar  vw_prime_bar;...
                       uw_prime_bar(i) vw_prime_bar w_prime_bar]; % m^2/s^2
    TKE(i) = .5*(u_prime_bar + v_prime_bar + w_prime_bar);
end
figure
plot(uw_prime_bar);
ylabel('u''w''bar (m^2/s^2)');
xlabel('timestep');
title('Reynolds Stress');

