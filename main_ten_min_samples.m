%% Ten minute turbulence vs C spatial and temporal variation
clc; clear; close all

load('data/S100882A003_LakeK_Sig2_31.mat')
% grab a ten minute burst (8 Hz shots * 512 shots/10min burst = 4096 rows)
n = 4096;
bottom_bin = 72;
lakebed = 71;

%for m=0:9
m = 6;
hour = m*n+1:(m+1)*n;
dir = ['./tenMinuteSamples/Sig31_hour_#' num2str(m)];
%mkdir(dir);

time = Data.BurstHR_Time(hour);
flowTime = datetime(time,'ConvertFrom','datenum');

% average amplitudes into 1 second ensembles
amp_dB = Data.BurstHR_AmpBeam5(hour,1:bottom_bin);
f = 8; %Hz
ten = length(amp_dB)/f;
dB = zeros(ten, bottom_bin);
for i=0:ten-1
    dB(i+1,:) = mean(amp_dB(i*8+1:i*8+8,:),1);
    time_avg(i+1,:) = time(i*8+1,:);
end

w = Data.BurstHR_VelBeam5(hour,1:bottom_bin);
%depth = .1 + (0:1:bottom_bin)*.04 - .02;
depth = (1:double(Data.BurstHR_NCells(1)))*Config.BurstHR_CellSize + ...
    Config.Burst_BlankingDistance - Config.BurstHR_CellSize/2; % dB range, in m
depth = -depth(1:bottom_bin);

% Despike data
w_despiked = w;
spikestd = 1;  % set higher for weaker despiking
for i = 1:bottom_bin
    % Pull out column
    w_temp=w_despiked(:,i);
    wspike = find( abs(w_temp) > (spikestd*nanstd(w_temp) + abs(nanmean(w_temp))) );
    w_temp(wspike) = NaN;
    % replace bad points with burst means
    w_temp(isnan(w_temp)) = nanmean(w_temp);
    
    % Put despiked data back into column
    w_despiked(:,i)=w_temp;
end

% initial timeseries
surface(time, depth, amp_dB', 'EdgeColor','None');
datetick('x',13,'keepticks')
dynamicDateTicks()
xlabel('time')
ylabel('depth, m')
title('amplitude')
colorbar
saveas(gcf,[dir '/initial_ts.jpg']); 

figure()
surface(time, depth, w_despiked', 'EdgeColor','None');
datetick('x',13,'keepticks')
dynamicDateTicks()
xlabel('time')
ylabel('depth, m')
title('vertical velocity, m/s')
colorbar; caxis([-.005 .005])
saveas(gcf,[dir '/w_despiked.jpg']); 

%% concentration calibration

% a = 750e-6; % particle size, m
% f = 1; % MHz
% sound = Data.BurstHR_Soundspeed(1:n); % m/s
% temp = Data.BurstHR_Temperature(1:n); % C
% c = zeros(ten, bottom_bin);
% T = zeros(ten, bottom_bin);
% for i=0:ten-1
%     c(i+1,:) = mean(sound(i*8+1:i*8+8,:),1);
%     T(i+1,:) = mean(temp(i*8+1:i*8+8,:),1);
% end
% k = 2*pi*f*1e6./c; % wave number, /m
% v = 1e-6; % viscosity of water, m^2/s
% rho_s = 1330; % kg/m^3, clay/silt/organics?
% rho_w = 1000; % kg/m^3
% 
% beta = ((pi*f*1e6)/v)^.5;
% delta = .5*(1 + 9/(beta*a));
% s = 9/(2*beta*a)*(1 + 2/(beta*a));
% sigma = rho_s/rho_w;
% e = 1;
% eta = ((e-1)/(3*e))^2 + 1/3*((sigma-1)/(2*sigma+1))^2;
% 
% % water and sediment/particle attenuation
% alpha_w = @(T) (55.9 - 2.37.*T + 4.77e-2.*T.^2 - 3.84e-4.*T.^3).*10^-3*f^2;
% alpha_s = @(k,M) ((k.^4.*a.^3)./(96.*rho_s) + (k.*(sigma-1)^2.*s)./(s*rho_s...
%     *(s^2 + (sigma + delta)^2))).*(20./log(10))*M;
% alpha = @(T,k,M) alpha_w(T) + alpha_s(k,M);
% 
% 
% M_calc = zeros(length(k),bottom_bin);
% A_corrected = zeros(length(k),bottom_bin);
% r = -depth;
% for i=1:length(k) % rows of Sv
%     for j=1:bottom_bin %length(r) % columns
%         M = .0001;
%         t = true;
%         while t==true
%             % Hoitink's claimed concentration from corrected amplitude relation
%             % M is prop to 10^(Sv/25)
%             if k(i)*a < 1 % Rayleigh scattering
%                 Sv = -10*log10((3*eta*k(i)^4)/rho_s*a^3*M);
%             else % ka>1
%             	Sv = -10*log10((3*eta*k(i)^4*a^3)/(2*rho_s*(2+3*(k(i)*a)^4))*M);
%             end
%             % Nortek's stated sonar equation to calculate corrected amplitude
%             Sv_actual = dB(i,j)*.43 + 20*log10(r(j)) + 2*alpha(T(i),k(i),M)*r(j);
%             
%             if abs(1-Sv/Sv_actual) < .001
%                 M_calc(i,j) = 1/M; % Univerting M from Hoitink Sv inversion
%                 A_corrected(i,j) = Sv;
%                 t = false;
%             else
%                 if Sv_actual > 100
%                     error(['M step too high ' num2str(a) ' ' num2str(rho_s)])
%                 end
% %                 fprintf('M = %d; ', M)
% %                 fprintf('Sv = %4.2f; ', Sv)
% %                 fprintf('Sv_actual = %4.2f\n', Sv_actual)
%                 M = M + .0001;
%             end
%         end
%     end
%     disp(i)
% end
% save([dir '/ten_minutes_of_data.mat'],'M_calc','A_corrected')

%%
load([dir '/ten_minutes_of_data.mat']);
% normalize concentration
normal_factor = max(M_calc,[],2); % set bottom conc = 1
normal_conc = M_calc./normal_factor; % normalize by bottom conc

% remove linear trend
background_conc = max(mink(normal_conc,2,1));
x = find(diff(background_conc)./max(background_conc)>.3);
if numel(x)>1
    x = min(x);
end
% setting background suspension = 0, i.e unknown C constant
final_cal_conc = normal_conc - [background_conc(1:x) zeros(1,length(background_conc)-x)];
final_avg_conc = mean(final_cal_conc,1);

% applying linear fix to corrected amplitude
background_dB = max(mink(A_corrected,2,1));
final_cor_amp = A_corrected - [background_dB(1:x) zeros(1,length(background_conc)-x)];

figure()
surface(time_avg, depth(1:lakebed), final_cal_conc(:,1:lakebed)','EdgeColor','None');
colorbar; caxis([0 .2])
title('Final Normalized Concentration')
datetick('x',13,'keepticks')
dynamicDateTicks()
xlim([min(time_avg),max(time_avg)]);
xlabel('time')
ylabel('depth, m')
saveas(gcf,[dir '/concentration.jpg']); 

figure()
surface(time_avg, depth(1:lakebed), final_cor_amp(:,1:lakebed)','EdgeColor','None');
colorbar; caxis([0 10])
title('Final Corrected Amplitude (-C dB)')
datetick('x',13,'keepticks')
dynamicDateTicks()
xlim([min(time_avg),max(time_avg)]);
xlabel('time')
ylabel('depth, m')
saveas(gcf,[dir '/corrected_amp.jpg']);

%% Vertical vertical TKE
% averaged per second
w_bar = mean(w_despiked(:));
w_prime_2 = (w_despiked - w_bar).^2;
k_w = nan*dB;
for i=0:ten-1
    k_w(i+1,:) = mean(w_prime_2(i*8+1:i*8+8,:));
end

figure()
surface(time_avg, depth(1:lakebed), k_w(:,1:lakebed)', 'EdgeColor','None');
datetick('x',13,'keepticks')
dynamicDateTicks()
title('Vertical TKE, m^2/s^2')
xlabel('time')
ylabel('depth, m')
colorbar; caxis([0 3e-5])
saveas(gcf,[dir '/reynolds_stress.jpg']);

%% tau vs C Regression
bin = 70;
x = k_w(:,bin);
y = final_cal_conc(:,bin);
figure('Renderer', 'painters', 'Position', [10 10 900 400])
yyaxis left
plot(time_avg, x);
ylabel('vertical TKE, m^2/s^2')
yyaxis right
plot(time_avg, y);
ylabel('normalized conc');
xlabel('time');
datetick('x',13,'keepticks')
dynamicDateTicks()
title(['Timeseries at z = ' num2str(depth(bin)) 'm']);
saveas(gcf,[dir '/ts.jpg']);

% figure()
% scatter(x, y);
% xlabel('Vertical vertical TKE, m^2/s^2')
% ylabel('Normalized Concentration')
% hold on
% regression(x, y, 1.96, time);

%% EOFs and PCs/ spatial temporal variation
[U,S,V] = svd(final_cal_conc(:,1:lakebed),0); % timeseries with time on row axis, so U is the PC and V is EOF

n_bins = size(final_cal_conc(:,1:lakebed),2);
var_exp = NaN(n_bins,1);
for i = 1:n_bins
	var_exp(i) = S(i,i)/sum(S(:));
end
figure()
plot(1:n_bins,var_exp)
xlabel('EOF number')
ylabel('variance explained')
saveas(gcf,[dir '/conc_variance.jpg']);

% first pattern
figure(); subplot(2,1,1)
plot(-depth(1:lakebed), -V(1,:));
xlabel('depth')
ylabel('spatial weight')
title('EOF 1')
axis([0 3 -.2 .2])

subplot(2,1,2)
plot(time_avg, -U(:,1))
xlabel('time')
datetick('x',13,'keepticks')
dynamicDateTicks()
ylabel('temporal weight')
title('PC 1')
saveas(gcf,[dir '/conc_pattern1.jpg']);

% second pattern
figure(); subplot(2,1,1)
plot(-depth(1:lakebed), V(2,:));
xlabel('depth')
ylabel('spatial weight')
title('EOF 2')
axis([0 3 -.2 .2])

subplot(2,1,2)
plot(time_avg, U(:,2))
xlabel('time')
datetick('x',13,'keepticks')
dynamicDateTicks()
ylabel('temporal weight')
title('PC 2')
saveas(gcf,[dir '/conc_pattern2.jpg']);

% third pattern
figure(); subplot(2,1,1)
plot(-depth(1:lakebed), V(3,:));
xlabel('depth')
ylabel('spatial weight')
title('EOF 3')
axis([0 3 -.2 .2])

subplot(2,1,2)
plot(time_avg, U(:,3))
xlabel('time')
datetick('x',13,'keepticks')
dynamicDateTicks()
ylabel('temporal weight')
title('PC 3')
saveas(gcf,[dir '/conc_pattern3.jpg']);

[U,S,V] = svd(k_w(:,1:lakebed),0); % timeseries with time on row axis, so U is the PC and V is EOF

n_bins = size(k_w(:,1:lakebed),2);
var_exp = NaN(n_bins,1);
for i = 1:n_bins
	var_exp(i) = S(i,i)/sum(S(:));
end
figure()
plot(1:n_bins,var_exp)
xlabel('EOF number')
ylabel('variance explained')
saveas(gcf,[dir '/tau_variance.jpg']);

% first pattern
figure(); subplot(2,1,1)
plot(-depth(1:lakebed), -V(1,:));
xlabel('depth')
ylabel('spatial weight')
title('EOF 1')
axis([0 3 -.2 .2])

subplot(2,1,2)
plot(time_avg, -U(:,1))
xlabel('time')
datetick('x',13,'keepticks')
dynamicDateTicks()
ylabel('temporal weight')
title('PC 1')
saveas(gcf,[dir '/tau_pattern1.jpg']);

% second pattern
figure(); subplot(2,1,1)
plot(-depth(1:lakebed), V(2,:));
xlabel('depth')
ylabel('spatial weight')
title('EOF 2')
axis([0 3 -.2 .2])

subplot(2,1,2)
plot(time_avg, U(:,2))
xlabel('time')
datetick('x',13,'keepticks')
dynamicDateTicks()
ylabel('temporal weight')
title('PC 2')
saveas(gcf,[dir '/tau_pattern2.jpg']);

% third pattern
figure(); subplot(2,1,1)
plot(-depth(1:lakebed), V(3,:));
xlabel('depth')
ylabel('spatial weight')
title('EOF 3')
axis([0 3 -.2 .2])

subplot(2,1,2)
plot(time_avg, U(:,3))
xlabel('time')
datetick('x',13,'keepticks')
dynamicDateTicks()
ylabel('temporal weight')
title('PC 3')
saveas(gcf,[dir '/tau_pattern3.jpg']);

%close all
%end
