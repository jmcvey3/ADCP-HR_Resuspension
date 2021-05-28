a = 750e-6;
rho_s = 1330;

% range = linspace(1,39,39);
% C = []; dB = []; TKE = []; t = [];
% for i=range
%      matfile = ['./adcpDataFiles/S100882A003_LakeK_Sig2_' num2str(i) '.mat'];
%      [conc, amp_cor, dir] = ABS_Correction(matfile, a, rho_s, false);
%      [tke, time] = Beam5_TKE(matfile);
%      
%      C = [C; conc];
%      dB = [dB; amp_cor];
%      TKE = [TKE; tke];
%      t = [t; time];
%      %fprintf('%4.2f; %4.2f\n', length(tke(:,1)), length(time));
% end

% two extra ensembles from sig 6 and 7 - not full hour data

% normal_factor = max(C,[],2); % set bottom conc = 1
% normal_conc = C./normal_factor; % normalize by bottom conc
% 
% background_conc = max(mink(normal_conc,2,1));
% %x = find(diff(background_conc)./max(background_conc)>.3);
% %if numel(x)>1
% %   x = min(x);
% %end
% x = 69;
% % setting background suspension = 0, i.e unknown C constant
% final_cal_conc = normal_conc - [background_conc(1:x) zeros(1,length(background_conc)-x)];
% 
% % applying linear fix to corrected amplitude
% background_dB = max(mink(dB,2,1));
% final_cor_amp = dB - [background_dB(1:x) zeros(1,length(background_conc)-x)];
% 
% save('./adcpDataFiles/S100882A004_LakeK_Sig2.mat', 'C', 'final_cal_conc', 'dB', 'final_cor_amp', 'TKE', 't');

%% Aug18-Nov 18 Surface plots of corrected amplitude and concentration
load('./adcpDataFiles/S100882A004_LakeK_Sig2.mat', 'C', 'final_cal_conc', 'dB', 'final_cor_amp', 'TKE', 't');

figure()
surface(flipud(final_cal_conc'),'EdgeColor','None');
colorbar
title('Final Normalized Concentration')
ylabel('Depth Bins');
xlabel('Time (hr)');

figure()
surface(flipud(final_cor_amp'),'EdgeColor','None');
colorbar
title('Final Corrected Amplitude (-C dB)')
ylabel('Depth Bins');
xlabel('Time (hr)');

figure()
surface(flipud(TKE'),'EdgeColor','None');
colorbar; caxis([0 .3e-4])
title('Vertical TKE (m^2/s^2)')
ylabel('Depth Bins');
xlabel('Time (hr)');

% Surface plots timeseries
figure()
bin = [60 64 68]; % depth bin to plot
subplot(3,1,1)
plot(t, final_cor_amp(:,bin));
xlabel('time, hr');
datetick('x');
ylabel('corrected amplitude, -C dB');
title('August-November');
legend('Bin 60','Bin 64','Bin 68');

subplot(3,1,2)
plot(t, final_cal_conc(:,bin));
xlabel('time, hr');
datetick('x');
ylabel("concentration, normalized");
legend('Bin 60','Bin 64','Bin 68');

subplot(3,1,3)
plot(t, TKE(:,bin));
xlabel('time');
datetick('x');
ylabel('TKE, m^2/s^2');
legend('Bin 60','Bin 64','Bin 68');

