matfiles = {'S100882A003_LakeK_Sig2_8' 'S100882A003_LakeK_Sig2_31'...
    'S100882A003_LakeK_Sig2_26'};

a = [15e-6 65e-6 100e-6 250e-6 500e-6 750e-6]; % particle diameters
rho_s = [1010 1050 1100 1200 1330]; % particle desities
% .75mm at 1330 looks to be the winner blue
% or .065mm at 1010 salmon red
% of .01mm at 1010 purple

%% Force solution for concentrations
% for k=6:length(a)
%     for j=1:length(rho_s)
%         % attenuations: using .1mm diameter and 1100kg/m3 for now
%         for i=1:length(matfiles)
%             matfile = [matfiles{i} '.mat'];
%             [conc, amp_cor, dir] = ABS_Correction(matfile, a(k), rho_s(j), true);
%             save([dir '\' matfiles{i} '_concentration.mat'], 'conc', 'amp_cor');
%         end
%     end
%     close all
% end

% calculate tke for each adcp datafile
for i=1:length(matfiles)
    matfile = [matfiles{i} '.mat'];
    [tke, time] = Beam5_TKE(matfile);
    save([matfiles{i} '_tke.mat'],'tke','time')
end

% normalize concentration data - set bottom = 1
for k=1:length(a)
    for j=1:length(rho_s)
        for i=1:length(matfiles)
            dir = ['.\calibrationFiles' matfiles{i} '_' num2str(a(k)) 'm' num2str(rho_s(j)) 'density'];
            %addpath(fullfile(['E:\Documents\UW\CEWA 600\' matfiles{i} '_' num2str(a(k)) 'm' num2str(rho_s(j)) 'density']));
            load([dir '\' matfiles{i} '_concentration.mat']);
            %conc_orig = load([matfiles{i} '_concentration.mat']);
            normal_factor = max(conc,[],2); % set bottom conc = 1
            normal_conc = conc./normal_factor; % normalize by bottom conc
            avg_normal = mean(normal_conc,1);
            save([dir '\' matfiles{i} '_concentration.mat'], 'conc', 'amp_cor', 'normal_conc', 'avg_normal');
            
            figure()
            surface(flipud(normal_conc'),'EdgeColor','None');
            colorbar
            title('Normalized Concentration')
            xlabel('Depth Bins');
            ylabel('Time (hr)');
            saveas(gcf, [dir '/Normalized_Conc.jpg']);
        end
    close all
    end
end

%% apply linear correction to normalized concentration and corrected
% applitude - it appears the beam spreading attenuation is too strong
for k=1:length(a)
    for j=1:length(rho_s)
        for i=1:length(matfiles)
            dir = ['.\calibrationFiles' matfiles{i} '_' num2str(a(k)) 'm' num2str(rho_s(j)) 'density'];
            load([dir '\' matfiles{i} '_concentration.mat']);
            
            background_conc = max(mink(normal_conc,2,1)); % second smallest in case of bad data?
            %plot(background_conc');
%             x = find(diff(background_conc)./max(background_conc)>.3);
%             if numel(x)>1
%                 x = min(x);
%             end
            x = 69;
            % setting background suspension = 0, i.e unknown C constant
            final_cal_conc = normal_conc - [background_conc(1:x) zeros(1,length(background_conc)-x)];
            final_avg_conc = mean(final_cal_conc,1);
                        
            % applying linear fix to corrected amplitude
            load([dir '\' matfiles{i} '_concentration.mat']);
            background_dB = max(mink(amp_cor,2,1));
            final_cor_amp = amp_cor - [background_dB(1:x) zeros(1,length(background_conc)-x)];
            
            save([dir '/' matfiles{i} '_calibrated_data.mat'], 'final_cal_conc', 'final_avg_conc', 'final_cor_amp');
            
            figure()
            surface(flipud(final_cal_conc'),'EdgeColor','None');
            colorbar
            title('Final Normalized Concentration')
            xlabel('Depth Bins');
            ylabel('Time (hr)');
            saveas(gcf, [dir '/Final_Cal_Conc.jpg']);
            
            figure()
            surface(flipud(final_cor_amp'),'EdgeColor','None');
            colorbar
            title('Final Corrected Amplitude (-C dB)')
            xlabel('Depth Bins');
            ylabel('Time (hr)');
            saveas(gcf, [dir '/Final_Cor_Amp.jpg']);
        end
        close all
    end
end

%% measuring uncertainty in particle size? -
% take means of each concentration dataset and coplot them

for i=1:length(matfiles)
    figure()
    for k=1:length(a)
        for j=1:length(rho_s)
            dir = ['.\calibrationFiles\' matfiles{i} '_' num2str(a(k)) 'm' num2str(rho_s(j)) 'density'];
            load([dir '\' matfiles{i} '_concentration.mat']);
            
            plot(1:length(avg_normal), avg_normal)
            hold on
        end
    end
    xlabel('Depth Bins (4cm each)');
    ylabel('Normalized Concentration');
end

for i=1:length(matfiles)
    figure()
    for k=1:length(a)
        for j=1:length(rho_s)
            dir = ['.\calibrationFiles\' matfiles{i} '_' num2str(a(k)) 'm' num2str(rho_s(j)) 'density'];
            load([dir '\' matfiles{i} '_calibrated_data.mat']);
            
            plot(1:length(final_avg_conc), final_avg_conc, 'DisplayName', [num2str(a(k)) 'um ' num2str(rho_s(j)) 'kg/m^3'])
            hold on
        end
    end
    title('Final Calibrated Concentration');
    xlabel('Depth Bins (4cm each)');
    ylabel('Normalized Concentration');
    axis([0 length(final_avg_conc) 0 .1]);
    %legend('Location','Northwest');
    saveas(gcf, [matfiles{i} '.jpg']);
end


%% Create individual timeseries of corrected amplitude (dB), normalized concentration, and TKE
% Matlab can't find directories it just created. That's dumb.
for k=1:length(a)
    for j=1:length(rho_s)
        for i=1:length(matfiles)
            dir = ['.\calibrationFiles' matfiles{i} '_' num2str(a(k)) 'm' num2str(rho_s(j)) 'density'];
            load([dir '\' matfiles{i} '_calibrated_data.mat']);
            load([matfiles{i} '_tke.mat']);

            figure()
            bin = [60 64 68]; % depth bin to plot
            subplot(3,1,1)
            plot(time, final_cor_amp(:,bin));
            xlabel('time, hr');
            datetick('x');
            ylabel('corrected amplitude, -C dB');
            title(matfiles{i});
            legend('Bin 60','Bin 64','Bin 68');

            subplot(3,1,2)
            plot(time, final_cal_conc(:,bin));
            xlabel('time, hr');
            datetick('x');
            ylabel("concentration, normalized");
            legend('Bin 60','Bin 64','Bin 68');

            subplot(3,1,3)
            plot(time, tke(:,bin));
            xlabel('time, hr');
            datetick('x');
            ylabel('TKE, m^2/s^2');
            legend('Bin 60','Bin 64','Bin 68');

            saveas(gcf, [dir '/TS.jpg']);
        end
        close all
    end
end
