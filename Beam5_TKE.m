function [TKE, Time_hourly] = Beam5_TKE(matfile)

load(matfile);

bins = 72;
w = Data.BurstHR_VelBeam5(:,1:bins);
Time = Data.BurstHR_Time;

FlowTime = datetime(Time,'ConvertFrom','datenum');
FlowTT = timetable(FlowTime,w(:,1)); % pick random bin for timetable
Flow_Hourly = retime(FlowTT,'hourly','mean');
Time_hourly = datenum(Flow_Hourly.FlowTime);

% Make z-scale for HR bins
z = double(Config.BurstHR_BlankingDistance+Config.BurstHR_CellSize.*(0:bins));
z_midbin = zeros(length(z)-1,1);
for j = 1:(length(z)-1)
    z_midbin(j)=(z(j)+z(j+1))/2;
end

% Despike data
w_despiked = w;
spikestd = 1;  % set higher for weaker despiking
for i = 1:bins
    % Pull out column
    w_temp=w_despiked(:,i);
    wspike = find( abs(w_temp) > (spikestd*nanstd(w_temp) + abs(nanmean(w_temp))) );
    w_temp(wspike) = NaN;
    % replace bad points with burst means
    w_temp(isnan(w_temp)) = nanmean(w_temp);
    
    % Put despiked data back into column
    w_despiked(:,i)=w_temp;
end
    
% Average over some period, find fluctuating velocity (w') and TKE
n = 4096; % Samples we're averaging over (samples per burst)
% max(Data.BurstHR_EnsembleCount)

% Calculate TKE
for i = 1:bins
    % Make w vectors divisble by number of samples we're averaging 
    ensembles = floor(length(w_despiked(:,i))/n); % ensembles = number of sampling intervals (hrs)
    w_temp2(:,i)=w_despiked(1:(n*ensembles),i); % cut # of samples so divisible by ensembles
    
    % Average trimmed w
    w_reshaped = reshape(w_temp2(:,i),n,ensembles);
    w_bar = mean(w_reshaped);

    for j = 1:ensembles
        w_prime(:,j) = w_reshaped(:,j)-w_bar(j);
    end
    
    % power spectra
    if i==63
        for j=1:ensembles
            [psd_w(:,j), f(:,j)] = pwelch(w_prime(:,j),[],[],[],16);
        end
        figure
        loglog(f, psd_w);
        xlabel('Frequency (Hz)');
        ylabel('PSD');
        title('w dir, bin 63');
    end
    
    % Calculate TKE
    w_prime_squared = w_prime.^2;
    TKE(:,i) = mean(w_prime_squared,1);
end
pwelch(w_prime,[],[],[],8);

Time_hourly = Time_hourly(1:ensembles);

end
