function [Amp, c, Temp, r, f] = load_LakeK_Sig_file(matfile)
% reads in necessary data from Signature1000 datafiles and averages by
% burst

load(matfile);

dB_received = Data.BurstHR_AmpBeam5;
sound = Data.BurstHR_Soundspeed;
temp = Data.BurstHR_Temperature;

% Average data by hour bursts
%sample_freq = Config.Burst_SamplingRate; % in Hz
n = Config.Burst_NSample;
hours = length(dB_received)/n; % number of hours/number of bursts
r = (1:double(Data.BurstHR_NCells(1)))*Config.BurstHR_CellSize + ...
    Config.Burst_BlankingDistance - Config.BurstHR_CellSize/2; % dB range, in m
f = Config.HeadFrequency/1000; % in MHz

Amp = zeros(hours,Data.BurstHR_NCells(1));
c = zeros(hours,1);
Temp = zeros(hours,1);
for i=0:hours-1
    Amp(i+1,:) = mean(dB_received(i*n+1:i*n+n,:),1);
    c(i+1,:) = mean(sound(i*n+1:i*n+n,:),1);
    Temp(i+1,:) = mean(temp(i*n+1:i*n+n,:),1);
end
end
    
    



