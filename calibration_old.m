% Nortek Signature1000
% Average 4096 profiles into one array of 107 bins
% sampling rate supposed to be at 8Hz, 1MHz transducer
clc; clear; %close all
%% Sonar Equation log transform (Hoitink)
% What the adcp records
% Sv = 2*alpha*R + Kc*(E-Er) + 10*log10((Tt*R^2)/(L*Pt)) + C;
% Sediment backout:
% Sv = 10*log10((3*eta*k^4)/rho_s*a^3*M) + [2*alpha*R] + C
%       M = sediment mass concentration (kg/m^3)
%       a = average particle radius (m)
%       eta = ((e-1)/(3*e))^2 + 1/3*((sigma-1)/(2*sigma+1))^2
%       e = ratio of elasticity of sediment to water
% Sv = volume backscatter strength (dB)
% R = range along central beam axis (m), adcp calculated
% alpha = attenuation coeff = alpha_w + alpha_s
%       alpha_s = ((k^4*a^3)/(96*rho_s) +
%         (k*(sigma-1)^2*s)/(s*rho_s*(s^2 + (sigma + delta)^2)))*20/ln(10)*M
%       k = wave number (m^-1)
%       s = 9/(2*beta*a)*(1 + 2/(beta*a))
%       sigma = rho_s/rho_w
%       rho = density of water/sediment (kg/m^3)
%       delta = .5*(1 + 9/(beta*a))
%       beta = (pi*f/v)^.5
%       f = emitted sigal (MHz)
%       v = kinematic viscosity of water (m^2/s)
% E = echo intensity (counts), adcp spec
% Er = recieved noise (counts), beam spec
% Kc = scale factor (dB/count), beam spec
% Tt = transducer temp (C), adcp spec
% L = transmit pulse length (m) - translate from (s)
% Pt = transmit power (W), adcp spec
% C = constant (dB)

% E = 100; % counts - VTamp (4096 in a burst?)
% Er = 55; % counts - Hnoise
% Kc = 1/.45; % dB/counts
% Tt = 17; % C
% Pt = 5; % W
% L = .04; % m, for a 4cm bin size
z_bottom = 72;

matfile = 'S100882A003_LakeK_Sig2_13.mat';
[Amp, c, T, r, f] = load_LakeK_Sig_file(matfile);

a = 100e-6; % particle size, m
%f = 1; % MHz
%c = 1485; % m/s
%T = 17; % C
k = 2*pi*f*1e6./c; % wave number, /m
v = 1e-6; % viscosity of water, m^2/s
rho_s = 1100; % kg/m^3, clay/silt/organics?
rho_w = 1000; % kg/m^3

beta = ((pi*f*1e6)/v)^.5;
delta = .5*(1 + 9/(beta*a));
s = 9/(2*beta*a)*(1 + 2/(beta*a));
sigma = rho_s/rho_w;
e = 1;
eta = ((e-1)/(3*e))^2 + 1/3*((sigma-1)/(2*sigma+1))^2;

% water and sediment/particle attenuation
alpha_w = @(T) (55.9 - 2.37.*T + 4.77e-2.*T.^2 - 3.84e-4.*T.^3).*10^-3*f^2;
alpha_s = @(k,M) ((k.^4.*a.^3)./(96.*rho_s) + (k.*(sigma-1)^2.*s)./(s*rho_s...
    *(s^2 + (sigma + delta)^2))).*(20./log(10))*M;
alpha = @(T,k,M) alpha_w(T) + alpha_s(k,M);

% Cannot truly solve unless C and a and rho_s are known.
C = 90;
M_calc = zeros(length(k),length(r));
A_corrected = zeros(length(k),length(r));
for i=1:length(k) % rows of Sv
    for j=1:72%length(r) % columns
        M = .01;
        t = true;
        while t==true
            % Hoitink's claimed concentration from corrected amplitude relation
            % Sv = 10*log10((3*eta*k(i)^4)/rho_s*a^3*M)
            % M is prop to 10^(Sv/25), i.e. Sv comes out negative
            if k(i)*a < 1 % Rayleigh scattering
                Sv = -10*log10((3*eta*k(i)^4*a^3)/rho_s*M);
                % C depends completely on a and rho_s
                % M = linspace(.01,1000,10000);
                % plot(M,10*log10((3*eta*k(i)^4*a^3)/rho_s*M)+C)
                % hold on; plot(M, Amp(i,j)*.43 + 20*log10(r(j)) + 2*alpha(T(i),k(i),M)*r(j))
                % assuming its the first intersection
            else % ka>1
            	Sv = -10*log10((3*eta*k(i)^4*a^3)/(2*rho_s*(2+3*(k(i)*a)^4))*M);
            end
            % Nortek's stated sonar equation to calculate corrected amplitude
            % Sv = Amp(i,j)*.43 + 20*log10(r(j)) + 2*alpha(T(i),k(i),M)*r(j)
            Sv_actual = Amp(i,j)*.43 + 20*log10(r(j)) + 2*alpha(T(i),k(i),M)*r(j);
            
            if abs(1-Sv/Sv_actual) < .001
                M_calc(i,j) = 1/M; % for some reason, with the negatives above, M becomes inverted
                A_corrected(i,j) = Sv;
                t = false;
            else
                fprintf('M = %d; ', M)
                fprintf('Sv = %4.2f; ', Sv)
                fprintf('Sv_actual = %4.2f\n', Sv_actual)
                M = M + .01;
            end
        end
    end
end
            
figure()
surface(flipud(M_calc(:,1:z_bottom)'),'EdgeColor','None');
colorbar
title('Calc Concentration')
xlabel('Depth Bins');
ylabel('Time (hr)');

figure()
surface(flipud(A_corrected(:,1:z_bottom)'),'EdgeColor','None');
colorbar
title('Corrected Amplitude (dB)')
xlabel('Depth Bins');
ylabel('Time (hr)');

%% Attenuation Equation (Thorne, 1991)
% % Sv = k0*((tau*c*M)^.5/r*psi)*exp(-2*r*(alpha_w+alpha_s))
% % Sv = ensemble averaged rms signal (dB)
% % k0 = constant
% % t = pulse length (s)
% % c = speed of sound in water (m/s)
% % M = partical mass concentration (kg/m^3)
% % r = range (m)
% % psi = spreading loss
% %     = 1 for farfield
% %     = (1 + 1.35*z + 2.5*z)^3.2/(1.35*z + 2.5*z)^3.2 for nearfield
% %       z = r/r_crit, r_crit = pi*d_transducer^2*k
% % alpha_w = attenuation coeff of water (T = [C], f = [MHz])
% %         = (55.9 - 2.37*T + 4.77e-2*T^2 - 3.84e-4*T^3)*10^-3*f^2
% % alpha_s = attenuation coeff of sediment
% %         = (1/r)*integral from 0 to r of zeta*M*dr - depth averaged
% %       zeta = (1/rho_w*a)*(beta*(k*a)^4/(1 + (k*a)^2 + (4/3)*beta*(k*a)^4))
% %       a = particle equivalent sphere radius
% %       beta = (yk^2 + yp^2/3)/6
% %       yk = (k_s - k_w)/k_w
% %       yp = 3*(rho_s - rho_w)/(2*rho_s + rho_w)
% %       k = compressibility of water and sediment(s)
% %       rho = density of water and sediment(s)
% 
% matfile = 'S100882A003_LakeK_Sig2_8.mat';
% [Sv_actual, c, T, r, f] = load_LakeK_Sig_file(matfile);
% z_bottom = 71;
% 
% % For winter calibration
% % Sv_actual = []; c = []; T = []; r = []; f = [];
% % for i=141:145
% %     matfile = ['S100882A004_LakeK3_' num2str(i) '.mat'];
% %     [Sv, C, t, r, f] = load_LakeK_Sig_file(matfile);
% %     Sv_actual = [Sv_actual; Sv];
% %     c = [c; C];
% %     T = [T; t];
% % end
% 
% a = 15e-6; % effective particle scattering diameter
% k = 2*pi*f*1e6./c; % wave number, /m
% tau = .04./c; % pulse length, s
% rho_s = 1100; % particle density
% rho_w = 1000; % water density
% bin_size = .04; % m
% blanking_dis = .1; % m
% d_transducer = .04; %m, tranducer diameter
% %r = (1:107)*bin_size + blanking_dis; % array of bin ranges
% r_crit = pi*d_transducer^2*k; % critical beam distance
% z = r./r_crit;
% 
% % Beam attentuation
% psi = zeros(size(z));
% for i=1:length(r_crit)
%     for j=1:length(r)
%         if z(i,j)<1
%             psi(i,j) = (1 + 1.35.*z(i,j) + 2.5.*z(i,j)).^3.2./(1.35.*z(i,j) + 2.5.*z(i,j)).^3.2;
%         else
%             psi(i,j) = 2;
%         end
%     end
% end
% 
% % Water attenuation
% alpha_w = @(T) (55.9 - 2.37.*T + 4.77e-2.*T.^2 - 3.84e-4.*T.^3).*10^-3*f^2;
% 
% % Sediment and water attenuation
% k_s = 5e-10; % Pa
% k_w = 5e-10; % Pa
% yp = 3*(rho_s - rho_w)/(2*rho_s + rho_w);
% yk = (k_s - k_w)/k_w;
% beta = (yk^2 + yp^2/3)/6;
% if k(1)*a < 1 %Rayleigh scattering
%     zeta = beta/(rho_w.*a)*(k.*a).^4;
% else
%     zeta = (1./rho_w.*a).*(beta.*(k.*a).^4./(1 + (k.*a).^2 + (4/3).*beta.*(k.*a).^4));
% end
% alpha_s = @(zeta,SSC) zeta*SSC; %(1/r(end)).*trapz(zeta.*SSC, 2);
% alpha = @(T,zeta,SSC) alpha_w(T) + alpha_s(zeta,SSC);
% 
% k0 = 8;
% %SSC = ((Sv_actual./(k0.*exp(-2.*r.*alpha))).*(r./psi)).^2./(tau.*c);
% SSC_calc = zeros(length(k),length(r));
% for i=1:length(k) % rows of Sv
%     for j=45:72 %length(r) % columns
%         SSC = .0001;
%         t = true;
%         while t==true
%             Sv = k0*((tau(i)*c(i)*SSC)^.5/r(j)*psi(i,j))*exp(-2*r(j)*alpha(T(i),zeta(i),SSC));
%             
%             if (1-Sv/Sv_actual(i,j) < .001)
%                 SSC_calc(i,j) = SSC;
%                 t = false;
%             else
%                 fprintf('SSC = %d; ', SSC)
%                 fprintf('Sv = %4.2f\n', Sv)
%                 SSC = SSC + .0001;
%             end
%         end
%     end
% end
% 
% figure()
% surface(SSC(:,1:z_bottom),'EdgeColor','None')
% 
% %% Not so sure this is accurate considering the papers claim a M = 10^(Sv/25) relationship
% % Removes exponential function - counteracts what papers say about exp relation
% if strcmp(matfile, 'S100882A004_LakeK3_145.mat') == 0
%     winter_calibration = load(['winter_calibration' num2str(a) '.txt']);
%     SSC2 = SSC2 - winter_calibration;
% end
% 
% figure()
% surface(SSC2(:,1:z_bottom),'EdgeColor','None')
% 
% % Winter exponential zero calibration
% if strcmp(matfile, 'S100882A004_LakeK3_145.mat') == 1
%     winter_cal = mean(SSC2);
%     dlmwrite(['winter_calibration' num2str(a) '.txt'], winter_cal);
% end

%%
% sediment ranges (10-500um, 1010-1330kg/m3) - ip

% differences between attentuation types (water, sediment, beam)
% lit review lake scatterers (plankton, sediment, density microstructure) - ip
% timeseries, check with adv
