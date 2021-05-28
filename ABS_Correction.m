function [M_calc, A_corrected, dir] = ABS_Correction(matfile, a, rho_s, saveon)

z_bottom = 72;
if saveon==true
    dir = ['./calibrationFiles' matfile(1:end-4) '_' num2str(a) 'm' num2str(rho_s) 'density'];
    mkdir(dir);
else
    dir = 'None';
end

%matfile = 'S100882A003_LakeK_Sig2_8.mat';
[Amp, c, T, r, f] = load_LakeK_Sig_file(matfile);

%a = 500e-6; % particle size, m
%f = 1; % MHz
%c = 1485; % m/s
%T = 17; % C
k = 2*pi*f*1e6./c; % wave number, /m
v = 1e-6; % viscosity of water, m^2/s
%rho_s = 1100; % kg/m^3, clay/silt/organics?
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

% Engineered concentration Sv so that it can only intersect the sonar Sv at
% one point, otherwise it's not possible to solve these two equations as
% they rely heavily on C, a, and rho_s
M_calc = zeros(length(k),z_bottom);
A_corrected = zeros(length(k),z_bottom);
for i=1:length(k) % rows of Sv
    for j=1:z_bottom %length(r) % columns
        %M = .01;
        if rho_s<1050
            M = .001;
        elseif rho_s<1100
            M = .0008;
        elseif rho_s<1200
            M = .0005;
        else
            M = .0001;
        end
        t = true;
        while t==true
            % Hoitink's claimed concentration from corrected amplitude relation
            % M is prop to 10^(Sv/25)
            if k(i)*a < 1 % Rayleigh scattering
                Sv = -10*log10((3*eta*k(i)^4)/rho_s*a^3*M);
            else % ka>1
            	Sv = -10*log10((3*eta*k(i)^4*a^3)/(2*rho_s*(2+3*(k(i)*a)^4))*M);
            end
            % Nortek's stated sonar equation to calculate corrected amplitude
            Sv_actual = Amp(i,j)*.43 + 20*log10(r(j)) + 2*alpha(T(i),k(i),M)*r(j);
            
            if abs(1-Sv/Sv_actual) < .001
                M_calc(i,j) = 1/M; % Univerting M from Hoitink Sv inversion
                A_corrected(i,j) = Sv;
                t = false;
            else
                if Sv_actual > 100
                    error(['M step too high ' num2str(a) ' ' num2str(rho_s)])
                end
%                 fprintf('M = %d; ', M)
%                 fprintf('Sv = %4.2f; ', Sv)
%                 fprintf('Sv_actual = %4.2f\n', Sv_actual)
                %M = M + .01;
                if rho_s<1050
                    M = M + .001;
                elseif rho_s<1100
                    M = M + .0008;
                elseif rho_s<1200
                    M = M + .0005;
                else
                    M = M + .0001;
                end
            end
        end
    end
    disp(i)
end

if saveon==true
    figure()
    surface(flipud(M_calc(:,1:z_bottom)'),'EdgeColor','None');
    colorbar
    title('Calc Concentration')
    xlabel('Depth Bins');
    ylabel('Time (hr)');
    saveas(gcf, [dir '/Calc_Conc.jpg']);

    figure()
    surface(flipud(A_corrected(:,1:z_bottom)'),'EdgeColor','None');
    colorbar
    title('Corrected Amplitude (dB)')
    xlabel('Depth Bins');
    ylabel('Time (hr)');
    saveas(gcf, [dir '/Amp_Cor.jpg']);
end

end