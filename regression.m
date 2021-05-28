function y_reg = regression(x, y, t_bound, Time)
% Linear Regression
% y = a + bx
% b = Sxy/Sxx
% a = ybar - b*xbar
% Sxy = sum(x_i*y_i) - (1/n)*(sum(x_i)*sum(y_i))
% Sxx = sum(x_i**2) - (1/n)*(sum(x_i)**2)
% ybar = (1/n)*sum(y_i)
% xbar = (1/n)*sum(x_i)
n = length(x);
xbar = nanmean(x);
ybar = nanmean(y);
Sxy = nansum(x.*y) - (1./n).*(nansum(x)*nansum(y));
Sxx = nansum(x.^2) - (1./n).*(nansum(x).^2);
Syy = nansum(y.^2) - (1./n).*(nansum(y).^2);
b = Sxy/Sxx;
a = ybar - b*xbar;

y_reg = a + b.*x;
plot(x, y_reg, 'g');
R_2 = Sxy^2/(Sxx*Syy);
fprintf('R: %4.4f\n', R_2^.5)
fprintf('R^2: %4.4f\n\n', R_2);
%fprintf('Only 12 percent of the variance is accounted for?\n\n');

% b confidence bounds
%t_bound = 1.96; %stats.t.ppf(1-.05/2,n-2)

% uncertainty
y_err = y - y_reg;
sse = nansum(y_err.^2);
st_err = sqrt(sse/(n-2));

s_b = sqrt(st_err^2/(nansum((x-xbar).^2)));
b_lower = b - t_bound*s_b;
b_upper = b + t_bound*s_b;
t_test = b/s_b;
fprintf('t_test: %4.2f\n', t_test);
fprintf('t_stat to beat: %4.2f\n\n', t_bound);
%fprintf('slope stat is still ridiculously high\n\n');

fprintf('lower 95th slope: %4.3f\n', b_lower)
fprintf('regression slope: %4.3f\n', b)
fprintf('upper 95th slope: %4.3f\n\n', b_upper)

p_x = linspace(nanmin(x), nanmax(x), n);
sigma_Ep = sqrt(st_err^2.*(1 + 1/n + (n.*(p_x-xbar).^2)./(n*nansum(x.^2)-(nansum(x.^2).^2))));
p_y = a + b.*p_x;
y_lower = p_y - t_bound*sigma_Ep;
y_upper = p_y + t_bound*sigma_Ep;
plot(p_x, y_lower,'r', p_x, y_upper, 'r')
legend('data','linear regression','95% bounds');

% figure()
% plot(Time, y, Time, y_reg)
% xlabel('Time')
% ylabel('Conc from Linear Regression')
% legend('data','linear regression');
end