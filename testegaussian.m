clear
close all
dt = 1e-8;
t = 1*dt:dt:500*dt;

Ap = 1;
tau = 10*dt;
t0 = 200*dt;


p = -Ap*sqrt(2*exp(1)/tau^2)*(t-t0).*exp(-((t-t0)/tau).^2);

dp = zeros(1,500);

for k = 2:499
    dp(k) = (p(k+1)-p(k-1))/(t(k+1)-t(k-1));
end

plot(t,p)

figure
plot(t,abs(dp))