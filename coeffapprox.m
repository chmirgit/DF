clear all
close all

n = 20000; %number of time steps
M = 4; %order of the filter >= 4

%generating white noises and signals
sigma_v1 = 0.42;
sigma_v2 = 0.72;

v1 = sqrt(sigma_v1)*randn(1,n);
v1 = v1 - mean(v1);

v2 = sqrt(sigma_v2)*randn(1,n);
v2 = v2 - mean(v2);

%u(n) = ucoeff * u(n:n-2) + v1(n)
u = zeros(n,1);
ucoeff = [-0.87 -0.22 -0.032];
u(1) = v1(1);
u(2) = ucoeff(1)*u(1) + v1(2);
u(3) = ucoeff(1)*u(2) + ucoeff(2)*u(1) + v1(2);
for i=4:n
    u(i) = ucoeff * u(i:-1:i-2) + v1(i);
end

%s(n) = scoeff * u(n:n-3);
s = zeros(n,1);
scoeff = [-0.13 0.67 -0.18 0.39];
s(1) = scoeff(1)*u(1);
s(2) = scoeff(1)*u(2) + scoeff(2)*u(1);
s(3) = scoeff(1)*u(3) + scoeff(2)*u(2) + scoeff(3)*u(1);
for i=4:n
    s(i) = scoeff * u(i:-1:i-3);
end


%x(n) = xcoeff * x(n-1:n-3) + v2(n)
x = zeros(n,1);
xcoeff = [-0.57 -0.16 -0.08];
x(1) = v2(1);
x(2) = xcoeff(1)*x(1) + v2(2);
x(3) = xcoeff(1)*x(2) + xcoeff(2)*x(1) + v2(2);
for i=4:n
    x(i) = xcoeff * x(i:-1:i-2) + v2(i);
end

%desired signal
d = x + s;




%autocorr vector ru
%matrix is derived from given signal descriptions
sm = [1     0.87  0.22  0.032;...
      0.87  1.22  0.032 0    ;...
      0.22  0.902 1     0    ;...
      0.032 0.22  0.87  1    ];

ru = sm \ [sigma_v1; 0; 0; 0];

%autocorr matrix for filter of M coefficients
p = zeros(M,1);

if M>4
    for k=5:M
        ru = [ru ; ucoeff * ru(k-1:-1:k-3)];
    end
end

%calculating crosscorrelation
for k=1:M
    idx = k:-1:k-3;
    for l=1:length(idx)
        if idx(l)<=0
           idx(l) = 2-idx(l);
        end
    end
    p(k) = scoeff * ru(idx);
end
Ru = toeplitz(ru);
wo = Ru \ p; %optimal wiener coefficients

%% LMS
tic
mulms = 0.0001;
wlms = zeros(M,1);
ylms = zeros(n,1);
elms = zeros(n,1);
Jlms = zeros(n,1);

for i=(M+1):n
    ylms(i) = wlms' * u(i:-1:(i-M+1));
    elms(i) = d(i) - ylms(i);
    wlms = wlms + mulms * elms(i) * u(i:-1:(i-M+1));
    
    Jlms(i) = elms(i)^2;
end
figure(1)
semilogy(Jlms(M+1:end));
xlabel('time step n');
ylabel('e^2(n)');
legend({'mu = 0.0001'});
title('LMS');
timelms = toc;

%% N-LMS
tic
munlms = 0.0001;
alpha = 0.1;
wnlms = zeros(M,1);
ynlms = zeros(n,1);
enlms = zeros(n,1);
Jnlms = zeros(n,1);

for i=(M+1):n
    U = u(i:-1:(i-M+1));
    ynlms(i) = wnlms' * U;
    enlms(i) = d(i) - ynlms(i);
    dinom = alpha + (U' * U);
    wnlms = wnlms + (munlms/dinom)*U*enlms(i);
    Jnlms(i) = enlms(i)^2;
end

figure(2)
semilogy(Jnlms(M+1:end));
xlabel('time step n');
ylabel('e^2(n)');
legend({'mu = 0.0001 | alpha = 0.1'});
title('N-LMS');
timenlms = toc;

%% RLS
tic
delta = 1/1000;
lambda = 1;
wrls = zeros(M,1);
yrls = zeros(n,1);
erls = zeros(n,1);
Jrls = zeros(n,1);

P = (1/delta) * eye(M,M);
for i=(M+1):n
    yrls(i) = wrls' * u(i:-1:(i-M+1));
    k = ( (lambda^-1)*P*u(i:-1:i-M+1) / (1+(lambda^-1)*u(i:-1:i-M+1)'*P*u(i:-1:i-M+1)));
    erls(i) = d(i) - yrls(i);
    wrls = wrls + k * erls(i);
    P = (lambda^-1)*P - (lambda^-1)*k*u(i:-1:i-M+1)'*P;
    
    Jrls(i) = erls(i)^2;
end

figure(3)
semilogy(Jrls(M+1:end));
xlabel('time step n');
ylabel('e^2(n)');
legend('delta = 1/1000');
title('RLS');
timerls = toc;
