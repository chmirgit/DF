clear all
close all

load('speakerA.mat');
load('speakerB.mat');

n = length(u); %number of time steps
M = 6600; %order of the filter >= 4

%% LMS
tic
mulms = 0.001;
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


timerls = toc;