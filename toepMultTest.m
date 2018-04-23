% Plots times for multiplications and prints mean 
% error throughout the trials

clear ;
close all;

n = 8000;
u = randn(n,1);

stdmult = 0;
toepmult = 0;
err = 0;

for m = 500:500:8000
    w = randn(m,1);
    T = toeplitz(u(m:n), u(m:-1:1));
    
    % standard matrix multiplixation
    tic
    y = T*w;
    stdmult = [stdmult, toc];
    
    % toepMult
    tic
    yf = toepMult(u,w);
    toepmult = [toepmult, toc];
    
    err = [err, norm(yf-y)/norm(y)];
end

s1 = subplot(2,1,1);
plot(stdmult);
title(s1, 'standard');
ylabel(s1,'time');
xlabel(s1,'*500 w elements');

s2 = subplot(2,1,2);
plot(toepmult);
title(s2, 'toeplitz');
ylabel(s2,'time');
xlabel(s2,'*500 w elements');

meanErr = mean(err(2:end))