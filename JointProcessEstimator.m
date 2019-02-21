clear all
close all

load('music.mat');

%s = signal

M = 100; %order of the filer
sigma = var(s);
delta = 100; %lags

%autocorrelation vector
[r,delta] = autocorr(s,delta);

%unnormalized
unr = r*sigma;

R = toeplitz(unr,unr');

%find wo
%here we rearrange the equations to solve for the vector:
%wap = [Pm , am1 ... amm] = [Pm , -w1 -w2 ...]
Rm = R;
Rm(:,1) = [-1 ; zeros(M,1)];
eq = -unr;
wap = Rm \ eq;
wo = -wap(2:end);

%reflection coefficients
[ai,g,L,D] = LevinsonDurbin_iterative(M,r);
woi = -ai(2:end);

%filter parameters compared to L-D iterative
dev = norm(wo-woi);

tic
[f, b] = LatticeFilter(s, g);
latticeFilterTime = toc;

%find P using the vector am = [1 -w1 -w2 ...]
P = R * [1 ; -wo];
gamma = L * P ./ D;

tic
y = jointProcessFilter(s,g,gamma);
JPFilterTime = toc;

J = (y-s) .^ 2;
semilogy(J)
title('error')
xlabel('time');