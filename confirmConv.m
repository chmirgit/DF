% Test of convolution theorem
% computing the convolutio of 2 random vectors  in 4 ways

clear ;
close all;

n = 4; %length of x
m = 3; %length of y

x = randn(n,1);
y = randn(m,1);

%% Matlab convolution function
c1 = conv(x,y)

%% Convolution using Toeplitz multiplication
% The toeplitz matrix Y comes from the vector y

r = [y(1); zeros(length(x)-1,1)]; %row vector of Toep (conv length = n+m-1)
c = [y; zeros(length(x)-1,1)]; %column vector of Toep
v = [flipud(r(2:end)) ; c]; %setup toepMult input

c2 = toepMult(v,x)

%% Convolution by extending Y to be circulant
%to get the circulant take v from previous section and shift down rows
%result will be toeplitz with column vector = v

l = length(v);
newv = [v(2:end) ; v]; %setup for new toepMult
newx = [x ; zeros(l-n,1)]; %fill zeros to match dimension
c3 = toepMult(newv,x)
c3 = flipud(c3(end:-1:m+n-1)); %some extra 0 are attached????

%% Convolution using FFT
%padding zeros to fix lengths
x = [x ; zeros(m-1,1)];
y = [y ; zeros(n-1,1)];

c4 = ifft(fft(x).*fft(y))