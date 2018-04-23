% Fast calculate product T*w where T is
% Toeplitz matrix constructed from u according to:
% toeplitz( u(m:n), u(m:-1:1) )
% with complexity O(nlogn) due to fft

% inputs
% u vector n x 1
% w vector of length m x 1. m<=n

function y = toepMult(u,w)

n = length(u);
m = length(w);

%setup vectors
a = u(m:n);
l = length(a);

b = u(m:-1:1);
b = b.';

c = [a; 0; fliplr(b(2:end)).'];
%end setup vectors

%multiplication
p = ifft(fft(c).*fft([w; zeros(l,1)]));

%pick the correct values
y = p(1:l);