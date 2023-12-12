N = 8;
M = 4000;
theta = linspace(-1,1,M)';
u = zeros(N,M);
f = zeros(N,1);


phi = -pi / 6;
for i = 1:1:N
    for m = 1:1:M
        u(i,m) = 1/N*exp(1j*pi*(i-1)*theta(m));
        f(i) = exp(1j*pi*(i-1)*sin(phi));
    end

end

g = u'*f;

plot(theta,abs(g))

