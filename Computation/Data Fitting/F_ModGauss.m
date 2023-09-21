function [Gaussian] = F_ModGauss(x, a, b, c)

Gaussian = [];

for ix = x
    x_ = a*exp(-((ix-b)/c)^2);
    Gaussian = [Gaussian, x_];

end

