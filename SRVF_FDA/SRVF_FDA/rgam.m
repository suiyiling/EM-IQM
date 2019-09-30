% subroutine for generating random time warping

function gam = rgam(N, sigma, num, factor)

% N: number of discretized bins
% sigma: weight on the covariance
% num: number of warping functions

gam = zeros(num,N);
for k = 1:num
    epsilon = normrnd(0,sigma);
    if factor == 0
        points = linspace(0+epsilon,1+epsilon,N);
        tmp = points.^2;
    elseif factor == -1
        points = linspace(0+epsilon,1+epsilon,N);
        tmp = (1-points).^2;
    else
        points = linspace(0+epsilon,factor*pi+epsilon,N);
        tmp = sin(points)+1;
    end
    gam(k,:) = cumtrapz(tmp)/trapz(tmp);
end
