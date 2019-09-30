function [fn,qn,q0,fmean,mqn,gam,psi,stats] = time_warping(f,t,lambda,option)
% default options
% option.parallel = 0; % turns offs MATLAB parallel processing (need
% parallel processing toolbox)
% option.closepool = 0; % determines wether to close matlabpool
% option.smooth = 0; % smooth data using standard box filter
% option.sparam = 25; % number of times to run filter
% option.showplot = 1; % turns on and off plotting

if nargin < 3
    lambda = 0;
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
elseif nargin < 4
    option.parallel = 0;
    option.closepool = 0;
    option.smooth = 0;
    option.sparam = 25;
    option.showplot = 1;
end

% time warping on a set of functions
if option.parallel == 1
    if matlabpool('size') == 0
        matlabpool open
    end
end
%% Parameters

fprintf('\n lambda = %5.1f \n', lambda);

a = size(t,1);
if (a ~=1)
    t = t';
end
binsize = mean(diff(t));
[M, N] = size(f);
f0 = f;

if option.smooth == 1
    for r = 1:option.sparam
        for i = 1:N
            f(2:(M-1),i) = (f(1:(M-2),i)+2*f(2:(M-1),i) + f(3:M,i))/4;
        end
    end
end

if option.showplot == 1
    figure(1); clf;
    plot(t, f, 'linewidth', 1);
    title('Original data', 'fontsize', 16);
    pause(0.1);
end

%% Compute the q-function of the plot
addpath(genpath('bspline_tools'))
fy = zeros(M,N);
for ii = 1:N
    y = Bspline(f(:,ii),3);
    ydiff = diff(y);
    fy(:,ii) = ydiff(1:length(t))/binsize;
    f(:,ii) = y(1:length(t));
end
q = fy./sqrt(abs(fy)+eps);

%% Set initial using the original f space
fprintf('\nInitializing...\n');
mnq = mean(q,2);
dqq = sqrt(sum((q - mnq*ones(1,N)).^2,1));
[~, min_ind] = min(dqq);
mq = q(:,min_ind);
mf = f(:,min_ind);

gam = zeros(N,size(q,1));
parfor k = 1:N
    q_c = q(:,k,1)'; mq_c = mq';
    [G,T] = DynamicProgrammingQ2(mq_c/norm(mq_c),t,q_c/norm(q_c),t,t,t,lambda);
    gam0 = interp1(T,G,t);
    gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
end
gamI = SqrtMeanInverse(gam);
gamI_dev = gradient(gamI, 1/(M-1));
mf = interp1(t, mf, (t(end)-t(1)).*gamI + t(1))';
mq = gradient(mf, binsize)./sqrt(abs(gradient(mf, binsize))+eps);

%% Compute Mean
fprintf('Computing Karcher mean of %d functions in SRVF space...\n',N);
ds = inf;
MaxItr = 20;
qun = zeros(1,MaxItr);
for r = 1:MaxItr
    fprintf('updating step: r=%d\n', r);
    if r == MaxItr
        fprintf('maximal number of iterations is reached. \n');
    end
    
    % Matching Step
    clear gam gam_dev;
    % use DP to find the optimal warping for each function w.r.t. the mean
    gam = zeros(N,size(q,1));
    gam_dev = zeros(N,size(q,1));
    parfor k = 1:N
        q_c = q(:,k,1)'; mq_c = mq(:,r)';
        [G,T] = DynamicProgrammingQ2(mq_c/norm(mq_c),t,q_c/norm(q_c),t,t,t,lambda);
        gam0 = interp1(T,G,t);
        gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
        gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));
        f_temp(:,k) = interp1(t, f(:,k,1), (t(end)-t(1)).*gam(k,:) + t(1))';
        q_temp(:,k) = gradient(f_temp(:,k), binsize)./sqrt(abs(gradient(f_temp(:,k), binsize))+eps);
    end
    q(:,:,r+1) = q_temp;
    f(:,:,r+1) = f_temp;
    
    ds(r+1) = sum(simps(t, (mq(:,r)*ones(1,N)-q(:,:,r+1)).^2)) + ...
        lambda*sum(simps(t, (1-sqrt(gam_dev')).^2));
    
    % Minimization Step
    % compute the mean of the matched function
    mq(:,r+1) = mean(q(:,:,r+1),2);
    
    qun(r) = norm(mq(:,r+1)-mq(:,r))/norm(mq(:,r));
    if qun(r) < 1e-2 || r >= MaxItr
        break;
    end
end

% last step with centering of gam
r = r+1;
parfor k = 1:N
    q_c = q(:,k,1)'; mq_c = mq(:,r)';
    [G,T] = DynamicProgrammingQ2(mq_c/norm(mq_c),t,q_c/norm(q_c),t,t,t,lambda);
    gam0 = interp1(T,G,t);
    gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
    gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));
end
gamI = SqrtMeanInverse(gam);
gamI_dev = gradient(gamI, 1/(M-1));
mq(:,r+1) = interp1(t, mq(:,r), (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
for k = 1:N
    q(:,k,r+1) = interp1(t, q(:,k,r), (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
    f(:,k,r+1) = interp1(t, f(:,k,r), (t(end)-t(1)).*gamI + t(1))';
    gam(k,:) = interp1(t, gam(k,:), (t(end)-t(1)).*gamI + t(1));
end

%% Aligned data & stats
fn = f(:,:,r+1);
qn = q(:,:,r+1);
q0 = q(:,:,1);
mean_f0 = mean(f0, 2);
std_f0 = std(f0, 0, 2);
mean_fn = mean(fn, 2);
std_fn = std(fn, 0, 2);
mqn = mq(:,r+1);
fmean = mean(f0(1,:))+cumtrapz(t,mqn.*abs(mqn));

fgam = zeros(M,N);
for ii = 1:N
    fgam(:,ii) = interp1(t, fmean, (t(end)-t(1)).*gam(ii,:) + t(1));
end
var_fgam = var(fgam,[],2);

stats.orig_var = trapz(t,std_f0.^2);
stats.amp_var = trapz(t,std_fn.^2);
stats.phase_var = trapz(t,var_fgam);

gam = gam.';
[~,fy] = gradient(gam,binsize,binsize);
psi = sqrt(fy+eps);

if option.showplot == 1
    figure(2); clf;
    plot((0:M-1)/(M-1), gam, 'linewidth', 1);
    axis square;
    title('Warping functions', 'fontsize', 16);
    
    figure(3); clf;
    plot(t, fn, 'LineWidth',1);
    title(['Warped data, \lambda = ' num2str(lambda)], 'fontsize', 16);
    
    figure(4); clf;
    plot(t, mean_f0, 'b-', 'linewidth', 1); hold on;
    plot(t, mean_f0+std_f0, 'r-', 'linewidth', 1);
    plot(t, mean_f0-std_f0, 'g-', 'linewidth', 1);
    title('Original data: Mean \pm STD', 'fontsize', 16);
    
    figure(5); clf;
    plot(t, mean_fn, 'b-', 'linewidth', 1); hold on;
    plot(t, mean_fn+std_fn, 'r-', 'linewidth', 1);
    plot(t, mean_fn-std_fn, 'g-', 'linewidth', 1);
    title(['Warped data, \lambda = ' num2str(lambda) ': Mean \pm STD'], 'fontsize', 16);
    
    figure(6); clf;
    plot(t, fmean, 'g','LineWidth',1);
    title(['f_{mean}, \lambda = ' num2str(lambda)], 'fontsize', 16);
end
if option.parallel == 1 && option.closepool == 1
    matlabpool close
end