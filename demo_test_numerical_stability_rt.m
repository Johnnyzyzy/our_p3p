close all;clear;clc;
addpath(genpath('P3P_tool'))
addpath others
addpath(genpath('My_P3P_cubic'))

solvers = {@(m,X)p3p_gao(m,X), @(m,X)p3p_kneip(m,X),@(m,X)p3p_banno(m,X),...
    @(m,X)p3p_ke(m,X,1), @(m,X)p3p_nakano_bmvc2019(m,X), ...
    @(m,X)p3p_lambdatwist(m,X,1), @(m,X)my_cubic_p3p(m,X,0,1,0)}; %  
methods = {'Gao', 'Kneip', 'Banno', 'Ke ', 'Nakano ','LambdaTwist', 'Ours'};
Color = {'m','k','c','r',[255,187,0]/255,'b','g'};
LineStyle = {'-','-','-','-','-','-','-'};
% rng(1) % 

npt= 3; nl = 0; width= 640; height= 480; f= 800; 
num = 1e5;
num_solvers = length(solvers);
err_list_r = zeros(num_solvers,num);
err_list_t = zeros(num_solvers,num);
fail_flag = zeros(num_solvers,num);
num_outputs = zeros(num_solvers,num);
for i = 1:num
    Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])]; % ordinary
    Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-0.05 0.05]); xrand(1,npt,[4 8])]; % thin-flat 1
    Xc= [xrand(1,npt,[-0.05 0.05]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])]; % thin-flat 2
    
    d = Xc(3,:);
    t= mean(Xc,2);
    R= rodrigues(randn(3,1));
    XXw= R\(Xc-repmat(t,1,npt));
    % projection
    xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
    xxn= xx+randn(2,npt)*nl;
    xx_normalized = xxn/f;
    
    normalized_pts2d = [xx_normalized; ones(1,npt)];
    d = d.*vecnorm(normalized_pts2d);
    normalized_pts2d = normalized_pts2d./vecnorm(normalized_pts2d);
    pts3d = XXw;
    
    % run solvers
    for j=1:num_solvers
        [R_est, t_est] = solvers{j}(normalized_pts2d, pts3d);
        R_est = real(R_est);  t_est = real(t_est);
        if isempty(R_est) || isempty(t_est)
            fail_flag(j,i) = 1; % record failures
            continue
        else
            [R_err,fail_flag_r] = calc_R_err(R, R_est);
            [t_err,fail_flag_t] = calc_t_err(t, t_est);
            if fail_flag_r || fail_flag_t
                fail_flag(j,i) = 1; % record failures
                continue
            else
                err_list_r(j,i) = min(R_err);
                err_list_t(j,i) = min(t_err);
                num_outputs(j,i) = length(R_err);
            end
        end
    end
    showpercent(i,num);
end
fail_index = logical(sum(fail_flag));
err_list_r(:,fail_index) = [];err_list_t(:,fail_index) = [];num_outputs(:,fail_index) = []; % delete the columns where some of the solvers fail
% [sum(num_outputs,2),sum(fail_flag,2)]

% draw histogram of rotation error
figure('color','w','position',[100 100 540 480]);
hold all;
box on;
xx = linspace(-18,0,100);
for i = 1:num_solvers
    y_temp = ksdensity(log10(err_list_r(i,:)),xx);
    y_temp = y_temp./sum(y_temp); % avoid y > 1
    plot(xx,y_temp,'Color',Color{i},'LineStyle',LineStyle{i},'LineWidth',3);  
    hold on;
end
xlim([-18 0]);
xtick= -18:2:0;
set(gca,'xtick',xtick);
title('Histogram of Rotation Error');
xlabel('Log_{10} Rotation Error');
ylabel('Frequency');
legend(methods);


% draw histogram of translation error
figure('color','w','position',[600 100 540 480]);
hold all;
box on;
xx = linspace(-18,0,100);
for i = 1:num_solvers
    y_temp = ksdensity(log10(err_list_t(i,:)),xx);
    y_temp = y_temp./sum(y_temp);
    plot(xx,y_temp,'Color',Color{i},'LineStyle',LineStyle{i},'LineWidth',3);  
    hold on;
end
xlim([-18 0]);
xtick= -18:2:0;
set(gca,'xtick',xtick);
title('Histogram of Translation Error');
xlabel('Log_{10} Translation Error');
ylabel('Frequency');
legend(methods);


%% helper functions
function [R_err,fail_flag] = calc_R_err(R_gt, R_est)
    fail_flag = 0;
    num_sols = size(R_est,3);
    R_err = zeros(1,num_sols);
    for i=1:num_sols
        if sum(sum(isinf(R_est(:,:,i)))) || sum(sum(isnan(R_est(:,:,i)))) % reject inf or nan
            R_err(i) = inf;
            continue
        else
            R_err(i) = norm(R_gt-R_est(:,:,i), 1);
        end
    end
    
    R_err(isinf(R_err)) = [];
    if isempty(R_err)
       fail_flag = 1;
    end
end

function [t_err,fail_flag] = calc_t_err(t_gt, t_est)
    fail_flag = 0;
    num_sols = size(t_est,2);
    t_err = zeros(1,num_sols);
    for i=1:num_sols
        if sum(sum(isinf(t_est(:,i)))) || sum(sum(isnan(t_est(:,i)))) % reject inf or nan
            t_err(i) = inf;
            continue
        else
            t_err(i) = norm(t_gt - t_est(:,i),1);
        end
    end
    
    t_err(isinf(t_err)) = [];
    if isempty(t_err)
       fail_flag = 1;
    end
end