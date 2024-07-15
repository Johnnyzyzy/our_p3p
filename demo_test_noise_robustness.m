close all; clear; clc;
addpath(genpath('P3P_tool'))
addpath others
addpath(genpath('My_P3P_cubic'))

solvers = {@(m,X)p3p_gao(m,X), @(m,X)p3p_kneip(m,X),@(m,X)p3p_banno(m,X),...
    @(m,X)p3p_ke(m,X,1), @(m,X)p3p_nakano_bmvc2019(m,X), ...
    @(m,X)p3p_lambdatwist(m,X,1), @(m,X)my_cubic_p3p(m,X,0,1,0)};
methods = {'Gao', 'Kneip', 'Banno', 'Ke', 'Nakano ','LambdaTwist', 'Ours'};
marker= {'>', 'o', '+','*','^','x','s'};
color= {'m','k','c','r',[255,187,0]/255,'b','g'};
markerfacecolor=  {'n','n','n','n','n','n','n'};
linestyle= {'-','-','-','-','-','-','-'};
% rng(2)

npt= 3;width= 640;height= 480;f= 800; 
nls = 0:1:5;

num = 5000;
num_solvers = length(solvers);
num_noise_level = length(nls);
err_list_r_all = cell(1,num_noise_level); err_list_t_all = cell(1,num_noise_level);
fail_flag = zeros(num_solvers,num,num_noise_level);
mean_r = zeros(num_solvers,num_noise_level); median_r = mean_r; mean_t = mean_r; median_t = mean_r;
for i = 1:num_noise_level
    nl = nls(i);
    err_list_r_temp = zeros(num_solvers,num); err_list_t_temp = zeros(num_solvers,num);
    for j = 1:num
        Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
            
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
        for k=1:num_solvers
            [R_est, t_est] = solvers{k}(normalized_pts2d, pts3d);
            if isempty(R_est) || isempty(t_est)
                fail_flag(k,j,i) = 1; % record failures
                continue
            else
                [R_err,fail_flag_r] = calc_R_err(R, R_est);
                [t_err,fail_flag_t] = calc_t_err(t, t_est);
                if fail_flag_r || fail_flag_t
                    fail_flag(k,j,i) = 1; % record failures
                    continue
                else
                    [~,ind_min] = min(R_err+t_err); 
                    err_list_r_temp(k,j) = R_err(ind_min); % min(R_err)  R_err(ind_min)
                    err_list_t_temp(k,j) = t_err(ind_min); % min(t_err)  t_err(ind_min)
                end
            end
        end
        showpercent(j,num);    
    end
    fprintf('\n');
    % delete the columns where some of the solvers fail
    fail_index = logical(sum(fail_flag(:,:,i)));
    err_list_r_temp(:,fail_index) = [];  err_list_t_temp(:,fail_index) = []; 
    err_list_r_all{i} = err_list_r_temp;  err_list_t_all{i} = err_list_t_temp;  
    % record mean and median errors
    mean_r(:,i) = mean(err_list_r_temp,2);
    median_r(:,i) = median(err_list_r_temp,2);
    mean_t(:,i) = mean(err_list_t_temp,2);
    median_t(:,i) = median(err_list_t_temp,2);
end

% draw mean error
figure('color','w','position',[100 100 560 420]); box on; hold on;
for i= 1:num_solvers
    plot(nls,mean_r(i,:),'marker',marker{i},'color',color{i},...
        'markerfacecolor',markerfacecolor{i},'displayname',methods{i}, ...
        'LineWidth',2,'MarkerSize',10,'LineStyle',linestyle{i});
end
title('Mean Rotation Error','FontSize',12,'FontName','Arial');
ylabel('Absolute Error','FontSize',11);
xlabel('Image Noise (pixel)','FontSize',11);
legend('Location','NorthWest');

figure('color','w','position',[800 100 560 420]); box on; hold on;
for i= 1:num_solvers
    plot(nls,mean_t(i,:),'marker',marker{i},'color',color{i},...
        'markerfacecolor',markerfacecolor{i},'displayname',methods{i}, ...
        'LineWidth',2,'MarkerSize',10,'LineStyle',linestyle{i});
end
title('Mean Translation Error','FontSize',12,'FontName','Arial');
ylabel('Absolute Error','FontSize',11);
xlabel('Image Noise (pixel)','FontSize',11);
legend('Location','NorthWest');


