clear;clc;
addpath(genpath('P3P_tool'));addpath others;addpath(genpath('My_P3P_cubic'));

solvers = {@(m,X)p3p_gao(m,X), @(m,X)p3p_kneip(m,X),@(m,X)p3p_banno(m,X),...
    @(m,X)p3p_ke(m,X,1), @(m,X)p3p_nakano_bmvc2019(m,X), ... 
    @(m,X)p3p_lambdatwist(m,X,1),@(m,X)my_cubic_p3p(m,X,0,1,1)};
methods = {'Gao', 'Kneip', 'Banno',  'Ke', 'Nakano', 'LambdaTwist','Ours'};        
npt= 3;
nl = 0;
width= 640;
height= 480;
f= 800;  

num = 1e5;
num_solvers = length(solvers);
time_list = zeros(num_solvers,num);
for i = 1:num
    Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
    d = Xc(3,:);
    t= median(Xc,2);
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
    for j = 1:num_solvers
        tic
        [R_est, t_est] = solvers{j}(normalized_pts2d, pts3d);
        time_list(j,i) = toc;
    end

    showpercent(i,num);
end

% runing time compare
% median(time_list,2)'  
mean_time = mean(time_list,2)';
for i=1:num_solvers
    disp([ methods{i}, '   --> mean runtime (microseconds): ', num2str(mean_time(i)*1e6)])
end


