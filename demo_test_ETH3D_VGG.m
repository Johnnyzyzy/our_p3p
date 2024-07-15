clear;clc;close all;warning off;format long g;
addpath(genpath('./P3P_tool'))
addpath(genpath('./My_P3P_cubic'))
addpath ./others
addpath ./experiments

%% ETH3D dataset (high_resolution, undistorted) choose scene:
%% you can choose one scene or several scenes
%% 'facade' data is larger than 25mb and can't be uploaded to github directly, we will upload it later
% scene_choose = {'courtyard','delivery_area','electro','kicker',...
%                  'meadow', 'office', 'pipes', 'playground', 'relief',...
%                  'relief_2', 'terrace', 'terrains'}; % 'facade',
scene_choose = {'courtyard','electro','kicker',...
                 'office', 'pipes', 'terrace'}; % 'facade',
% scene_choose = {'courtyard'};

%% VGG dataset choose scene:
%% you can choose one scene or several scenes
% scene_choose = {'MH', 'Corridor', 'MC1', 'MC2', 'MC3', 'UL', 'WC'};
% scene_choose = {'WC'};

%%
solvers = {@(m,X)p3p_gao(m,X), @(m,X)p3p_kneip(m,X),@(m,X)p3p_banno(m,X),...
    @(m,X)p3p_ke(m,X,1), @(m,X)p3p_nakano_bmvc2019(m,X), ...
    @(m,X)p3p_lambdatwist(m,X,1), @(m,X)my_cubic_p3p(m,X,0,1,0)};
methods = {'Gao', 'Kneip', 'Banno', 'Ke', 'Nakano ','LambdaTwist', 'Ours'};

num_methods = length(methods);
all_mean_R = [];  all_mean_t = [];  all_fails = [];  all_mean_time = [];
for ind_scene = 1:length(scene_choose) % 
    eval(['load ETH3D_hreso_undist_',scene_choose{ind_scene}]); % for ETH3D
%     eval(['load VGG_',scene_choose{ind_scene}]); % for VGG
    
    num_test = 500;
    num_frame = length(imagePoints);
    
%     sum = 0;  
%     for jjj = 1:num_frame
%         sum = sum + length(imagePoints{jjj});
%     end
%     average_point = sum/num_frame

%     rng(3)

    err_R = zeros(num_frame,num_test,num_methods);  err_t = zeros(num_frame,num_test,num_methods);
    time_list = zeros(num_frame,num_test,num_methods);
    fail_record = zeros(num_frame,num_methods);
    for ind_frame = 1:num_frame
        ind_frame
        Rt = Rts(:,:,ind_frame);  tt = tts(:,ind_frame);  A = As(:,:,ind_frame);
        xxd = imagePoints{ind_frame};  
        xxd = A\[xxd;ones(1,size(xxd,2))]; 
        XXw = WorldPoints{ind_frame};  
        num_points = size(imagePoints{ind_frame},2);

        % filter out points with large reprojection error (only for ETH3D)
%         XXw_err = WorldPoints_error{ind_frame}; ind_r = XXw_err < 2; 
        ind_r = true(1,num_points); % no filtering
        XXw_accu = XXw(:,ind_r); xxd_accu = xxd(:,ind_r); 
       
        normalized_pts2d = xxd_accu./vecnorm(xxd_accu);
        pts3d = XXw_accu;

        % num_test trails for each frame
        % if any method fails, run again
        num_points_accu = size(XXw_accu,2);
        ind_test = 1;
        while(ind_test <= num_test)
            ind_points = randperm(num_points_accu,3);       
            fail_flag = 0;

            for ind_method = 1:num_methods
                tic
                [Re, te] = solvers{ind_method}(normalized_pts2d(:,ind_points), pts3d(:,ind_points));
                time_temp = toc;
                if isempty(Re) || isempty(te)
                    fail_flag=1; fail_record(ind_frame,ind_method) = fail_record(ind_frame,ind_method)+1;
                    continue
                end       

                num_sols = size(te,2);
                [err_R_temp,fail_flag_R]=calc_R_err(Rt, Re);  
                [err_t_temp,fail_flag_t]=calc_t_err(tt, te);
                if fail_flag_R || fail_flag_t
                    fail_flag=1; fail_record(ind_frame,ind_method) = fail_record(ind_frame,ind_method)+1;
                    continue
                end       
                [~,ind_min] = min(sum([err_R_temp;err_t_temp]));
                
                % failures
                if err_R_temp(ind_min) > 1e3
                    fail_flag=1; fail_record(ind_frame,ind_method) = fail_record(ind_frame,ind_method)+1;
                    continue
                end
                
                err_R(ind_frame,ind_test,ind_method) = err_R_temp(ind_min);  
                err_t(ind_frame,ind_test,ind_method) = err_t_temp(ind_min);
                time_list(ind_frame,ind_test,ind_method) = time_temp;

            end

            if fail_flag == 1 % if any method fails, do not record but run again
               fail_flag = 0;
               continue
            end
            ind_test = ind_test + 1;
        end
    end
    mean_err_Rt = [mean( permute(sum(err_R,2),[1 3 2]) ./ num_test );...
                                            mean( permute(sum(err_t,2),[1 3 2]) ./ num_test )]';
%     mean_time = [mean( permute(sum(time_list,2),[1 3 2]) ./ num_test )]';
%     fails = sum(fail_record)';
    
    all_mean_R = [all_mean_R mean_err_Rt(:,1)];
    all_mean_t = [all_mean_t mean_err_Rt(:,2)];
%     all_mean_time = [all_mean_time mean_time * 1e6]; % microseconds
%     all_fails = [all_fails fails];    
end
all_mean_R
all_mean_t