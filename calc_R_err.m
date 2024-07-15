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
