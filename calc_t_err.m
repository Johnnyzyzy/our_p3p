function [t_err,fail_flag] = calc_t_err(t_gt, t_est)
    fail_flag = 0;
    num_sols = size(t_est,2);
    t_err = zeros(1,num_sols);
    for i=1:num_sols
        if sum(sum(isinf(t_est(:,i)))) || sum(sum(isnan(t_est(:,i)))) % reject inf or nan
            t_err(i) = inf;
            continue
        else
            t_err(i) = norm(t_gt-t_est(:,i), 1);
        end
    end
    
    t_err(isinf(t_err)) = [];
    if isempty(t_err)
       fail_flag = 1;
    end
end