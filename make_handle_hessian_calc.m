function hessian_calc = make_handle_hessian_calc(i, C, type_diffusion),


    function [obj, gradient, hessian] = hessian_calc_inner(a_hat, need_h),
    %%put a flag here to see if it needs to compute hessian!!!!

        num_nodes = 1024;
        obj = 0;
        gradient = zeros(num_nodes,1);
        %hessian = zeros(1024,num_nodes);
        hessian = sparse([], [], [], 1024, 1024);
        a_hat(a_hat<0)=0;


        obj = obj + linear_term.'*a_hat;

        gradient = gradient + linear_term;

        for u = 1:size(log_term,1),
            coeff = log_term(u,:);
            if sum(coeff~=0) == 0,
                continue
            end
            var = coeff*a_hat;
            obj = obj - log(var);
            gradient = gradient - coeff.'/(var);
            if need_h,
                ff = find(coeff);
                coeffsp = sparse(ff*0+1,ff,coeff(ff),1, 1024);
                hessian = hessian + (coeffsp.'*coeffsp)/((var)^2);
            end
        end
    end
[linear_term, log_term] = pre_hessian_calc(i, C, type_diffusion);
hessian_calc = @hessian_calc_inner;
        
        
      
end