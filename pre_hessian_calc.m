


function [linear_term, log_term] = pre_hessian_calc,
%%put a flag here to see if it needs to compute hessian!!!!
%%save the coefficient as sparse matrix and use the formula!!!
%     disp("enter");
%     disp(a_hat);
    global A_potential
    global A_bad
    global C
    global type_diffusion
    global i
    
    num_nodes = 1024;
    linear_term = zeros(num_nodes,1);
    log_term = zeros(size(C,1),num_nodes);
    
    for j=1:num_nodes,
        if (A_potential(j,i) > 0)
            linear_term(j) = A_potential(j,i) + A_bad(j,i);
        end
    end

    for c=1:size(C, 1),
        idx = find(C(c,:)~=-1); % used nodes
        [val, ord] = sort(C(c, idx));

        idx_ord = idx(ord);
        cidx = find(idx_ord==i);
        if (~isempty(cidx) && cidx > 1)
            if (strcmp(type_diffusion, 'exp')) 
                log_term(c,idx_ord(1:cidx-1)) = 1;
            end
        end
    end

     %disp(obj);
%         disp(grad);
%gradient = gradient/1000;
%disp(gradient(2));
%disp(a_hat(2));
%     for u = 1:1024,
%         if (a_hat(u) == 0)
%             gradient(u) = 0;
%             continue;
%         end
        
%         a_up = a_hat(:);
%         a_up(u) = a_up(u) + 0.00000001;
%         a_low = a_hat(:);
%         a_low(u) = a_low(u) - 0.00000001;
%         
%         gradient(u) = (estimate_network_copy(a_up) - estimate_network_copy(a_low))/0.00000002;
%     end
%     disp("grad");
%     disp(gradient);
    
end