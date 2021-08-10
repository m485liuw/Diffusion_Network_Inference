disp('=== start newton method === ');
tic
%%To study: From where to where convex?? piecewise convex

% use bisection on breakpoints in line 92


network = 'kronecker-core-periphery-n1024-h10-r0_01-0_25-network.txt';
cascades = 'kronecker-core-periphery-n1024-h10-r0_01-0_25-1000-cascades.txt';
horizon = 10;
type_diffusion = 'exp';
num_nodes = 1024;
A = create_adj_matrix(network, num_nodes);
C = create_cascades(cascades, num_nodes);
X_hat = zeros(size(A));

total_obj = 0;

example = matfile('saveAp.mat');
A_potential = example.A_potential;
example2 = matfile('saveAb.mat');
A_bad = example2.A_bad;
example3 = matfile('savecas.mat');
num_cascades = example3.num_cascades;

n = 1024;
l   = zeros(n,1);
u = inf(n,1);
threshold = eps;

for i = 1:n,
    fprintf('\n --- i = %d --- \n', i);
    x_0 = ones(n,1)*1;
    x_0(A_potential(:,i)==0) = 0;
    hessian_calc = make_handle_hessian_calc(i, C, type_diffusion);
    
    if (num_cascades(i)==0)
        X_hat(:,i) = 0;
        continue;
    end


    [ob0, grad0, scrap] = hessian_calc(x_0,false);
    MAXITER =  200000;
    x = x_0;
    x_last = x;
    for iter=1:MAXITER,
        disp("iter:" + int2str(iter));
        [obj, grad, hessian] = hessian_calc(x, true);
        pgrad = grad .* (grad < 0 | x > sqrt(eps));
        if norm(pgrad) < sqrt(eps) * norm(grad0)
            break
        end
        hessian = hessian + (eps^.66 * norm(hessian,Inf))*eye(n); %not useful
        %defien binding set
        %B = find(x <= sqrt(eps) & (grad>0));
        B = find(x <= norm(x) * 5 * eps & grad > 0);
        ismem = ismember(1:n,B);
        F = find(~ismem);
        %if length of F == 0  --> terminate!!!!
        if isempty(F),
            break;
        end 
               
        %principal submatrix
        hessian = hessian(F,F);
        S = inv(hessian);
        %line search
        %disp("enter line search");
        temp_h = hessian\grad(F);
        h = zeros(n,1);
        h(F) = temp_h;
        % Must have a descent direction
        assert(grad' * h >= 0)
        
        t = x./h;
        subset = find(h > 0);  
        pos_t1 = t(subset);
        pos_t_ind0 = (1:n)';        
        pos_t_ind1 = pos_t_ind0(subset);
        [pos_t2,ord] = sort(pos_t1);
        pos_t_ind2 = pos_t_ind1(ord);
        pos_t = [0; pos_t2; Inf];
        pos_t_ind = [0;pos_t_ind2];      
        
        begin_j = 1;
        end_j = length(pos_t);
        special_j = 0;
        num_j_it = 0;
        while begin_j + 1 < end_j
            num_j_it = num_j_it + 1;
            j = fix((end_j+begin_j)/2);
            %disp("line searching index" + int2str(j));
            bp = pos_t(j);
            xtest = x - bp * h;
            assert(abs(xtest(pos_t_ind(j))) <= (norm(x) + norm(h) * abs(bp)) * 1e-15);
            xtest(pos_t_ind(2:j)) = 0;  
            [obj_test, grad_test, ignore] = hessian_calc(xtest,false);
            %fprintf('intv = (%d %d %d) bp = %e norm(x) = %e obj_test = %e\n', ...
            %    begin_j, j, end_j, bp, norm(x), obj_test);
   
            if isinfornan(obj_test)
                end_j = j;
                continue
            end
            hleft = h;
            hleft(pos_t_ind(2 : j-1)) = 0;
            hright = h;
            hright(pos_t_ind(2 : j)) = 0;
            phi_prime_l = -grad_test' * hleft;
            phi_prime_r = -grad_test' * hright;
            assert(~isinfornan(phi_prime_l) && ~isinfornan(phi_prime_r))
            if phi_prime_l <= 0 && phi_prime_r >= 0
                special_j = j;
                break
            end
            if phi_prime_l >= 0
                end_j = j;
            else
                begin_j = j;
            end
        end
        if special_j > 0
            %fprintf('special_j = %d\n', special_j);
            bp = pos_t(special_j);
            zind = pos_t_ind(2:special_j);
            num_bis_it = 0;
        else            
            %disp("enter consec bisection");
            bp_begin_j = pos_t(begin_j);
            h_begin_j = h(:);
            h_begin_j(pos_t_ind(2:begin_j)) = 0;
            [obj_begin_j, grad_begin_j, scrap] = hessian_calc(x-bp_begin_j*h,false);                             
            phi_prime_begin_j = grad_begin_j'*h_begin_j*(-1);
            assert(phi_prime_begin_j <= 0)
            bp_end_j = pos_t(end_j);
            if isinfornan(bp_end_j)
                for src = 1 : 100
                    bp_end_j = bp_begin_j + 2^src;
                    [obj_end_j, grad_end_j, scrap] = hessian_calc(x - bp_end_j*h, false);
                    phi_prime_end_j = grad_end_j'*h_begin_j*(-1);
                    if phi_prime_end_j >= 0 || isinfornan(phi_prime_end_j)
                        break
                    end
                    assert(src < 100)
                end
            else                
                [obj_end_j, grad_end_j, scrap] = hessian_calc(x - bp_end_j*h, false);
                phi_prime_end_j = grad_end_j'*h_begin_j*(-1);
            end
                
            assert(phi_prime_end_j>=0 || isinfornan(phi_prime_end_j));
            
            delta = bp_end_j - bp_begin_j;
            a1 = obj_begin_j;
            a2 = obj_end_j;
            b1 = phi_prime_begin_j;
            b2 = phi_prime_end_j;
            
            if begin_j == 1 && bp_end_j >= 1
                bp = 1;
                num_bis_it = -4; %pure newton's method
            elseif isinfornan(phi_prime_end_j) || b2 > 1e6 * abs(b1)
                %fprintf('a = %e b = %e scale1 = %e scale2 = %e obj1 = %e obj2 = %e\n', ...
                %    a, b, scale1, scale2, obj_begin_j, obj_end_j);
                %[a,b, num_bis_it] = bisection(a,b,scale1,scale2,threshold,x,h,hessian_calc);
                %bp = (a+b) / 2;
                
                [obj_begin_j, grad_begin_j, hessian_begin_j] = ...
                    hessian_calc(x-bp_begin_j*h,true);
                c1 = h_begin_j' * hessian_begin_j * h_begin_j / 2;
                z = -c1 * delta^2;
                xx = b1 + z / delta;
                yy = a1 - z * log(delta);
                bp = bp_end_j - z / xx;
                num_bis_it = -3;
                assert(bp >= bp_begin_j && bp < bp_end_j)
            elseif (b1 + b2) / 2 > (a2 - a1) / delta + (abs(a1) + abs(a2)) * 1e-6 / delta
                   
                kwbar = @(wbar)(b1 + (b1 - b2) * wbar + (a1 - a2) / delta - ...
                    (b1 - b2) * wbar .* (1 + wbar) .* log(1 + 1./wbar));
                ub = 1;
                while kwbar(ub) < 0
                    ub = ub * 2;
                    if ub > 1e6
                        keyboard
                    end
                end
                lb = 0;
                itcount = 0;
                while ub - lb > max(1e-13, 1e-13 * lb)
                    mid = (lb + ub) / 2;
                    if kwbar(mid) >= 0
                        lb = mid;
                    else
                        ub = mid;
                    end
                    if itcount > 100
                        keyboard
                    end
                    itcount = itcount + 1;
                end
                wbar = mid;
                w = wbar * delta;
                z = (b1 - b2) * w * (delta + w) / delta;
                xx = b1 + z / (delta + w);
                y = a1 - z * log(delta + w);
                bp = bp_end_j + w - z / xx;
                assert(bp >= bp_begin_j && bp <= bp_end_j)
                num_bis_it = -2;
            elseif (b1 + b2) / 2 < (a2 - a1) /delta - (abs(a1) + abs(a2)) * 1e-6 / delta
                hwbar = @(wbar)(-b1 * wbar + (1 + wbar) * b2 + (a1 - a2) / delta + ...
                    (b1 - b2) * wbar .* (1 + wbar) .* log(1 + 1./wbar));
                ub = 1;
                while hwbar(ub) > 0
                    ub = ub * 2;
                    if ub > 1e6
                        keyboard
                    end
                end
                lb = 0;
                itcount = 0;
                while ub - lb > max(1e-13, 1e-13*lb)
                    mid = (lb + ub) / 2;
                    if hwbar(mid) >= 0
                        lb = mid;
                    else
                        ub = mid;
                    end
                    if itcount > 100
                        keyboard
                    end
                    itcount = itcount + 1;
                end
                wbar = mid;
                w = wbar * delta;
                z = (b1 - b2) * w * (delta + w) / delta;
                xx = b1 - z / w;
                y = a1 - z * log(w);
                bp = bp_begin_j - w - z / xx;
                assert(bp >= bp_begin_j && bp <= bp_end_j)
                num_bis_it = -1;
                
                
                
                
            else
                a = bp_begin_j;
                b = bp_end_j;
                scale1 = phi_prime_begin_j;
                scale2 = phi_prime_end_j;
                [a,b, num_bis_it] = bisection(a,b,scale1,scale2,threshold,x,h,hessian_calc);
                bp = (a+b) / 2;
            end
            
            zind = pos_t_ind(2:begin_j);
        end
            
        x_last = x;
        x = x - bp*h;
        x(zind)=0;
        x(A_potential(:,i)==0) = 0;
        assert(all(x >= 0))
        [obje, grade, scrap] = hessian_calc(x, false);
        pgrade = grade .* (grade < 0 | x > sqrt(eps));
        
        fprintf('at end of iteration, norm(pgrad) = %e, num_j_it = %d num_bis_it = %d bp = %e\n', ...
           norm(pgrade), num_j_it, num_bis_it, bp);
        %fprintf('at end of iteration norm(x) = %e obj = %e norm(grad) = %e, norm(pgrad) = %e\n', ...
        %    norm(x), obje, norm(grade), norm(pgrade));
    end
    X_hat(:,i) = x;
    
end
toc


