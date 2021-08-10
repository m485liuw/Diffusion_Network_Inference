function [a,b,numit] = bisection(a,b,scale1,scale2,threshold,x,h,hessian_calc)
numit = 0;
 while abs(a*abs(scale1)-b*abs(scale2))>threshold || isinfornan(scale2),
     numit = numit + 1;
                c=(a+b)/2;
                if c == a || c == b;
                    break;
                end
                [obj_c,grad_c,hessian_c] = hessian_calc(x-c*h,false);
                h_c = h(:);
                h_c(x-c*h < 0) = 0;
                phi_prime_c = grad_c.'*h_c*(-1);
                [obj_a,grad_a,hessian_a] = hessian_calc(x-a*h,false);
                h_a = h(:);
                h_a(x-a*h < 0) = 0;
                phi_prime_a = grad_a.'*h_a*(-1);
                if phi_prime_c*phi_prime_a < 0 || isinfornan(phi_prime_c),
                    b = c;
                    scale2 = phi_prime_c;
                else
                    a = c;
                    scale1 = phi_prime_c;
                end
 end
