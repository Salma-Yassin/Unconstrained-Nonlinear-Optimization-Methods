function [soln , fmin , k] = marquardt(x0 , f , alpha , c1 , c2 , eta,l_option)

grad = gradient(f);
A_exp = hessian(f);
variables =symvar(f);

g_f = subs(grad,variables,x0);

% initialization 
k = 0;
s = - g_f;
x = x0;
tol = vpa(g_f'*g_f);
f_old = vpa(subs(f,variables,x0));
I = eye(length(variables));

while( tol > eta)
    
    k = k + 1;
    
    % find the optimal step length lambda 
    A = vpa(subs(A_exp,variables,x));
    s = - inv([A + alpha * I ]) * g_f;
    syms lambda_sym;
    f_lambda = subs(f , variables , x + lambda_sym * s');
    
    % get the optimal step length
    if l_option == 1
        lambda = golden_section(f_lambda, 0 , 10 , eta);
    elseif l_option == 2 
        lambda =fibonacci_section( f_lambda , 0 , 10 , 0.001);
    elseif l_option == 3
        lambda = quadratic_interpolation(f_lambda , 0.01 , 0.001);
    else 
        lambda = cubic_interpolation(f_lambda, 0.01 , 0.001);
    end
 
    % Find the point xi = xi-1 + lambda si-1 
   
    xi = x + lambda * s';
    fi =  subs(f,variables,xi);
    
    while(fi > f_old)
        alpha = alpha*c2;
        s = - [A + alpha * I ] \ g_f;
        xi = x + lambda * s';
        fi =  subs(f,variables,xi);
    end 
    
    alpha = alpha * c1;
    
    % Find grad f
    g_f = subs(grad,variables,xi);
    tol = g_f'*g_f; 
    x = xi;
    f_old = fi;
    
end 

soln = x;
fmin = vpa(subs(f , variables , soln));

end