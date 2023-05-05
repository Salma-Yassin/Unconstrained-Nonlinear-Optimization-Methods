function [soln, fmin , k] = quasi_newton(x0 , f , eta ,l_option)


grad = gradient(f);
A_exp = hessian(f); % used only for intialization step 
variables =symvar(f);

% intialize the variables
x = x0;
g_f_old = subs(grad,variables,x);

% intialize the value of B to be the Hessian matrix 
B =  vpa(subs(A_exp,variables,x));
f_old =  subs(f,variables,x);
tol = g_f_old'*g_f_old ;
k = 0

while(tol> eta)
    
    % iteration number 
    k = k + 1;
    
    s = - inv(B) * g_f_old;
    
    % find optimal step size  
    syms lambda_sym
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
    
    x = x + lambda * s';
    g_f = subs(grad,variables,x);
    y = g_f - g_f_old;
    
    % the BFGS update 
    B = B -(B*(s*s')*B)/(s'*B*s) + (y*y')/(y'*s);
    tol = g_f'*g_f ;
    g_f_old = g_f;
    
end 

soln = x;
fmin = vpa(subs(f , variables , soln));

end 