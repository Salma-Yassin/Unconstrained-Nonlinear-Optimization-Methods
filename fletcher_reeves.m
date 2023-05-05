function [soln , fmin , k] = fletcher_reeves(x0,f,eta)

grad = gradient(f);
A_exp = hessian(f);
variables =symvar(f);

% 2. Set the first search direction
g_f = subs(grad,variables,x0);

% initialization 
k = 0;
s = -g_f;
x = x0;
tol = g_f'*g_f;  

while( tol > eta)
    
    k = k + 1;
    % find the optimal step length lambda 
    A = subs(A_exp,variables,x);
    num = tol;
    dem = s' * A * s;
    lambda = vpa(num/dem);
    
    % Find the point xi = xi-1 + lambda si-1 
    x = x + lambda * s';
    
    % Find grad f
    g_f_old = g_f;
    g_f = subs(grad,variables,x);
    
    % calculate beta 
    beta = vpa((g_f'*g_f)/(g_f_old'*g_f_old));
    
    % update s
    s = - g_f + beta * s;
    
    tol = g_f'*g_f ;
    
end 

soln = x;
fmin = vpa(subs(f , variables , soln));
end 