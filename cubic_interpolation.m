function [lambda] = cubic_interpolation(f , t0 , eta)
% Get the variable used
variable = symvar(f);
grad = gradient(f);


% need to find A and B
A = 0;
B = t0;

g_f_b = subs(grad,variable,B);

% find the value of B at which the derivative is non-negative 
while (g_f_b <= 0)
    
    A = t0;
    t0 = t0 * 2;
    B = t0;
    g_f_b = subs(grad,variable,B);
    
end

fA = vpa(subs(f , variable , A));
fB = vpa(subs(f , variable , B));

g_f_a = subs(grad,variable,A);

Z = 3 * ( fA - fB ) / B + g_f_a + g_f_b;
Q = (Z^2 - g_f_a * g_f_b)^(1/2);

% discarding the negative value
lambda_1 = A + (g_f_a + Z + Q )/(g_f_a + g_f_b + 2*Z)*(B - A);

df_lambda = subs(grad,variable,lambda_1);



while( abs(df_lambda) > eta)
    
    if df_lambda > 0
        B = lambda_1;
        g_f_b = subs(grad,variable,B);
        fB = vpa(subs(f , variable , B));
    else
        A = lambda_1;
        g_f_a = subs(grad,variable,A);
        fA = vpa(subs(f , variable , A));
    end
    
    Z = 3 * ( fA - fB ) / B + g_f_a + g_f_b;
    Q = (Z^2 - g_f_a * g_f_b)^(1/2);
    
    % discarding the negative value
    lambda_1 = A + (g_f_a + Z + Q )/(g_f_a + g_f_b + 2*Z)*(B - A);
    df_lambda = subs(grad,variable,lambda_1) ;
    
end
lambda = lambda_1;

end