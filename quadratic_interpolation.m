function [lambda] = quadratic_interpolation(f , t0 , eta)

% get the variable used
variable = symvar(f);

% first fit

A = 0;

fA = vpa(subs(f , variable , A));
f1 = vpa(subs(f , variable , t0));

if f1 > fA
    fC = f1;
    t0 = t0 / 2;
    fB = vpa(subs(f , variable , t0));
else
    while(f1<fA)
        fB = f1;
        f2 = vpa(subs(f , variable , 2*t0));
        
        if f1 < f2
            fC = f2;
            break
        elseif f1 > f2
            f1 = f2;
            t0 = 2*t0;
        end
    end
    
end

B = t0;
C = 2*t0;

a = fA;
b = (4*fB - 3*fA - fC)/(2*t0);
c = (fC + fA - 2 * fB)/(2 * t0^2);

h = a + b * variable + c * variable^2;

lambda = -b/(2*c);

f_lambda = vpa(subs(f , variable , lambda));
h_lambda = vpa(subs(h , variable , lambda));

tol =abs((f_lambda - h_lambda)/(f_lambda));


% re-fit 
while (tol > eta )
    if lambda > B && f_lambda <fB
        A = B;
        B = lambda;
    elseif lambda > B && f_lambda > fB
        C = lambda;
    elseif lambda < B && f_lambda < fB
        C = B;
        B = lambda;
    else
        A = lambda;
    end
    
    fA = vpa(subs(f , variable , A));
    fB = vpa(subs(f , variable , B));
    fC = vpa(subs(f , variable , C));
    
    a = (fA*B*C*(C-B) + fB*C*A*(A - C) + fC*A*B*(B - A))/((A-B)*(B -C)*(C -A));
    b = (fA*(B^2 - C^2) + fB*(C^2 - A^2) + fC*(A^2 - B^2))/((A-B)*(B -C)*(C -A));
    c = -(fA * (B - C) + fB * (C - A) + fC * (A - B))/((A-B)*(B -C)*(C -A));
    
    lambda = -b/(2*c);
    h = a + b * variable + c * variable^2;
    
    f_lambda = vpa(subs(f , variable , lambda));
    h_lambda = vpa(subs(h , variable , lambda));
    
    tol =abs((f_lambda - h_lambda)/(f_lambda));
    
end

end