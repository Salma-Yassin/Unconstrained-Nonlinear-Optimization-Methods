function [x] = golden_section( f , a , b , eta)

% get the variable used 
variable = symvar(f);

while( (b - a) > eta)
    
    d = 0.618*(b - a);
    x1 = a + d;
    x2 = b - d;
    
    f1 = vpa(subs(f , variable , x1));
    f2 = vpa(subs(f , variable , x2));
    
    if (f1 > f2)
        b = x1;
    elseif (f1<f2)
        a = x2;
    else
        a = x2;
        b = x1;
    end
   
end 

x =(x1 + x2)/2;
end 