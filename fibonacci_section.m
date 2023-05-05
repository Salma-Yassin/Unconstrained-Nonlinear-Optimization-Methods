function [x] = fibonacci_section( f , a , b , eta)

% get the variable used 
variable = symvar(f);

% Intialize the array 
fn = ceil((b - a)/eta);
fn_arr = [1 1];
last = fn_arr(end);

while(fn > last)
    
   last = fn_arr(end) + fn_arr(end - 1);
   fn_arr = [fn_arr last];
   
end


n = length(fn_arr);
k = n - 2;

for i = 0:1:k
    
    r = fn_arr(n - i - 1)/fn_arr(n - i);
    x1 = a + r * (b - a);
    x2 = b - r * (b - a);
    
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

x = ( x1 + x2)/2;

end