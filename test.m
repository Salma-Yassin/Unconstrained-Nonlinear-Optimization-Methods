%% Rosenbrock’s parabolic valley function (Fletcher_Reeves)

clear
clc
syms x1 x2

f = 100*(x2 - x1^2)^2 + (1 - x1)^2;
x0 = [1.2 1.0];
eta = 0.01;

tic 
[soln,fmin,k] = fletcher_reeves(x0,f,eta)
toc 


%% Powell’s quartic function (Fletcher_Reeves)

clear
clc

syms x1 x2 x3 x4

f = (x1 + 10 * x2)^2 + 5*(x3 - x4)^2 + (x2 - 2* x3)^4 + 10*(x1 - x4)^4;

x0 = [3.0, -1.0, 0.0, 1.0];

eta = 0.01 ;

tic 
[soln,fmin,k] = fletcher_reeves(x0,f,eta)
toc 

%% Rosenbrock’s parabolic valley function (Marquardt)

clear
clc
syms x1 x2

f = 100*(x2 - x1^2)^2 + (1 - x1)^2;

x0 = [1.2 1.0];

eta = 0.01;
alpha = 10^4;
c1 = 0.25;
c2 =2;
l_option = 1;

tic 
[soln,fmin,k] = marquardt(x0 , f , alpha , c1 , c2 , eta , l_option)
toc 




%% Rosenbrock’s parabolic valley function (Marquardt)

clear
clc

syms x1 x2 x3 x4

f = (x1 + 10 * x2)^2 + 5*(x3 - x4)^2 + (x2 - 2* x3)^4 + 10*(x1 - x4)^4;
x0 = [3.0, -1.0, 0.0, 1.0];

eta = 0.01;
alpha = 10^4;
c1 = 0.25;
c2 = 2;
l_option = 4;


tic 
[soln,fmin,k] = marquardt(x0 , f , alpha , c1 , c2 , eta , l_option)
toc 

%% Rosenbrock’s parabolic valley function (Quashi Newton)
clear
clc
syms x1 x2

f = 100*(x2 - x1^2)^2 + (1 - x1)^2;
x0 = [1.2 1.0];
eta = 0.01;
l_option = 1;

tic 
[soln,fmin, k] = quasi_newton(x0,f,eta,l_option)
toc


%% Powell’s quartic function (Quashi Newton)

clear
clc

syms x1 x2 x3 x4

f = (x1 + 10 * x2)^2 + 5*(x3 - x4)^2 + (x2 - 2* x3)^4 + 10*(x1 - x4)^4;
x0 = [3.0, -1.0, 0.0, 1.0];
eta = 0.01; 
l_option = 1;

tic 
[soln,fmin, k] = quasi_newton(x0,f,eta,l_option)
toc 



%% Fabonacci Section 

clear
clc

syms lambda 

f = 25600 * lambda^4 + 25600 * lambda^3 + 6416 * lambda^2 + 16 * lambda + 4
a = 0
b = 0.1
eta = 0.001

[x] = fibonacci_section(f, a , b , eta)

%% Golden Section 

clear
clc

syms lambda 

f = 25600 * lambda^4 + 25600 * lambda^3 + 6416 * lambda^2 + 16 * lambda + 4
a = 0
b = 0.1
eta = 0.01

[x] = golden_section(f, a , b , eta)


%% quadrtic interpolation 
clear
clc

syms lambda 

f = 25600 * lambda^4 - 25600 * lambda^3 + 6416 * lambda^2 - 16 * lambda + 4
t0 = 0.001
eta = 0.001

[lambda] = quadratic_interpolation(f , t0 , eta)


%% Cubic interpolation 
clear
clc

syms lambda 

f = 25600 * lambda^4 - 25600 * lambda^3 + 6416 * lambda^2 - 16 * lambda + 4
t0 = 0.001
eta = 0.001

[lambda] = cubic_interpolation(f ,t0 , eta)






