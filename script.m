% Uses Routh Hurwitz method to evaluate the range of the unkown variable
% for a stable, unstable and marginally stable system.
% Initialize
clear ; close all; clc

% Take in k1, k2, k3, k4
kvector = input('Enter value of [k1 k2 k3 k4] = ');

% The root polynomial looks like as
% a*s^4 + b*s^3 c*s^2 + d*s^1 + e + K = 0

% Now we will commpute values of these coeff

a=1;
b=0;
c=0;
d=0;
e=1;

for i = 1:4
    b = b + kvector(i);
    e = e*kvector(i);
    for j = i+1:4
        c = c + kvector(i)*kvector(j);
        for k = j+1:4
            d = d + kvector(i)*kvector(j)*kvector(k);
        end
    end
end

% Now we will compute the RH table
% a  c  e+K
% b  d  0
% x  y  0
% z  0  0

% K is a symbol as we have to compute its range
K = sym('K');

% we need to ensure b is not zero else throw an error
if(b==0)
    error('b is zero');
end
    
% using algo for RH table for row 3
x = (b*c-a*d)/b;
y = (e+K);


% using algo for RH table for row 4
if(x==0)
    error('x is zero');
else 
    z = (x*d-y*b)/x;
end

% column 1 of the RH table
col1 = [a b x z];

% column products vector
col_pro = [a*b b*x x*z];

% The first 3 elements in the col are independent of K, lets see if the
% system is unstable independent of K
unstable = 0;
if (sign(col_pro(1)) == -1 || sign(col_pro(2)) == -1)
    unstable = 1;
    disp('System is unstable independent of K')
elseif (z==0)
    disp('System is marginally stable independent of K')
else
    disp('System stability depends on K')
end

% Now we have made sure that that the system's stability depends on K
% To achieve stability col_pro(3) must be the same sign as the other 2
% This  will give us a constraint on K

Kroot = double(solve(col_pro(3), K));

% Due to nature of the problem the equation is linear and we get 3 cases.
% K=Kroot, K > Kroot, and K < Kroot

% if K=Kroot the system is marginally stable
String = ['System  is marginally stable for K=', num2str(Kroot)];
disp(String)

% Lets check for K > Kroot
K = Kroot + 1;
if (sign(col_pro(3)) == -1 && unstable == 0)
    String = ['System  is unstable for K > ', num2str(Kroot), ' and stable for K < ', num2str(Kroot)];
    disp(String)
else
    String = ['System  is stable for K > ', num2str(Kroot), ' and unstable for K < ', num2str(Kroot)];
    disp(String)
end


    






