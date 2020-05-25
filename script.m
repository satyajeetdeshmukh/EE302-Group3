% Initialize
clear ; close all; clc

% Take in k1, k2, k3, k4
kvector = input('Enter value of [k1 k2 k3 k4] = ');

% The root polynomial looks like as
% $as^{4} + bs^{3} cs^{2} + ds^{1} + e + K = 0$

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

% using algo for RH table
x = (b*c-a*d)/b;
y = (e+K);
z = (x*d-y*b)/x;


% column 1 of the RH table
% col1 = [a b x z];

% column products vector
col_pro = [a*b b*x x*z];

% only 3rd product depends upon K. We can check if we got a unstable pole
% from the first 2
unstable = 0;
if (sign(col_pro(1)) == -1 || sign(col_pro(2)) == -1)
    unstable = 1;
    disp('System is unstable')
else
    disp('Needs further analysis with K')
end








