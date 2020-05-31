% Uses Routh Hurwitz method to evaluate the range of the unkown variable
% for a stable, unstable and marginally stable system.
% Initialize
clear; close all; clc

% % Take in k1, k2, k3, k4
% kvector = input('Enter value of [k1 k2 k3 k4] = ');
% 
% % The root polynomial looks like as
% % a*s^4 + b*s^3 c*s^2 + d*s^1 + e + K = 0
% 
% % Now we will commpute values of these coeff
% 
% a=1;
% b=0;
% c=0;
% d=0;
% e=1;
% 
% for i = 1:4
%     b = b + kvector(i);
%     e = e*kvector(i);
%     for j = i+1:4
%         c = c + kvector(i)*kvector(j);
%         for k = j+1:4
%             d = d + kvector(i)*kvector(j)*kvector(k);
%         end
%     end
% end

% Now we will compute the RH table
% a  c  e+K
% b  d  0
% x  y  0
% z  0  0
% w  0  0

a=1;
b=0;
c=6;
d=4;
e=1;

% K is a symbol as we have to compute its range
syms K
assume(K,{'real', 'positive'})

if(b==0)
    flag_b=1;
else
    flag_b=0;
end

if(d==0)
    flag_d=1;
else
    flag_d=0;
end


% The first element of any row of the Routh array is zero.
if(flag_b==1 && flag_d==0)
    disp('b is zero, d is non-zero')
    syms b
    assume(K,{'real', 'positive'})
    flag_1 = 1;
else
    flag_1 = 0;
end
    
% All the elements of any row of the Routh array are zero.
if(flag_b==1 && flag_d==1)
    % Auxiliary equation: a*s^4 + c*s^2 + e + K
    % diff: 4*a*s^3 + 2*c*s
    b = 4*a;
    d = 2*c;
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

% using algo for RH table for row 5
if(z==0)
    error('z is zero');
else 
    w = y;
end

% column 1 of the RH table
col1 = [a b x z w];

if(flag_1 == 1)
    col1 = limit(col1,b,0,'right');
end

% The first 3 elements in the col are independent of K, lets see if the
% system is unstable independent of K
unstable = 0;
if (sign(col1(1)) ~= sign(col1(2)) || sign(col1(2)) ~= sign(col1(3)))
    unstable_poles = 0;
    for i = 1:4
        if(sign(col1(i)) ~= sign(col1(i+1)))
            unstable_poles = unstable_poles + 1;
        end
    end
     
    error('System is unstable independent of K')
else
    eqns = [sign(col1(3)) == sign(col1(4)) sign(col1(4)) == sign(col1(5))];
    S = solve(eqns,K,'ReturnConditions',true,'IgnoreAnalyticConstraints',true);
    S.K
    S.parameters
    disp('System is stable for K=x :')
    pretty(S.conditions)

end


    






