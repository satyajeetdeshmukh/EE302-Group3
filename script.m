% Uses Routh Hurwitz method to evaluate the range of the unkown variable
% for a stable, unstable and marginally stable system.
% Initialize
clear; close all; clc

%%Take Inputs

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
% 1) a  c  e+K
% 2) b  d  0
% 3) x  y  0
% 4) z  0  0
% 5) w  0  0

a=1;
b=0;
c=1;
d=0;
e=0;

% K is a symbol as we have to compute its range
syms K
assume(K,{'real', 'positive'})

%% ROW 2

% set flags to 1 if elements in row 2 are zero
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
    syms b
    assume(b,{'real', 'positive'})
    flag_take_limit_b = 1;
else
    flag_take_limit_b = 0;
end
    
% All the elements of second row of the Routh array are zero.
if(flag_b==1 && flag_d==1)
    % Auxiliary equation: a*s^4 + c*s^2 + e + K
    % diff: 4*a*s^3 + 2*c*s
    flag_row2_zero = 1;
    b = 4*a;
    d = 2*c;
else
    flag_row2_zero = 0;
end

%% ROW 3

% using algo for RH table for row 3
x = (b*c-a*d)/b;
y = (e+K);

% set flags to 1 if elements in row 3 are zero
if(x==0)
    flag_x=1;
else
    flag_x=0;
end

% The first element of row 3 of the Routh array is zero.
if(flag_x==1)
    syms x
    assume(x,{'real', 'positive'})
    flag_take_limit_x = 1;
else
    flag_take_limit_x = 0;
end


%% ROW 4

% Computing row 4
z = (x*d-y*b)/x;

if(z==0)
    flag_z=1;
else
    flag_z=0;
end

% All the elements of 4th row of the Routh array are zero.
if(flag_z==1)
    % Auxiliary equation: x*s^2 + y
    % diff: x*s
    z=x;
    flag_row4_zero = 1;
else
    flag_row4_zero = 0;
end

%% ROW 5

w = y;

%% Apply limits

% column 1 of the RH table
col1 = [a b x z w];

if(flag_take_limit_b == 1)
    col1 = limit(col1,b,0,'right');
end

if(flag_take_limit_x == 1)
    col1 = limit(col1,x,0,'right');
end

%% Check if system is unstable independent of K

% The first 3 elements in the col are independent of K, lets see if the
% system is unstable independent of K
unstable_independent = 0;
if (sign(col1(1))*sign(col1(2))==-1 || sign(col1(2))*sign(col1(3))==-1)
    unstable_independent = 1;
    disp('System is unstable independent of K')
end

%% Values of K for marginally stable system
if (unstable_independent == 0)
   
    % ROW 2 ZERO
    if (flag_row2_zero==1)
        % equation is a even polynomial
        % for system to be marginally stable no sign change should happen
        disp(col1)
        eqns = [sign(col1(1)) == sign(col1(2)) sign(col1(2)) == sign(col1(3)) sign(col1(3)) == sign(col1(4)) sign(col1(4)) == sign(col1(5))];
        S = solve(eqns,K,'ReturnConditions',true,'IgnoreAnalyticConstraints',true);
        S.K
        S.parameters
        disp('Equation is a even polynomial and system is marginally stable for K=x :')
        pretty(S.conditions)
    end

    %ROW 4 ZERO
    if (flag_row4_zero==1)
        % Second order even polynomial is a factor
        % for system to be marginally stable no sign change should happen
        disp(col1)
        eqns = [sign(col1(3)) == sign(col1(4)) sign(col1(4)) == sign(col1(5))];
        S = solve(eqns,K,'ReturnConditions',true,'IgnoreAnalyticConstraints',true);
        S.K
        S.parameters
        disp('Second order even polynomial is a factor and System is marginally stable for K=x :')
        pretty(S.conditions)
    end
    
end


%% Find conditions on K for system to be stable.

if (unstable_independent == 0)

    if (flag_row2_zero==0 && flag_row4_zero==0)
        eqns = [sign(col1(3)) == sign(col1(4)) sign(col1(4)) == sign(col1(5))];
        [K, param, cond] = solve(eqns,K,'ReturnConditions',true,'IgnoreAnalyticConstraints',true);
        disp('System is stable for K=x :')
        pretty(cond)
    end

end

    






