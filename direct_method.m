% Initialize
clear ; close all; clc
% 
% % Take in k1, k2, k3, k4
% kvector = input('Enter value of [k1 k2 k3 k4] = ');
% 
% % The root polynomial looks like as a*s^4 + b*s^3 c*s^2 + d*s^1 + e + K = 0
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

syms K
assume(K,{'real', 'positive'})
coefvct = [1 1 3 K]; 
x = roots(coefvct);     
x_real = real(x);
eqn1 = K > 0;
eqn2 = imag(K) == 0;
eqn3 = x_real<0;
eqns = [eqn3(1) eqn3(2) eqn3(3)];
S = solve(eqns,K,'ReturnConditions',true,'IgnoreAnalyticConstraints',true,'MaxDegree', 3);
S.K
S.parameters
pretty(S.conditions)



