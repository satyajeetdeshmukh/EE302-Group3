% Initialize
clear ; close all; clc

% Take in k1, k2, k3, k4
kvector = input('Enter value of [k1 k2 k3 k4] = ');

%% The root polynomial looks like as
% $as^{4} + bs^{3} cs^{2} + ds^{1} + e + K = 0$
%
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






