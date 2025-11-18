function [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in, t, XA, h, BT_struct)
A = BT_struct.A;
B = BT_struct.B;   
C = BT_struct.C;

s = size(A,1);
m = numel(XA);
K = zeros(m, s);
num_evals = 0;

for i = 1:s
    if i==1
        argX = XA;
    else
        argX = XA + h * (K(:,1:i-1) * A(i,1:i-1)');
    end
    ti = t + C(i)*h;
    K(:,i) = rate_func_in(ti, argX);
    num_evals = num_evals + 1;
end

if size(B,1) == 1
    XB1 = XA + h*(K * B(:));
    XB2 = XB1;
else
    XB1 = XA + h*(K * B(1,:).');
    XB2 = XA + h*(K * B(2,:).');
end
