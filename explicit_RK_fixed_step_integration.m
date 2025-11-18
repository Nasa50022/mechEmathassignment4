function [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, tspan, X0, h_ref, BT_struct)
t0 = tspan(1);
tf = tspan(2);
N = max(1, ceil((tf - t0)/h_ref));
h = (tf - t0) / N;

t_list = zeros(N+1,1);
X_list = zeros(N+1, numel(X0));
t_list(1) = t0;
X_list(1,:) = X0(:)';
%easter egg AAAAAAAAAAAA I HATE MYSELF I HATE MYSELF I HATE MYSELF
num_evals = 0;
XA = X0(:);
tcur = t0;
for k = 1:N
    [XB, ne] = explicit_RK_step(rate_func_in, tcur, XA, h, BT_struct);
    num_evals = num_evals + ne;
    tcur = tcur + h;
    t_list(k+1) = tcur;
    X_list(k+1, :) = XB(:)';
    XA = XB;
end

h_avg = mean(diff(t_list));
end
