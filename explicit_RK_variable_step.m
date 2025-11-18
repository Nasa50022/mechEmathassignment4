function [XB, num_evals, h_next, redo] = explicit_RK_variable_step(rate_func_in, t, XA, h, BT_struct, p, error_desired)
safety = 0.9;
alpha = 5;    % maximum growth factor
[XB1, XB2, num_evals] = RK_step_embedded(rate_func_in, t, XA, h, BT_struct);
err_est = norm(XB1 - XB2, inf);
XB = XB1;
if err_est == 0
    factor = alpha;
else
    factor = safety * (error_desired / err_est)^(1/p);
    factor = min(factor, alpha);
end
h_next = factor * h;
redo = (err_est > error_desired);
end
