function [t_list,X_list,h_avg, num_evals, step_failure_rate] = explicit_RK_variable_step_integration(rate_func_in, tspan, X0, h_ref, BT_struct, p, error_desired)
t0 = tspan(1); tf = tspan(2);
tcur = t0;
XA = X0(:);
h = min(h_ref, tf - tcur);
t_list = tcur;
X_list = XA.';
num_evals = 0;
attempted = 0;
failed = 0;

while tcur < tf - 1e-14
    h = min(h, tf - tcur);
    attempted = attempted + 1;
    [XB, ne, h_next, redo] = explicit_RK_variable_step(rate_func_in, tcur, XA, h, BT_struct, p, error_desired);
    num_evals = num_evals + ne;
    if redo
        failed = failed + 1;
        h = h_next;
        if h < 1e-14
            warning('your code it broke you are stupid, I too am in this error message');
            break;
        end
        continue
    else
        tcur = tcur + h;
        XA = XB;
        t_list(end+1,1) = tcur;
        X_list(end+1,:) = XA.';
        h = h_next;
    end
end

h_avg = mean(diff(t_list));
step_failure_rate = failed / max(1, attempted);
end
