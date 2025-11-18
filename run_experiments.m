clear; close all; clc;
% does it work sure, is it right? god help me.
G = 6.67430e-11;           
m_sun = 1.98847e30; 
m_planet = 1;
orbit_params.G = G;
orbit_params.m_sun = m_sun;
orbit_params.m_planet = 1;

AU = 1.495978707e11;
    a = AU; e = 0.6;
mu = G * m_sun;
r0 = a*(1-e);
v0 = sqrt( mu * (1+e) / (a*(1-e)) );

V0 = [r0; 0; 0; v0];

T = 2*pi*sqrt(a^3 / mu);
tspan = [0, 0.5*T];

Bogacki = struct();
Bogacki.C = [0, 1/2, 3/4, 1];
Bogacki.A = [0,0,0,0; 1/2,0,0,0; 0,3/4,0,0; 2/9,1/3,4/9,0];
Bogacki.B = [2/9,1/3,4/9,0; 7/24,1/4,1/3,1/8];

DormandPrince = struct();
DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
DormandPrince.A = [0,0,0,0,0,0,0;
1/5,0,0,0,0,0,0;
3/40,9/40,0,0,0,0,0;
44/45,-56/15,32/9,0,0,0,0;
19372/6561,-25360/2187,64448/6561,-212/729,0,0,0;
9017/3168,-355/33,46732/5247,49/176,-5103/18656,0,0;
35/384,0,500/1113,125/192,-2187/6784,11/84,0];
DormandPrince.B = [35/384,0,500/1113,125/192,-2187/6784,11/84,0;
5179/57600,0,7571/16695,393/640,-92097/339200,187/2100,1/40];
%HAHAHAHAAHHAHAHAHA THANK YOU STACK EXCHANGE
HeunEuler = struct();
HeunEuler.C = [0,1];
HeunEuler.A = [0,0;1,0];
HeunEuler.B = [1/2,1/2;1,0];

 BT = DormandPrince; 
p = 5; 

t_ref = linspace(tspan(1), tspan(2), 500);
Vref = compute_planetary_motion(t_ref, V0, orbit_params);

interp_ref = @(tq) interp1(t_ref, Vref, tq, 'spline');

h_list = [T/200, T/400, T/800, T/1600, T/3200];
errors_fixed = zeros(size(h_list));
evals_fixed = zeros(size(h_list));
for i = 1:numel(h_list)
    h = h_list(i);
    [tlist, Xlist, havg, ne] = explicit_RK_fixed_step_integration(@(tt,XX) gravity_rate_func(tt,XX,orbit_params), tspan, V0, h, BT);
    V_end = Xlist(end,:)';
    V_truth_end = interp_ref(tlist(end))';
    errors_fixed(i) = norm(V_end - V_truth_end, inf);
    evals_fixed(i) = ne;
    fprintf('Fixed h=%.3e: error=%.3e, evals=%d, steps=%d\n', h, errors_fixed(i), ne, numel(tlist)-1);
end

error_desired_list = [1e6, 1e5, 1e4, 1e3, 1e2]; 
errors_adaptive = zeros(size(error_desired_list));
evals_adaptive = zeros(size(error_desired_list));
havg_adaptive = zeros(size(error_desired_list));
failure_rates = zeros(size(error_desired_list));
for i = 1:numel(error_desired_list)
    ed = error_desired_list(i);
    h_ref = T/200;
    [tlistA, XlistA, havg, ne, failrate] = explicit_RK_variable_step_integration(@(tt,XX) gravity_rate_func(tt,XX,orbit_params), tspan, V0, h_ref, BT, p, ed);
    V_end = XlistA(end,:)';
    V_truth_end = interp_ref(tlistA(end))';
    errors_adaptive(i) = norm(V_end - V_truth_end, inf);
    evals_adaptive(i) = ne;
    havg_adaptive(i) = havg;
    failure_rates(i) = failrate;
    fprintf('Adaptive error_des=%g: global_err=%.3e, evals=%d, havg=%.3e, failrate=%.3f\n', ed, errors_adaptive(i), ne, havg, failrate);
end

figure; hold on; grid on;
loglog(h_list, errors_fixed, 'o-','DisplayName','Fixed RK');
loglog(havg_adaptive, errors_adaptive, 's-','DisplayName','Adaptive RK');
xlabel('Average step size (s)'); ylabel('Global error (inf-norm)');
title('Global error vs avg step size (fixed vs adaptive)');
legend('Location','best');
%easter egg https://tenor.com/cmOuxNjrpEk.gif
figure; hold on; grid on;
loglog(evals_fixed, errors_fixed, 'o-','DisplayName','Fixed RK');
loglog(evals_adaptive, errors_adaptive, 's-','DisplayName','Adaptive RK');
xlabel('Function evaluations'); ylabel('Global error (inf-norm)');
title('Global error vs function evaluations');

figure; semilogx(havg_adaptive, failure_rates, 'o-'); grid on;
xlabel('Average step size (s)'); ylabel('Step failure rate'); title('Failure rate vs avg step size (adaptive)');
fprintf('\nFixedstep summary:\n');
T1 = table(h_list', errors_fixed', evals_fixed', 'VariableNames', {'h','global_error','num_evals'});
disp(T1);

fprintf('\nAdaptive summary:\n');
T2 = table(error_desired_list', havg_adaptive', errors_adaptive', evals_adaptive', failure_rates', ...
    'VariableNames', {'error_desired','avg_h','global_error','num_evals','failure_rate'});
disp(T2);


figure;
subplot(2,1,1);
plot(tlistA, XlistA(:,1)/AU, 'b-', tlistA, XlistA(:,2)/AU, 'r-');
xlabel('Time (s)'); ylabel('Position (AU)');
legend('x','y'); grid on;
title('Position vs time (adaptive)');

subplot(2,1,2);
plot(tlistA, XlistA(:,3)/1e4, 'b-', tlistA, XlistA(:,4)/1e4, 'r-');
xlabel('Time (s)'); ylabel('Velocity (×10^4 m/s)');
legend('vx','vy'); grid on;
title('Velocity vs time (adaptive)');


r = sqrt(XlistA(:,1).^2 + XlistA(:,2).^2);
h_list = diff(tlistA);
 figure;
semilogy(r(1:end-1)/AU, h_list, 'o-');
xlabel('Distance from sun (AU)');
ylabel('Step size (s)');
title('Stepsize clustering vs distance');
grid on;



x = XlistA(:,1); y = XlistA(:,2);
vx = XlistA(:,3); vy = XlistA(:,4);
r = sqrt(x.^2 + y.^2);
v2 = vx.^2 + vy.^2;
mu = orbit_params.G * orbit_params.m_sun;
E = 0.5*v2 - mu./r;            
 H = x.*vy - y.*vx;             

figure;
subplot(2,1,1);
plot(tlistA, (E - E(1))/abs(E(1)));
xlabel('Time (s)'); ylabel('\DeltaE / E_0');
title('Relative change in energy'); grid on;

subplot(2,1,2);
plot(tlistA, (H - H(1))/abs(H(1)));
xlabel('Time (s)'); ylabel('\DeltaH / H_0');
title('Relative change in angular momentum'); grid on;



test_BT = Bogacki; 
test_fun = @(t,x) -x;
x0 = 1;
t0 = 0;
x_true = @(t) exp(-t);

h_vals = logspace(-4, -1, 8);
lte_vals = zeros(size(h_vals));
for i = 1:length(h_vals)
    h = h_vals(i);
    [xB, ne] = explicit_RK_step(test_fun, t0, x0, h, test_BT);
    lte_vals(i) = abs(x_true(h) - xB);
end

figure; loglog(h_vals, lte_vals, 'o-');
xlabel('Step size h'); ylabel('Local truncation error');
title('Local truncation error scaling (x''=-x)');
grid on;


p_est = polyfit(log(h_vals), log(lte_vals), 1);
fprintf('\nEstimated local truncation error order (slope): %.2f\n', -p_est(1));
