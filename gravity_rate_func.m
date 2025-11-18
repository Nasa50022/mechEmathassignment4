function dVdt = gravity_rate_func(~, V, orbit_params)
x = V(1); y = V(2); vx = V(3); vy = V(4);
r = sqrt(x^2 + y^2);
mu = orbit_params.G * orbit_params.m_sun; 

ax = - mu * x / (r^3);
ay = - mu * y / (r^3);

dVdt = [vx; vy; ax; ay];
end
