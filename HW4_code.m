clear all; close all; clc;

%% Useful Equations

% r_ddot = -mu/r^3 * r + a_d;
% rdot_hat = rdot ir + r fdot it
% rdot = r^2 fdot e sin(f) / p


%% Initial Orbital Elements

a = 7300; % km
e = 0.05;
i = 42; % degrees
Omega = 0; % degrees
omega = 45; % degrees
M0 = 0; % degrees


%% Initial Conditions

G = 6.674 * 10^(-20); % km^3/kg-s^2
m1 = 5.97219 * 10^(24); % kg
mu = G * m1;

% Conversion matrix to Cartesian
theta0 = omega;
ON = [cosd(theta0) sind(theta0)*cosd(i) sind(theta0)*sind(i);...
    -sind(theta0) cosd(theta0)*cosd(i) cosd(theta0)*sind(i);...
    0 -sind(i) cosd(i)];
NO = ON';

% Position IC's
p = a*(1-e^2);
r_O0 = [p/(1+e);0;0];
r_N0 = NO * r_O0;

% Velocity IC's
f0 = 0;
r = r_O0(1);
f_dot = sqrt(mu * p) / r^2;
v_O0 = [0;r*f_dot;0];
v_N0 = NO * v_O0;


%% Implementation

tfinal = 10 * (2*pi*sqrt(a^3/mu));
tspan = [0 tfinal];
x0 = [r_N0;v_N0];

[t,x] = ode89(@vdp1,tspan,x0);

r_earth = 6378; % km
[X,Y,Z] = sphere;
X2 = X * r_earth;
Y2 = Y * r_earth;
Z2 = Z * r_earth;

figure(1)
hold on
plot3(x(:,1),x(:,2),x(:,3))
surf(X2,Y2,Z2,'FaceAlpha',0.1,'FaceColor','c')
xlim([-8000 8000])
ylim([-8000 8000])
zlim([-8000 8000])
view(35,10)
hold off
title('10-Orbit Non-Linear Numerical Simulation')
xlabel('X_{EFI} (km)')
ylabel('Y_{ECI} (km)')
zlabel('Z_{ECI} (km)')


%% Data Conversion to Classical Orbital Elements

h = cross(r_N0,v_N0);
e_i = cross(v_N0,h);

a = [];
e = [];
Omega = [];
inc = [];
omega = [];
M0 = [];
for i = 1:length(x)
    
    % Semi-major axis, a
    r_i = sqrt(x(i,1)^2 + x(i,2)^2 + x(i,3)^2);
    v_i = sqrt(x(i,4)^2 + x(i,5)^2 + x(i,6)^2);
    a(i) = (2/r_i - v_i^2/mu)^(-1);

    % Eccentricity, e
    r_t = [x(i,1);x(i,2);x(i,3)];
    v_t = [x(i,4);x(i,5);x(i,6)];
    h_vec = cross(r_t,v_t);
    e_vec = cross(v_t,h_vec)/mu - r_t/r_i;
    e(i) = norm(e_vec);

    % Angles, Omega omega and i, from DCM
    i_e = e_vec/e(i);
    i_h = h_vec/norm(h_vec);
    i_p = cross(i_h,i_e);
    DCM = [i_e(1) i_e(2) i_e(3);...
        i_p(1) i_p(2) i_p(3);...
        i_h(1) i_h(2) i_h(3)];
    Omega(i) = atand(-DCM(3,1)/DCM(3,2));
    inc(i) = acosd(DCM(3,3));
    omega(i) = atand(DCM(1,3)/DCM(2,3));

    % Mean anomaly
    M0(i) = acos(dot(e_vec,e_i)/(norm(e_vec)*norm(e_i)));

end


%% Comparison of numerical and analytical techniques for e(t)

a_init = 7300; % km
e_init = 0.05;
i_init = 42; % degrees
Omega_init = 0; % degrees
omega_init = 45; % degrees
M0_init = 0; % degrees

e0 = [Omega_init;i_init;omega_init;a_init;e_init;M0_init];
[t_a,e_bold] = ode89(@analytical,tspan,e0);
[t_bar,e_bar] = ode89(@ebar,tspan,e0);

% Plot Omega, i, and omega
figure(2)
tiledlayout(3,1)
nexttile
hold on
plot(t,Omega)
plot(t_a,real(e_bold(:,1)))
plot(t_bar,e_bar(:,1))
title('Longitude of Ascending Node')
ylabel('(^\circ)')
legend('Numerical','Analytical','Analytical Mean')
nexttile
hold on
plot(t,inc)
plot(t_a,real(e_bold(:,2)))
plot(t_bar,e_bar(:,2))
title('Inclination')
ylabel('(^\circ)')
nexttile
hold on
plot(t,omega)
plot(t_a,real(e_bold(:,3)))
plot(t_bar,e_bar(:,3))
title('Argument of Periapsis')
ylabel('(^\circ)')
xlabel('Time (seconds)')

% Plot a, e, and M0
figure(3)
tiledlayout(3,1)
nexttile
hold on
plot(t,a)
plot(t_a,real(e_bold(:,4)))
plot(t_bar,e_bar(:,4))
title('Semi-Major Axis')
ylabel('(km)')
legend('Numerical','Analytical','Analytical Mean')
nexttile
hold on
plot(t,e)
plot(t_a,real(e_bold(:,5)))
plot(t_bar,e_bar(:,5))
title('Eccentricity')
nexttile
hold on
plot(t,M0)
plot(t_a,real(e_bold(:,6)))
plot(t_bar,e_bar(:,6))
title('Mean Anomaly')
ylabel('(^\circ)')
xlabel('Time (seconds)')


%% Functions

function dxdt = vdp1(t,x)

    J2 = 1082.63*10^(-6);
    
    G = 6.674 * 10^(-20); % km^3/kg-s^2
    m1 = 5.97219 * 10^(24); % kg
    mu = G * m1;
    
    r_x = x(1); r_y = x(2); r_z = x(3);
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    r_eq = 6371; % km
    
    a_J2_coeff = -3/2*J2*(mu/r^2)*(r_eq/r)^2;
    a_J2_comps = [(1-5*(x(3)/r)^2)*(x(1)/r);...
        (1-5*(x(3)/r)^2)*(x(2)/r);(3-5*(x(3)/r)^2)*(x(3)/r)];
    a_J2 = a_J2_coeff * a_J2_comps;
    
    v1 = x(4);
    v2 = x(5);
    v3 = x(6);
    
    dxdt = [v1;v2;v3;-mu/r^3*r_x+a_J2(1);...
        -mu/r^3*r_y+a_J2(2);...
        -mu/r^3*r_z+a_J2(3)];

end

function debolddt = analytical(t,x)

    G = 6.674 * 10^(-20); % km^3/kg-s^2
    m1 = 5.97219 * 10^(24); % kg
    mu = G * m1;

    Omega = x(1);
    i = x(2);
    omega = x(3);
    a = x(4);
    e = x(5);
    M0 = x(6);

    J2 = 1082.63*10^(-6);
    r_eq = 6371; % km

    n = sqrt(mu/a^3);
    b = a*sqrt(1-e^2);

    p = a*(1-e^2);
    eta = sqrt(1-e^2);

    term1 = M0 + n*t;
    E_pre = 0:0.01:100;
    tolerance = 10^(-2);
    for j = 1:length(E_pre)
        term2 = E_pre(j) - e*sin(E_pre(j));
        if abs(term1-term2) <= tolerance
            E = E_pre(j);
            break
        end
    end

    f = rad2deg(2*atan(sqrt((1+e)/(1-e))*tan(E/2)));
    theta = omega + f;
    r = a * (1 - e*cos(f));

    dOdt = -3*J2*n* a^2/(b*r) * (r_eq/r)^2 *...
        sind(theta)^2 * cosd(i) * 180/pi;
    didt = -3/4*J2*n* a^2/(b*r) * (r_eq/r)^2 *...
        sind(2*theta) * sind(2*i) * 180/pi;
    dodt = 3/2*J2*n* p/(r^2*e*eta^3) * (r_eq/r)^2 *...
        (2*r*e*cosd(i)^2*sind(theta)^2 - ...
        (p+r)*sind(f)*sind(i)^2*sind(2*theta) + ...
        p*cosd(f)*(1-3*sin(i)^2*sind(theta)^2)) * 180/pi;
    dadt = -3*J2*n * a^4/(b*r^2) * (r_eq/r)^2 *...
        (e*sind(f)*(1-3*sind(theta)^2*sind(i)^2) +...
        p/r*sind(2*theta)*sind(i)^2);
    dedt = -3/2*J2*n * a^2/(b*r) * (r_eq/r)^2 *...
        (p/r*sind(f)*(1-3*sind(theta)^2*sind(i)^2) +...
        (e+cosd(f)*(2+e*cosd(f)))*sind(2*theta)*sind(i)^2);
    dM0dt = 3/2*J2*n * p/(r^2*e*eta^2) * (r_eq/r)^2 *...
        ((p+r)*sind(f)*sind(i)^2*sind(2*theta) +...
        (2*r*e-p*cosd(f))*(1-3*sind(i)^2*sind(theta)^2));

    debolddt = [dOdt;didt;dodt;dadt;dedt;dM0dt];

end

function debolddt = ebar(t,x)

    G = 6.674 * 10^(-20); % km^3/kg-s^2
    m1 = 5.97219 * 10^(24); % kg
    mu = G * m1;

    Omega = x(1);
    i = x(2);
    omega = x(3);
    a = x(4);
    e = x(5);
    M0 = x(6);

    J2 = 1082.63*10^(-6);
    r_eq = 6371; % km

    n = sqrt(mu/a^3);
    p = a*(1-e^2);

    dOdt = -3/2*J2*n * (r_eq/p)^2 * cosd(i) * 180/pi;
    didt = 0;
    dodt = 3/4*J2*n * (r_eq/p)^2 * (5*cosd(i)^2 - 1) * 180/pi;
    dadt = 0;
    dedt = 0;
    dM0dt = 3/4*J2*n * (r_eq/p)^2 * sqrt(1-e^2) * (3*cosd(i)^2 - 1);

    debolddt = [dOdt;didt;dodt;dadt;dedt;dM0dt];

end
