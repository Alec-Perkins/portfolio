
% Clear initials
close all;
clf;
clear variables;
clc;
clear figure;


% Given values:
m1 = 40;
m2 = 20;
k1 = 200;
k2 = 100;
k3 = 250;

% Time values:
t = 0; tFinal = 10.05; dt = 0.05;

alpha = 0.6;
beta = 0.4;

% Equations of amplitude

A = [ (k1 + k2) / m1 , -k2 / m1 ; -k2 / m2 , ( k2 + k3 )/ m2 ];

% Find eigenvalues and eigenvectors

[V, D] = eig(A);

lambda1 = D(1,1);

lambda2 = D(2,2);

omega1 = sqrt(lambda1);
omega2 = sqrt(lambda2);


% Find Amplitudes:

amp1lam1 = V(1,1);
amp2lam1 = V(2,1);

amp1lam2 = V(1,2);
amp2lam2 = V(2,2);

while t < tFinal
    if t+dt > tFinal
        dt = tFinal - t;
    end

   

    % Equations of motion for mode 1:
    mode1x1 = amp1lam1*sin(omega1*t);
    mode1x2 = amp2lam1*sin(omega1*t);

    % Equations of motion for mode 2:
    mode2x1 = amp1lam2*sin(omega2*t);
    mode2x2 = amp2lam2*sin(omega2*t);

    % Combined equations of general motion
    x1 = alpha*amp1lam1*sin(omega1*t) + beta*amp1lam2*sin(omega2*t);
    x2 = alpha*amp2lam1*sin(omega1*t) + beta*amp2lam2*sin(omega2*t);

    

    % Plot motion of each mass in the first mode
    figure(1)
    plot(t, mode1x1, 'b.'); hold on;
    plot(t, mode1x2, 'r.'); 
    xlabel('Time(s)'); 
    ylabel('Position(m)');
    s = sprintf('Mode 1 Positions at Time = %2.2f', t);
    title(s); legend('x1','x2');
    xlim([0 10]); ylim([-1 1]); fixfig; pause(dt);

    % Plot motion of each mass in the second mode
    figure(2)
    plot(t, mode2x1, 'b.'); hold on;
    plot(t, mode2x2, 'r.'); 
    xlabel('Time(s)'); 
    ylabel('Position(m)');
    s = sprintf('Mode 2 Positions at Time = %2.2f', t);
    title(s); legend('x1','x2');
    xlim([0 10]); ylim([-1 1]); fixfig; pause(dt);

    % Plot the combined motion of the modes
    figure(3)
    plot(t, x1, 'b.'); hold on;
    plot(t, x2, 'r.');
    xlabel('Time(s)');
    ylabel('Position(m)');
    s = sprintf('Superposition of Modes at Time = %2.2f (a = %2.1f, b = %2.1f)',t,alpha,beta);
    title(s); legend('x1','x2');
    xlim([0 10]); ylim([-2 2]); fixfig; pause(dt);

    % Prepare for next time step
    t = t + dt;

end























