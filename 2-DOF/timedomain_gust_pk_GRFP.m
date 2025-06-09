  %% Computing forces from Theodorsens equations
close all
clc;
clear;
format long

% %Establishing data from the system structure in a structure (dimensional data)
%Span needs to be larger than the chord for the theory to apply (analogous
%to lifting line)

%Parameters included here
%c --> Non-dimensional distance from midchord to aileron hinge NOT INCLUDED IF JUST 2DOF
%b --> Semichord of the wing (used to non-dimensionalize all distances)
%a --> Non-dimensional distance from midchord to elastic axis of the wing
% x_alpha --> Non-dimensional distance from axis of rotation (elastic axis) to the center
% of gravity of the airfoil
%EI --> Bending rigidity of the wing
%GJ --> Torsional rigidity
%s --> semispan of the wing 
%M --> Mass per unit area of the wing

%Data from Wright and cooper to validate against it (value of a has been
%adapted to match that data, since it is given in terms of xf)
Struct = struct('b', 1, 's', 7.5, 'a', -0.2, 'EI', 4*10^7, 'GJ', 8*10^6, 'x_alpha', 0.04, 'M', 400);

%Parametrizing the structure with N points to integrate along the span
N = 500;
dy = Struct.s/N;
dx = Struct.b/N;
y = 0:dy:Struct.s;

%Parametrizing the chord of the structure
x = [0:dx:2*Struct.b];

%Creating a grid for double numerical integration
[X,Y] = meshgrid(x,y);

% %Bending mode shape (grid for double integration, vector for single
% integration)
 bending_deflection = (y/Struct.s).^2;
 bending = (Y/Struct.s).^2;

 %Wagner's model
 Aero = struct('eps1',0.0455,'eps2',0.3, 'C1', 0.165,'C2',0.335,'rho',1.225);
 Aero.phi_0 = 1-Aero.C1-Aero.C2;


% %Torsional mode shape
torsion = (Y/Struct.s).*(X-(Struct.b+Struct.b*Struct.a)); %We need the distance wrt elastic axis for the structure
torsion_deflection = (y/Struct.s); %For the aerodynamic loads we only need the deflection
R = R_parameters(bending_deflection,torsion_deflection,y);

%Computing structural matrices (assuming no structural damping)
[A_struct, inv_A_struct,  C_struct, omega_nat] = matrices_struct(Struct,bending,bending_deflection,torsion,torsion_deflection,y,x); %frequency given in rad/s


%Establishing the velocity range for the reference frequency values
V=linspace(1,300,300);

%Velocity points for which system identification will be used
Velocity_Time = [250, 252, 254];


%Preallocating the output variables
omega= zeros(length(V),length(A_struct));
damp = zeros(length(V),length(A_struct));

for i = 1:length(V)
    [Ma,Ca,Ka,W, Evol_aero] = AerodynamicMatricesTime(Struct,Aero,V(i),R);
    [SM] = SystemMatrix(A_struct,C_struct,Ma,Ca,Ka,W,Evol_aero);
    eval = eig(SM);
    index=1;
    for j = 1:length(eval)
        if imag(eval(j))>0
            lambdas(index) = eval(j);
            index = index+1;
        end
    end
    omega(i,:) = imag(lambdas);
    damp(i,:) = real(lambdas(:))./abs(lambdas(:));

end

omega = omega/(2*pi);
damp = damp*100;

load("FreqResults.mat")
figure(1)
hold on
title("System in time and frequency domain")
subplot(2,1,1)
hold on
h2 = plot(V_freq, omega_freq(:,1), "Color","r",LineWidth=2);
hold on
h1 = plot(V,omega(:,1), "k--", LineWidth=2);
hold on

 plot(V_freq, omega_freq(:,2), "Color","r",LineWidth=2);
  plot(V,omega(:,2),"k--",LineWidth=2);
 hold on
% plot(V,omega(:,3),'b');

% xlabel('$V/(\omega_\alpha \, b)$', 'Interpreter','latex')
xlabel('$V$ [m/s]', 'Interpreter','latex')
legend('Time domain', 'Frequency domain', 'Interpreter','latex');
xlim([0 V(end)])
% ylim ([0 8])
grid on
grid minor
ax = gca;
% ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '-';
ax.GridAlpha = 0.5;
% ax.Layer = 'top';
ylabel('$\omega$ [Hz]', 'Interpreter','latex')
pbaspect([3 1 1])
% title("System in time domain")

hold on
subplot(2,1,2)
hold on
plot(V_freq, damp_freq(:,1), "Color",  "red",LineWidth=2);
hold on
plot(V, damp(:,1), 'k--',LineWidth=2);
hold on
plot(V_freq, damp_freq(:,2), "Color",  "red",LineWidth=2);
hold on
plot(V, damp(:,2), 'k--',LineWidth=2);

% plot(V, damp(:,3), 'b');

% ylim([-20, 5])
xlim([0,V(end)])
grid on
grid minor
ax = gca;
% ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '-';
ax.GridAlpha = 0.5;
% ax.Layer = 'top';
yline(0, 'Color', 'k', 'LineWidth', 2); % Draw line for X axis.
xlabel('$V$ [m/s]', 'Interpreter','latex');
ylabel ('$\xi$ [\%]', 'Interpreter','latex');
pbaspect([3 1 1])


%Computing the time evolution of the system
%DDefining a time vector

% A VERY LONG TIME VECTOR GREATLY INCREASES COMPUTATIONAL TIME
fsample = 50;
dt = 1/fsample;
time = 1:dt:200;

%Defining the initial values for the DOFs:
alphadot0 = 0;
hdot0 = 0;
alpha0 = 0;
h0 = 0.5;
lambda1_0 = 0;
lambda2_0 = 0;

%Arranging them in a vector
x0 = [alphadot0; hdot0; alpha0; h0; lambda1_0; lambda2_0];

%Preallocating output vectors
x = zeros(length(x0), length(time));
x(:,1) = x0;

for i =1:length(Velocity_Time)
    V_g =Velocity_Time(i);
    [Ma,Ca,Ka,W,Evol_aero] = AerodynamicMatricesTime(Struct, Aero,V_g,R);
    [SM, inv_Mae] = SystemMatrix(A_struct,C_struct,Ma,Ca,Ka,W,Evol_aero);
    
    %time evolution of the system with ode45
    [t,x] = ode45(@(t,x)TimeEvolution(t,SM,x),time,x0);
    
    % %Computing the time evolution of the system with ode45
    % for i =2:length(time)
    %     x(:,i) = TimeEvolution(t,SM, x(:,i-1));
    % end
    
    
    
    %Amplitude of the angle for the gust
    Amplitude = 0.02; % given in radians
    
    %Generating a gust
    theta = GustAngle(time, Amplitude,lambda1_0, lambda2_0, Evol_aero);
    FM_gust = Gust(Ma,Ca,Ka,W,Aero,theta);
    
    %Computing the time evolution of the system in presence of a gust
    [t,x_gust] = ode45(@(t,x)GustyTimeEvolution(t,SM,x,FM_gust,inv_Mae,time), time, x0);
    
    % theta = theta(3,:);
    % 
    % save('MySystem.mat',"theta   ","x_gust");
    % figure(20)
    % hold on
    % plot(t,FM_gust(1,:))
    % title ('Aerodynamic moment due to a random gust with $\theta_g =0.02$ rad and $U_\infty =100$ m/s', 'Interpreter', 'latex');
    % xlabel('time (s)')
    % ylabel ('$M_g$ [N·m]', Interpreter='latex')
    % 
    % figure(21)
    % hold on
    % plot(t,FM_gust(2,:))
    % plot(t,FM_gust(1,:))
    % title ('Aerodynamic force due to a random gust with $\theta_g =0.02$ rad and $U_\infty =100$ m/s', 'Interpreter', 'latex');
    % xlabel('time (s)')
    % ylabel ('$F_g$ [N]', Interpreter='latex')

    % FIGURES FOR DEBUGGING PURPOSES
    % figure(22)
    % hold on
    % plot(t,theta(3,:))
    % title ('Random gust angle given as an input to the system', 'Interpreter', 'latex');
    % xlabel('time (s)', Interpreter='latex')
    % ylabel ('$\theta_g$ [rad]', Interpreter='latex')
    % grid on
    % 
    % figure(23)
    % hold on
    % plot(t,x_gust(:,3))
    % title ('Angle of attack of the system $\alpha(t)$', 'Interpreter', 'latex');
    % xlabel('time (s)', Interpreter='latex')
    % ylabel ('$\alpha$ [rad]', Interpreter='latex')
    % grid on
    % 
    % figure(24)
    % hold on
    % plot(t,x_gust(:,4))
    % title ('Plunge of the wing $h(t)$', 'Interpreter', 'latex');
    % xlabel('time (s)', Interpreter='latex');
    % ylabel ('$h$ [m]', Interpreter='latex')
    % grid on
    
    % %Numerator and denominators of the transfer functions involved
    % [Hnum,Hden] = TransferFunction(A_struct, C_struct,Ma,Struct,Aero,V_g,R);
    % 
    % 
    displacements(i,:,:) = [x_gust(:,3), x_gust(:,4)];

    random_force_signal(i,:) = theta(3,:).';
end
numerator_order = 2;
denominator_order = 4;
save('SIDFiles.mat',"displacements","random_force_signal", ...
    "numerator_order","denominator_order", "omega", "damp","Velocity_Time");

%% Aeroelasticity and time simulation related functions
function [A_struct, inv_A_struct, C_struct, omega_nat] = matrices_struct(Struct, bending,bending_deflection,torsion,torsion_deflection,y,x)

%Returns the structural matrices of the system following Wright and
%Cooper's approach
    EI = Struct.EI;
    GJ = Struct.GJ;
    M = Struct.M;

    %Computing the coefficients of the structure integrating the mode
    %shapes

    %From the kinetic energy
    %Contribution to the inertia matrix from bending diff in kinetic energy
    A_struct11 = M*trapz(y,trapz(x,torsion.^2,2)); %Double integration
    A_struct13 = M*trapz(y,trapz(x,bending.*torsion,2)); %Contribution to combination of bending and torsion
    A_struct31 = A_struct13;
    A_struct33 = M* trapz(y,trapz(x, bending.^2,2) ); %Contribution to q_bdotdot

    A_struct = [A_struct11, A_struct13;
                A_struct31, A_struct33];
    
    %Its inverse
    inv_A_struct = inv(A_struct);
     
    %No structural damping being considered

    %Computing the mode shapes from the potential energy for the stiffness
    %matrix
    torsion_strain = gradient(torsion_deflection,y); %Diff uses right hand differences, so using gradient instead
    bending_strain = gradient(gradient(bending_deflection,y),y);
    C_33 = EI*trapz(y,(bending_strain).^2);
    C_11 = GJ*trapz(y,(torsion_strain).^2);

    %Structural stiffness matrix
    C_struct = [C_11, 0;
                0, C_33];

    %Finding the natural frequencies of the system
    lam = eig(inv_A_struct*C_struct);
    omega_nat = sqrt(lam);
end

function [R] = R_parameters(bending, torsion,y)

%This function obtains the R coefficients to multiply times the Q matrix
%and obtain the full flutter equations for the case of a finite wing with
%NO aileron applying strip theory

%There are empty coefficients to ensure a consistent defintion of the
%matrix throughout different codes

R = zeros(3);
%Integrals that we need to compute to find the coefficients in the
%aerodynamic loads matrix
int1 = torsion.^2;
R(1,1) = trapz(y, int1); %R11 in the matrix

int3 = bending.*torsion;
R(1,3) = trapz(y,int3); %Both R13 and R31
R(3,1) = trapz(y,int3); 

int6 = bending.^2;
R(3,3) = trapz(y, int6); %R33

end

function [Ma,Ca,Ka,W,Evol_aero] = AerodynamicMatricesTime(Struct,Aero,V,R)
%Computes aerodynamic loads in time domain with R.T. Jones approximation
    a = Struct.a;
    b= Struct.b;

    rho = Aero.rho;
    eps1 = Aero.eps1;
    eps2 = Aero.eps2;
    phi_0 = Aero.phi_0;
    C1 = Aero.C1;
    C2 = Aero.C2;

    %Aerodynamic inertia matrix (no circulatory part
    Ma = pi*rho*b^2*[-b*(1/8 + a^2)*R(1,1), b*a*R(1,3);
                     b*a*R(3,1), -1*R(3,3)];

    %Aerodynamic  damping matrix
    Ca_non = pi*rho*b*[-V*b*(0.5 -a)*R(1,1), 0;
                     -b*V*R(3,1),             0];
    Ca_circ = 2*pi*rho*V*b*phi_0*[ -b*(0.5+a)*b*(0.5-a)*R(1,1), b*(0.5+a)*R(1,3);
                                   -b*(0.5-a)*R(3,1), -1*R(3,3)];
    Ca = Ca_non+ Ca_circ;

    %Aerodynamic stiffness matrix
    Ka = 2*pi*rho*V*b*phi_0*[b*(0.5+a)*V*R(1,1), 0;
                                  -V*R(3,1),0];
    
    %Aerodynamic state matrix
    W = 2*pi*rho*V*b*[b*(0.5+a)*R(1,1), b*(0.5+a)*R(1,3);
                        -1*R(3,1),      -1*R(3,3)];
    
    %Computing the evolution of the aerodynamic states
    Evol11 = V/b*[ C1*eps1*b*(0.5-a), C1*eps1; 
                   C2*eps2*b*(0.5-a), C2*eps2];
    Evol12 = V/b*[C1*eps1*V, 0;
                  C2*eps2*V, 0];
    Evol13 = V/b*[-eps1, 0;
                  0, -eps2];
    %Arranging them in a matrix
    Evol_aero = [Evol11, Evol12, Evol13];
end

function [SM, inv_Mae] = SystemMatrix(A_struct,C_struct,Ma,Ca,Ka,W,Evol_aero)
    %Computing the matrix containing the eignvalues of the system
    Mae = A_struct-Ma;
    Kae = C_struct-Ka;
    inv_Mae = inv(Mae);

    SM11 = inv_Mae*Ca;
    SM12 = -inv_Mae*Kae;
    SM13 = inv_Mae*W;

    SM = [SM11, SM12, SM13;
          eye(2), zeros(2), zeros(2);
          Evol_aero];

end

function [xdot] = TimeEvolution(t,SM, x0)
    %x0 is a state vector containing the degrees of freedom at the
    %beginning of the timestep

    %Order of the variables passing through x0
    % alphadot0 = x0(1);
    % hdot0 = x0(2);
    % alpha = x0(3);
    % h = x0(4);
    % lambda1_0 = x0(5);
    % lambda2_0 = x0(6);

    %The evolution of these matrices is given by the system's matrix:
    xdot = SM*x0;
    
end

function [FM] = Gust(Ma,Ca,Ka,W,Aero,theta)
    %Takes as an input the time at which we are computing the gust, the
    %mean angle theta about which we are computing the gust and the
    %amplitude of the gust in radians

    rho = Aero.rho;
    % eps1 = Aero.eps1;
    % eps2 = Aero.eps2;
    phi_0 = Aero.phi_0;
    % C1 = Aero.C1;
    % C2 = Aero.C2;

    %Given that angle theta, computing its force and moment as we did for
    %aerodynamic forces, matrices need reshaping to get rid of plunge
    %dependent terms
    %Force and moment
        %Aerodynamic inertia matrix (no circulatory part)
    Ma =Ma(:,1);
    Ca = Ca(:,1);
    Ka = Ka(:,1);
    
    %Combining all these terms, we can compute the force and moment due to
    %the gust
    FM =zeros(2,length(theta(1,:)));
    for t=1:length(theta(1,:))
        FM(:,t) = Ma*theta(1,t) + Ca*theta(2,t) + Ka*theta(3,t) + W*[theta(4,t);theta(5,t)];
    end
end

function [x] = GustAngle( time, A, lambda1_0, lambda2_0, Evol_aero)
    %Generates a random angle of attack with a mean around the angle given
    %by theta_m. Following the implementation from fmamitrotta

    rng(11,"twister"); %random seed and generator for reproducibility

    %generating a random signal
    signal = randn(1,length(time));
    %generating angle theta from that signal and scaling it properly
    % if a different value for the average theta is wanted, it should be
    % shifted here
    theta = signal * A /sqrt(mean(signal.^2));
    %the derivatives of theta must also be obtained
    thetadot = gradient(theta);
    thetadotdot = gradient(thetadot);

    %Getting rid of the plunge depending terms in evolaero
    Evol_aero = [Evol_aero(1,1), Evol_aero(1,3), Evol_aero(1,5), Evol_aero(1,6);
                 Evol_aero(2,1), Evol_aero(2,3), Evol_aero(2,5), Evol_aero(2,6)];

    %Computing the lag states
    [t, lags] = ode45(@(t,lags)lagsevol(t,lags, thetadot, theta,time, Evol_aero),time,[lambda1_0; lambda2_0]);
    
    lags = lags';

    %Creating the output vector
    x = [thetadotdot; thetadot; theta; lags];

end

function [xdot] = GustyTimeEvolution(t,SM, x0, FM, inv_Mae, time)
    %Simulates the evolution of the system in the presence of a prescribed
    %gust

    %FM is the force due to a gust previously computed
    %Interpolating the force
    FM1_inter = interp1(time, FM(1,:), t);
    FM2_inter = interp1(time, FM(2,:), t);
    FM_inter = [FM1_inter;FM2_inter];

    gust = [inv_Mae*FM_inter; zeros(4,1)];

    xdot = SM*x0 + gust;
end

function dlagsdt = lagsevol(t,lags, thetadot, theta,time, Evol_aero)

    %Evolution of the lag states involved in the computation of the gust
    thetadot_inter = interp1(time,thetadot,t);
    theta_inter = interp1(time,theta,t);
    dlagsdt = Evol_aero*[thetadot_inter;theta_inter; lags(1); lags(2)];
end

