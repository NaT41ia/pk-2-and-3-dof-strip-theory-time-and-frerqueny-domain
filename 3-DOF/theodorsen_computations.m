% Theodorsen flutter speed analysis

%% Script to replicate the plots presented by Thedorosen in 1940

%{ 
This script implements the computations described by Theodorsen in https://ntrs.nasa.gov/api/citations/19930091762/downloads/19930091762.pdf
The model is described in terms of non-dimensional parameters, following
the description given in the paper above. 
 %}
clear; clc;
close all

% Data following definition from Theodorsen's paper 1940 
kappa = 0.25;
c = 0.6;
a = -0.4;
x_beta = 0.000000000000000000000000000001;
x_alpha = 0.2;
r_alpha = sqrt(0.25);
r_beta = sqrt(0.0012);
omega_h = sqrt(1/16);
omega_beta = omega_h*sqrt(3/2);
omega_alpha = 1;


% Compute d based on the given formula
xeh = c; % not chord but distance between elastic axis and location of the hinge
d = sqrt(1 - xeh^2);
p = -1/3 * d^3;

% Define each T parameter
T1 = -d * (2 + xeh^2) / 3 + xeh * acos(xeh);
T2 = xeh * d^2 - d * (1 + xeh^2) * acos(xeh) + xeh * (acos(xeh))^2;
T3 = -(1/8 + xeh^2) * (acos(xeh))^2 + xeh * d * acos(xeh) * (7 + 2 * xeh^2) / 4 - d^2 * (5 * xeh^2 + 4) / 8;
T4 = -acos(xeh) + xeh * d;
T5 = -d^2 - (acos(xeh))^2 + 2 * xeh * d * acos(xeh);
T6 = T2;
T7 = -(1/8 + xeh^2) * acos(xeh) + xeh * d * (7 + 2 * xeh^2) / 8;
T8 = -d * (1 + 2 * xeh^2) / 3 + xeh * acos(xeh);
T10 = d + acos(xeh);
T11 = acos(xeh) * (1 - 2 * xeh) + d * (2 - xeh);
T12 = d * (2 + xeh) - acos(xeh) * (1 + 2 * xeh);

% Defining the parameters as in Theodorsen's paper
% Coefficients for A
A_alpha1 = (r_alpha^2 / kappa) + (1 / 8 + a^2);
A_alpha2 = 1 / 2 - a;
A_alpha3 = 0;
A_beta1 = (r_beta^2 / kappa) - (T7 / pi) + (c - a) * ((x_beta / kappa) - (T1 / pi));
A_beta2 = (1 / pi) * (-2 * p - ((1 / 2) - a) * T4);
A_beta3 = (1 / pi) * (T4 + T10);
A_h1 = (x_alpha / kappa) - a;
A_h2 = 0;
A_h3 = 0;

% Coefficients for B
B_alpha1 = A_beta1;
B_alpha2 = (1 / pi) * (p - T1 - (1 / 2) * T4);
B_alpha3 = 0;
B_beta1 = (r_beta^2 / kappa) - (1 / pi^2) * T3;
B_beta2 = -(1 / (2 * pi^2)) * T4 * T11;
B_beta3 = (1 / pi^2) * (T5 - T4 * T10);
B_h1 = (x_beta / kappa) - (1 / pi) * T1;
B_h2 = 0;
B_h3 = 0;

% Coefficients for C
C_alpha1 = A_h1;
C_alpha2 = 1;
C_alpha3 = 0;
C_beta1 = B_h1;
C_beta2 = -(1 / pi) * T4;
C_beta3 = 0;
C_h1 = (1 / kappa) + 1;
C_h2 = 0;
C_h3 = 0;

% Solving first try

%Defining omegas (2DOF case)
OMEGA_alpha = 1;
OMEGA_beta = (omega_beta*r_beta/(omega_alpha*r_alpha))^2;
OMEGA_h = (omega_h/omega_alpha)^2 * 1/r_alpha^2;

%Defining the system coefficients
k = linspace(1/5,120,100000);
i = 1;
X1 = zeros(length(k),3);
X2 = ones(length(k),2);

while i < length(k)+1 %&& abs(X1(i,2) - X2(i))> 1e-8

    C_Theodorsen =  besselh(1,2,k(i))/(besselh(1,2,k(i)) + 1j*besselh(0,2,k(i)));
    F = real(C_Theodorsen);
    G = imag(C_Theodorsen);

    Ra_alpha = -A_alpha1 + (0.25-a^2)*(2*G)/k(i) - (0.5+a) *2*F/k(i)^2;
    % Ra_beta = -A_beta1 + (1 / k(i)^2) * A_beta3 + (1 / (k(i)*pi)) * (a + 0.5) * (T11 * G -2*(1/k(i))*T10*F);

    % Rb_alpha = -B_alpha1 - (1 / k(i)) * T12 / pi * ((0.5 - a) * G - (1 / k(i)) * F);
    Rb_beta = -B_beta1 + (1 / k(i)^2) * B_beta3 - (T12 / (k(i) *2 * pi^2)) * (T11 * G - 2 * T10 / k(i) * F);
    
    Rc_h = -C_h1 - 2*G/k(i);

    Ia_alpha = 1/k(i) * ((A_alpha2 -(0.5 + a )*2*G/k(i) - (0.25 - a^2 )* 2*F ));
    Ia_beta = 1/k(i) * (-(0.5 + a )*(T10*2*G/(pi*k(i))+ T11*F/pi) + A_beta2);

    Ib_alpha = 1/k(i)*(T12/(2*pi))*((2*G/k(i) + (0.5-a)*2*F)+B_alpha2);
    Ib_beta = 1/k(i)*(B_beta2 + T12*T10*2*G/(2*pi^2*k(i)) + (T12*T11*2*F)/(4*pi^2));
    
    Ic_h = 2*F/k(i);

    %Coefficients from 2dof cases

    %Case 1
    A1 = A_alpha1 * C_h1 - A_h1 * C_alpha1;

    B1 = A_alpha1 + (0.5+a)*C_alpha1 + (1/2 - a) * (-(1/2 + a) * C_h1 - A_h1);
    
    C1 = (A_h1 - A_alpha2) + (0.5+a)*(C_h1-C_alpha2);

    D1 = -(A_alpha2*C_h1 -A_h1 * C_alpha2);

    M1R = A1 + B1 *2*G/k(i) + C1*2*F/k(i)^2;

    M1I = 1/k(i)*(D1 + C1*2*G/k(i) - B1*2*F);

    % Case 2
    A2 = C_h1 * B_beta1 - C_beta1 * B_h1;
    A2_bar = -(C_h1 * B_beta3 - C_beta3 * B_h1);

    B2 = (B_beta1 - C_beta1 * T12 / (2 * pi))+ T11/(2*pi)*(C_h1 * T12 / (2 * pi) - B_h1 );
    B2_bar =  -(B_beta3 - T12/(2*pi)*C_beta3);
    
    C2 = - (B_beta2 - T12 / (2 * pi)*C_beta2) - T10/pi *(T12 / (2 * pi)*C_h1 - B_h1);

    D2 = - (C_h1*B_beta2 - C_beta2*B_h1);

    M2R = A2 + A2_bar*1/k(i)^2 +(B2 + B2_bar*1/k(i)^2)*2*G/k(i) + C2*2*F/(k(i)^2);

    M2I = 1/k(i)*(D2 +C2*2*G/k(i) - (B2 + B2_bar*1/k(i)^2)*2*F);

    %Case 3
    A3 = A_alpha1 * B_beta1 - A_beta1 * B_alpha1;
    A3_bar = -(A_alpha1 * B_beta3 - A_beta3 * B_alpha1) - (A_alpha2 * B_beta2 - A_beta2 * B_alpha2);

    B3 = (1/2 - a) * (-(1/2 + a) * B_beta1 - A_beta1 * T12 / (2 * pi))+ T11/(2*pi)*(A_alpha1 * T12 / (2 * pi) - B_alpha1 * -(1/2 + a) );
        M1 = [-(1/2 + a), A_beta3; T12 / (2 * pi),  B_beta3];
        M2 = [ A_alpha2, -(1/2 + a); B_alpha2, T12 / (2 * pi)];
        M3 = [-(1/2 + a),A_beta2;  T12 / (2 * pi) , B_beta2];
    B3_bar =  -(1/2 - a) * det(M1) - T10 / pi *det(M2) - det(M3);
    
        det1 = (-(1/2 + a) * B_beta2) - (T12/(2*pi)) * A_beta2;
        det2 = (T11/(2*pi)) * ((T12/(2*pi)) * A_alpha2 + (1/2 + a) * B_alpha2);
        det3 = (T10/pi) * ((T12/(2*pi)) * A_alpha1 + (1/2 + a) * B_alpha1);
        det4 = (-(1/2 + a) * B_beta1) - (T12/(2*pi)) * A_beta1;
    
    C3 = -(1/2 - a) * (det1) - det2 - det3 - det4;
    C3_bar = -(1/2 + a) * B_beta3 - A_beta3* T12/(2*pi);

    D3 = - (A_alpha1*B_beta2 - A_beta2*B_alpha1) - (A_alpha2*B_beta1 - B_alpha2*A_beta1);
    D3_bar = A_alpha2* B_beta3 - A_beta3*B_alpha2;

    M3R = A3 + A3_bar/k(i)^2 +(B3 + B3_bar/k(i)^2)*2*G/k(i) + (C3 + C3_bar/k(i)^2)*2*F/(k(i)^2);

    M3I =  1/k(i)*(D3 + D3_bar*1/k(i)^2 +(C3 + C3_bar*1/k(i)^2)*2*G/k(i) - (B3 + B3_bar*1/k(i)^2)*2*F);

    %Parameters of the three degree of freedom case
    R = -det([A_alpha1, A_beta1, A_h1; ...
              B_alpha1, B_beta1, B_h1; ...
              C_alpha1, C_beta1, C_h1]);
    R_bar = det([A_alpha1, A_beta3, A_h1; ...
                 B_alpha1, B_beta3, B_h1; ...
                 C_alpha1, C_beta3, C_h1]) + det([A_alpha2, A_beta2, A_h1; ...
                                                  B_alpha2, B_beta2, B_h1; ...
                                                  C_alpha2, C_beta2, C_h1]);

    S = -det([A_alpha1, A_beta1, -(0.5+a); ...
              B_alpha1, B_beta1, T12/(2*pi); ...
              C_alpha1, C_beta1, 1])  - (0.5-a)*det([-(0.5+a), A_beta1, A_h1; ...
                                                     T12/(2*pi), B_beta1, B_h1; ...
                                                     1, C_beta1, C_h1])  - T11/(2*pi)*det([A_alpha1,-(0.5+a), A_h1; ...
                                                                                           B_alpha1, T12/(2*pi), B_h1; ...
                                                                                           C_alpha1, 1, C_h1]);
    S_bar = det([A_alpha1, A_beta3, -(0.5+a); ...
              B_alpha1, B_beta3, T12/(2*pi); ...
              C_alpha1, C_beta3, 1])  + det([-(0.5+a), A_beta2, A_h1; ...
                                             T12/(2*pi), B_beta2, B_h1; ...
                                             1, C_beta2, C_h1]) +det([A_alpha2, A_beta2, -(0.5+a); ...
                                                                      B_alpha2, B_beta2, T12/(2*pi); ...
                                                                      C_alpha2, C_beta2, 1]) + T10/(pi)*det([A_alpha2,-(0.5+a), A_h1; ...
                                                                                                             B_alpha2, T12/(2*pi), B_h1; ...
                                                                                                             C_alpha2, 1, C_h1]) + (0.5-a)*  det([-(0.5+a), A_beta3, A_h1; ...
                                                                                                                                             T12/(2*pi), B_beta3, B_h1; ...
                                                                                                                                                  1, C_beta3, C_h1]);
    T = det([A_alpha1, A_beta2, -(0.5+a); ...
             B_alpha1, B_beta2, T12/(2*pi); ...
             C_alpha1, C_beta2, 1])  + det([-(0.5+a), A_beta1, A_h1; ...
                                             T12/(2*pi), B_beta1, B_h1; ...
                                             1, C_beta1, C_h1]) +det([A_alpha2, A_beta1, -(0.5+a); ...
                                                                      B_alpha2, B_beta1, T12/(2*pi); ...
                                                                      C_alpha2, C_beta1, 1]) + T10/(pi)*det([A_alpha1,-(0.5+a), A_h1; ...
                                                                                                             B_alpha1, T12/(2*pi), B_h1; ...
                                                                                                             C_alpha1, 1, C_h1]) + (0.5-a)*  det([-(0.5+a), A_beta2, A_h1; ...
                                                                                                                                                  T12/(2*pi), B_beta2, B_h1; ...
                                                                                                                                                  1, C_beta2, C_h1])+ T11/(2*pi)*det([A_alpha2,-(0.5+a), A_h1; ...
                                                                                                                                                                                    B_alpha2, T12/(2*pi), B_h1; ...
                                                                                                                                                                                    C_alpha2, 1, C_h1]);
    T_bar =  -det([-(0.5+a), A_beta3, A_h1; ...
                   T12/(2*pi), B_beta3, B_h1; ...
                    1, C_beta3, C_h1])- det([A_alpha2, A_beta3,-(0.5+a); ...
                                              B_alpha2,B_beta3, T12/(2*pi); ...
                                              C_alpha2, C_beta3,  1]);
    
    U =  det([A_alpha1, A_beta2, A_h1; ...
              B_alpha1, B_beta2, B_h1; ...
              C_alpha1, C_beta2, C_h1]) +det([A_alpha2, A_beta1, A_h1; ...
                                              B_alpha2, B_beta1, B_h1; ...
                                              C_alpha2, C_beta1, C_h1]);

    U_bar =-det([A_alpha2, A_beta3, A_h1; ...
                 B_alpha2, B_beta3, B_h1; ...
                 C_alpha2, C_beta3, C_h1]);

    %Real equation ( no structural damping)
    coeff1_X3 = OMEGA_alpha*OMEGA_beta*OMEGA_h;
    coeff1_X2 = OMEGA_alpha*OMEGA_beta*Rc_h + OMEGA_beta*OMEGA_h*Ra_alpha + OMEGA_h*OMEGA_alpha*Rb_beta;
    coeff1_X = OMEGA_alpha* M2R + OMEGA_beta*M1R + OMEGA_h * M3R;
    DR= R + R_bar/k(i)^2 + (S+S_bar/k(i)^2)*2*G/k(i) + (T+T_bar/k(i)^2)*2*F/k(i)^2;

    %Imaginary equation
    coeff2_X2 = OMEGA_alpha*OMEGA_beta*Ic_h + OMEGA_beta*OMEGA_h*Ia_alpha + OMEGA_h*OMEGA_alpha*Ib_beta;
    coeff2_X = OMEGA_alpha* M2I + OMEGA_beta*M1I + OMEGA_h * M3I;
    DI = 1/k(i)*(U + U_bar/k(i)^2 + (T+T_bar/k(i)^2)*2*G/k(i) - (S+S_bar/k(i)^2)*2*F);

    %Solving the equations
    X1(i,:) = roots([coeff1_X3, coeff1_X2, coeff1_X, DR]);
    X2(i,:) = roots([coeff2_X2, coeff2_X, DI]);

        i = i +1;
end

% hold on
% scatter(1./k,real(sqrt(X1(:,:))), "k.");
% scatter(1./k, sqrt(X2),'r.')
% grid on
% grid minor
% 
% 
% legend('Real part', 'Imaginary part')
% xlabel('$1/k$', Interpreter='latex');
% ylabel('$\sqrt{X}$', Interpreter='latex');
%% 
figure(1)
hold on
xlim([0 3]);
ylim = ([0 2]);
h1 = scatter(1./k, real(sqrt(X1(:,:))), 'k.', 'DisplayName', 'Real part');
h2 = scatter(1./k(250:100000), sqrt(X2(250:100000,:)), 'r.', 'DisplayName', 'Imaginary part');

% Force axis limits
axis([0 3 0 4.75]);

pbaspect([3 9.5 1]); 

grid on
grid minor

% Set legend correctly (assign legend only to the first handles)
legend([h1(1), h2(1)], 'Real part', 'Imaginary part', 'Location', 'best');
xlabel('$1/k$', 'Interpreter', 'latex');
ylabel('$\sqrt{X}$', 'Interpreter', 'latex');

hold off
%Compute flutter speed
%Assuming b constant at 1

%Numbers to modify to obtain the flutter speed
%Multiply times 1/k
%Divide by sqrtX

 fprintf("Vf =  r_alpha*omega_alpha/sqrt(kappa)*(1/k)/sqrt(X)");


 
