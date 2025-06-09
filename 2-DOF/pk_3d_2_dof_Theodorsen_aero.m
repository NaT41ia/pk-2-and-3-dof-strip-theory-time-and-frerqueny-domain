%% Computing forces from Theodorsens equations
close all
clc;
clear;
format long

% %Establishing data from the structure in a structure (dimensional data)
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
Struct = struct('b', 1, 's', 7.5, 'a', -0.04, 'EI', 2*10^7, 'GJ', 2*10^6, 'x_alpha', 0.04, 'M', 200);

%Tolerances for better pk convergence
Struct.tol = max(Struct.EI,Struct.GJ)/1e8;
%Aerodynamic data, density of the air
rho =  1.225; % in kg/m^3

%Computing structural matrices (assuming no structural damping)
[A_struct, inv_A_struct,  C_struct, omega_nat] = matrices_struct(Struct); %frequency given in rad/s


%Establishing the velocity range
V=linspace(1,120,100);
%Preallocating the output variables
omega= zeros(length(V),length(A_struct));
damp = zeros(length(V),length(A_struct));

%Iterating along the velocities
for j = 1:length(V)
    [omegas,dampings] = pk_theodorsen(V(j), Struct.b, omega_nat,A_struct, inv_A_struct, C_struct, Struct.a, Struct.s, Struct.M, rho, Struct.tol);
    for n = 1:length(omegas)
        if omegas(n) == 0
            omegas(n) = omega(j-1,n);
            dampings(n) = damp(j-1,n);
        end
    end
    omega(j,:) = omegas;
    damp(j,:) = dampings;
end

figure(1)
hold on
subplot(2,1,1)
plot(V,omega(:,1)/(2*pi),'-ksquare','MarkerIndices',1:10:length(omega(:,1)));
hold on
 plot(V,omega(:,2)/(2*pi),'-k^','MarkerIndices',1:10:length(omega(:,2)));
 hold on
% plot(V,omega(:,3),'b');

% xlabel('$V/(\omega_\alpha \, b)$', 'Interpreter','latex')
xlabel('$V$ [m/s]', 'Interpreter','latex')
legend('Torsion', 'Bending', 'Interpreter','latex');
xlim([0 V(end)])
grid on
grid minor
ax = gca;
% ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '-';
ax.GridAlpha = 0.5;
% ax.Layer = 'top';
ylabel('$\omega$ [Hz]', 'Interpreter','latex')
pbaspect([3 1 1])

hold on
subplot(2,1,2)
plot(V, -damp(:,1)*100, '-ksquare','MarkerIndices',1:10:length(damp(:,1)));
hold on
plot(V, -damp(:,2)*100, '-k^','MarkerIndices',1:10:length(damp(:,2)));
hold on

% plot(V, damp(:,3), 'b');

ylim([-5, 15])
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

function [omegas, dampings] = pk_theodorsen(V,b,omega_prev,A_struct, inv_A_struct, C_struct,a, s, M, rho, tol)
    omegas = zeros(1,length(A_struct));
    dampings = zeros(1,length(A_struct));
    for l = 1:length(A_struct)
        switch l
            case 1
                k_guess = omega_prev(1)*b/V;     
                k_max = k_guess*50;
                k_min = k_guess/5;
                inc = 0.00001;
                tol = tol;
                k0 = k_guess;
            case 2
                k_guess = omega_prev(2)*b/V;
                k_max = k_guess*2;
                k_min = k_guess/2;
                inc = 0.00001;
                tol = tol;
                k0 = k_guess;
            case 3
                k_guess =  omega_prev(3)*b/V;   
                k_max = k_guess*5;
                k_min = k_guess/2;
                inc = 3/10^(kappa/0.025);
                tol = tol;
                k0 = k_guess;
        end
        %Iterative procedure on k
        
        
        min_diff = 200;
        num = 0;
        r = 1;
        k_lambdas(r) = tol*k_guess +k_guess+10;
        cont_a = 0;
        cont_b = 0;

        %Evaluating the resulting errors from the limits of the interval
            %Maximum
            lambdas_pos_kmax = sol_flutter(b,M, k0, a, rho, V, A_struct, inv_A_struct, C_struct,s);
            k_lambdas_kmax = imag(lambdas_pos_kmax(1:end))*b/V;     
            [min_diff_kmax, r_kmax] = min(abs(k_lambdas_kmax(:)-k_guess),[], 'all');
            %Minimum
            lambdas_pos_kmin = sol_flutter(b,M, k0, a, rho, V, A_struct, inv_A_struct, C_struct,s);
            k_lambdas_kmin = imag(lambdas_pos_kmin(1:end))*b/V;     
            [min_diff_kmin, r_kmin] = min(abs(k_lambdas_kmin(:)-k_guess),[], 'all');
            
            min_diff_prev = min_diff;
          % while num<10 &&  min_diff>1e-3
        while num<2000 && (abs(k_lambdas(r)-k_guess)> tol )
            if (abs(min_diff) < abs(min_diff_prev) && k0 >0)  || k_lambdas(r) ==0
                k0 = k0 - inc;
            else
                  k0 = k0 +0.5*inc;
            end
            % k0 = k0 +inc;
            lambdas_pos = sol_flutter(b,M, k0, a, rho, V, A_struct, inv_A_struct, C_struct,s);
            % if abs(imag(lambdas_pos(1) - lambdas_pos(2)))<0.01 || abs(imag(lambdas_pos(2) - lambdas_pos(3)))<0.01 || abs(imag(lambdas_pos(1) - lambdas_pos(3)))<0.01
            %     break;
            % end
            k_lambdas = imag(lambdas_pos(1:end))*b/V;
                
            [~,r] = min(abs(k_lambdas(:)-k_guess),[], 'all');
            min_diff_prev= min_diff;
            [min_diff] = (k_lambdas(r)-k_guess);
            
            %THis is different in Theodorsen and Demasi
            omegas(l) = imag(lambdas_pos(r));
            dampings(l) = real(lambdas_pos(r))/abs(omegas(l)); 

            %Tracking number of iterations
            num = num+1;
        end
         if num == 2000
             disp('not converged');
             disp( V);
             disp(l);
         end
    end

end

function [A_struct, inv_A_struct, C_struct, omega_nat] = matrices_struct(Struct)

%Returns the structural matrices of the system
    b = Struct.b;
    a = Struct.a;
    % x_beta = Struct.x_beta;
    x_alpha = Struct.x_alpha;
    s = Struct.s;
    EI = Struct.EI;
    GJ = Struct.GJ;
    M = Struct.M;
    %Structural inertia matrix
    % A_struct = [M*s*2*b/5, s*M/4*(4*b^2/2 - 2*b*(b+a*b));
    %             s*M/4*(4*b^2/2 - 2*b*(b+a*b)), s*M/3*(8*b^3/3 - 4*b^2*(b+a*b) + (b+a*b)^2 *2*b)];

     A_struct = [s*M/3*(8*b^3/3 - 4*b^2*(b+a*b) + (b+a*b)^2 *2*b), s*M/4*(4*b^2/2 - 2*b*(b+a*b)) ;
                s*M/4*(4*b^2/2 - 2*b*(b+a*b)),M*s*2*b/5 ];

    %Its inverse
    inv_A_struct = inv(A_struct);
     
    %No structural damping being considered
    
    %Structural stiffness matrix
    C_struct = [GJ/s, 0;
                0,  4*EI/s^3];

    %Finding the natural frequencies of the system
    lam = eig(inv_A_struct*C_struct);
    omega_nat = sqrt(lam);
end

function [Ts] = T_parameters(c)
    % Compute d based on the given formula
    
    d = sqrt(1 - c^2);
    p = -1/3 * d^3;
    
    % Define each T parameter
    Ts(1) = -d * (2 + c^2) / 3 + c * acos(c);
    Ts(2) = c * d^2 - d * (1 + c^2) * acos(c) + c * (acos(c))^2;
    Ts(3) = -(1/8 + c^2) * (acos(c))^2 + c * d * acos(c) * (7 + 2 * c^2) / 4 - d^2 * (5 * c^2 + 4) / 8;
    Ts(4) = -acos(c) + c * d;
    Ts(5) = -d^2 - (acos(c))^2 + 2 * c * d * acos(c);
    Ts(6) = Ts(2);
    Ts(7) = -(1/8 + c^2) * acos(c) + c * d * (7 + 2 * c^2) / 8;
    Ts(8) = -d * (1 + 2 * c^2) / 3 + c * acos(c);
    Ts(9) = 0; %It is never used, just included here to maintain the numbering
    Ts(10) = d + acos(c);
    Ts(11) = acos(c) * (1 - 2 * c) + d * (2 - c);
    Ts(12) = d * (2 + c) - acos(c) * (1 + 2 * c);
end

function [ Q_re_dim, Q_im_dim] = matrices_aero(b, M, k0, a, rho, V,s)

% %Precomputing necessary quantities
%     % d = sqrt(1 - c^2);
%     % p = -1/3 * d^3;
%     % 
%     % Ts = T_parameters(c);
% 
    %All coefficients with Theodorsen's function
    C_Theodorsen =  besselh(1,2,k0)/(besselh(1,2,k0) + 1j*besselh(0,2,k0));
    % Arranging coefficients intTo matrices to solve following the
    % procedure from Demasi's book
    F =real(C_Theodorsen);
    G = imag(C_Theodorsen);
% 
%     %  %From 2017 NASA's paper
   
    %Expressions from strip theory (alpha and moment )
    %Essentially supressing division by M and multiplying times -b^2 and 1/3
    Q11_re = (2 * pi) * ((-(1/8 + a^2) * k0^2)*b^2 *1/3  - 2 * (a^2 - 1/4) * G * k0 * b^2 * 1/3 - 2 * (a + 1/2) * F*b^2 *1/3);
    Q11_im = (2 * pi) * (((1/2 - a) * k0)*b *1/3 + 2 * (a^2 - 1/4) * F * k0 *b^2 *1/3 - 2 * (a + 1/2) * G )*b*1/3;
%     % Q21_re  = (2 * pi / M) *(1/pi)* ((T7 + (c - a) * T1) * k0^2 - T12 * (1/2 - a) * G * k0 + T12 * F);
%     % Q21_im  = (2 / M) * ((p - T1 - T4 / 2) * k0 + T12 * (1/2 - a) * F * k0 + T12 * G);
% 

    %Expressions from strip theory (coefficient for alpha, lift equation)
    %Multiply times -b/4
      Q31_re = (2 * pi) * ((a * k0^2 * b * 1/4) - 2 * (1/2 - a) * G * k0 * b * 1/4 + 2 * F *b *1/4); 
      Q31_im = (2 * pi) * (k0 *b * 1/4 + 2 * (1/2 - a) * F * k0 *b *1/4 + 2 * G *b *1/4);

% 
    
    %Expressions from strip theory for Q13 (h and moment equation)
    %Multiply times  -b *1/4
    Q13_re = (2 * pi) * (a * k0^2 *b * 1/4 + 2 * (a + 1/2) * G * k0 *b *1/4);
    Q13_im = (2 * pi) * (-2 * (a + 1/2) * F * k0 *b *1/4);

%     Expressions from strip theory for Q33: h from lift equation
    % Multiply times -1/5
     Q33_re = (2 * pi) * (-k0^2 * 1/5 - 2 * G * k0 *1/5);
     Q33_im = (2 * pi) * (2 * F * k0 * 1/5);

%Arranging the matrix and accounting for the presence of the span
%Different order of the degrees of freedom now
    Q_re(:,:) = [Q11_re , Q13_re ;
                Q31_re,  Q33_re ];

    Q_im(:,:) = [Q11_im , Q13_im ;
                Q31_im , Q33_im];
%Matr
% % %Matrices are already dimensional, following derivation by NASA 2017
% 
    Q_re_dim = s*0.5*rho*V^2*Q_re;
    Q_im_dim = s*0.5*rho*V^2*Q_im;

end

function [lambdas_pos] = sol_flutter(b,M, k0, a, rho, V, A_struct, inv_A_struct, C_struct,s)
    lambdas_pos = zeros(1,length(A_struct));
    [Q_re_dim, Q_im_dim] =  matrices_aero(b,M, k0, a, rho, V,s);
    
    G =  [zeros(length(A_struct)), eye(length(A_struct));
          -inv_A_struct*(C_struct+ Q_re_dim),-inv_A_struct*Q_im_dim*b/(k0*V)];
    lambdas = eig(G);
    %Getting rid of negative imaginary parts
    index = 1;
    if imag(lambdas(:)) == 0
        return;
    end
    lambdas_pos = zeros(1,3);
    for i = 1:length(lambdas)
        if imag(lambdas(i)) > 0
            lambdas_pos(index) = lambdas(i);
            index = index +1;
        end
    end
end