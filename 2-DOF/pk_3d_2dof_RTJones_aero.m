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
Struct = struct('b', 1, 's', 7.5, 'a', -0.2, 'EI', 4*10^7, 'GJ', 8*10^6, 'x_alpha', 0.04, 'M', 400);

%Value of tolerances to modify convergence
Struct.tol = (max(Struct.EI, Struct.GJ)/1e8);
%Parametrizing the structure with N points to integrate along the span
N = 500;
dy = Struct.s/N;
dx = Struct.b/N;
y = [0:dy:Struct.s];

%Parametrizing the chord of the structure
x = [0:dx:2*Struct.b];

%Creating a grid for double numerical integration
[X,Y] = meshgrid(x,y);

% %Bending mode shape (grid for double integration, vector for single
% integration)
 bending_deflection = (y/Struct.s).^2;
 bending = (Y/Struct.s).^2;

% %Torsional mode shape
torsion = (Y/Struct.s).*(X-(Struct.b+Struct.b*Struct.a)); %We need the distance wrt elastic axis for the structure
torsion_deflection = (y/Struct.s); %For the aerodynamic loads we only need the deflection
R = R_parameters(bending_deflection,torsion_deflection,y);

%Aerodynamic data, density of the air
rho =  1.225; % in kg/m^3

%Computing structural matrices (assuming no structural damping)
[A_struct, inv_A_struct,  C_struct, omega_nat] = matrices_struct(Struct,bending,bending_deflection,torsion,torsion_deflection,y,x); %frequency given in rad/s


%Establishing the velocity range
V=linspace(1,300,250);
%Preallocating the output variables
omega= zeros(length(V),length(A_struct));
damp = zeros(length(V),length(A_struct));


%Iterating along the velocities
for j = 1:length(V)
    [omegas,dampings] = pk_theodorsen(V(j), Struct.b, omega_nat,A_struct, inv_A_struct, C_struct, Struct.a, Struct.s, Struct.M, rho,R, Struct.tol);
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
 title("System in frequency domain")
% plot(V,omega(:,3),'b');

% xlabel('$V/(\omega_\alpha \, b)$', 'Interpreter','latex')
xlabel('$V$ [m/s]', 'Interpreter','latex')
legend('Mode 2', 'Mode 1', 'Interpreter','latex');
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

hold on
subplot(2,1,2)
plot(V, damp(:,1)*100, '-ksquare','MarkerIndices',1:10:length(damp(:,1)));
hold on
plot(V, damp(:,2)*100, '-k^','MarkerIndices',1:10:length(damp(:,2)));
hold on

V_freq = V;
omega_freq  = omega/(2*pi);
damp_freq =   damp*100;

save("FreqResults.mat","V_freq", "omega_freq", "damp_freq" )
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
omegas = omegas/(2*pi);

function [omegas, dampings] = pk_theodorsen(V,b,omega_prev,A_struct, inv_A_struct, C_struct,a, s, M, rho, R, tol)
    omegas = zeros(1,length(A_struct));
    dampings = zeros(1,length(A_struct));
    for l = 1:length(A_struct)
        switch l
            case 1
                k_guess = omega_prev(1)*b/V;     
                k_max = k_guess*50;
                k_min = k_guess/5;
                inc = 0.001;
                tol = tol;
                k0 = k_guess;
            case 2
                k_guess = omega_prev(2)*b/V;
                k_max = k_guess*2;
                k_min = k_guess/2;
                inc = 0.43244234234;
                tol = tol;
                k0 = k_guess;
            case 3
                k_guess =  omega_prev(3)*b/V;   
                k_max = k_guess*5;
                k_min = k_guess/2;
                inc = 0.5;
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
            lambdas_pos_kmax = sol_flutter(b,M, k0, a, rho, V, A_struct, inv_A_struct, C_struct,s,R);
            k_lambdas_kmax = imag(lambdas_pos_kmax(1:end))*b/V;     
            [min_diff_kmax, r_kmax] = min(abs(k_lambdas_kmax(:)-k_guess),[], 'all');
            %Minimum
            lambdas_pos_kmin = sol_flutter(b,M, k0, a, rho, V, A_struct, inv_A_struct, C_struct,s,R);
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
            lambdas_pos = sol_flutter(b,M, k0, a, rho, V, A_struct, inv_A_struct, C_struct,s,R);
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
         % disp(num)
         if num == 2000
             disp('not converged');
             disp( V);
             disp(l);
         end
    end

end

function [A_struct, inv_A_struct, C_struct, omega_nat] = matrices_struct(Struct, bending,bending_deflection,torsion,torsion_deflection,y,x)

%Returns the structural matrices of the system
    b = Struct.b;
    a = Struct.a;
    % x_beta = Struct.x_beta;
    x_alpha = Struct.x_alpha;
    s = Struct.s;
    EI = Struct.EI;
    GJ = Struct.GJ;
    M = Struct.M;

    %Computing the coefficients of the structure integrating the mode
    %shapes

    %From the kinetic energy
    A_struct = zeros(2);
    %Contribution to the inertia matrix from bending diff in kinetic energy
    A_struct11 = M*trapz(y,trapz(x,torsion.^2,2)); %Double integration
    A_struct13 = M*trapz(y,trapz(x,bending.*torsion,2)); %Contribution to combination of bending and torsion
    A_struct31 = A_struct13;
    A_struct33 = M* trapz(y,trapz(x, bending.^2,2) ); %Contribution to q_bdotdot

    A_struct = [A_struct11, A_struct13;
                A_struct31, A_struct33];
    
    %Structural inertia matrix
     % A_struct = [s*M/3*(8*b^3/3 - 4*b^2*(b+a*b) + (b+a*b)^2 *2*b), s*M/4*(4*b^2/2 - 2*b*(b+a*b)) ;
     %            s*M/4*(4*b^2/2 - 2*b*(b+a*b)),M*s*2*b/5 ];

     %Computing the coefficients of these matrices from the mode shapes

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
    % C_struct = [GJ/s, 0;
    %             0,  4*EI/s^3];
    C_struct = [C_11, 0;
                0, C_33];

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

function [ Q_re_dim, Q_im_dim] = matrices_aero(b, M, k0, a, rho, V,s, R)

%Computing theodorsen's function from Jones' approximation
    C_approx_Theodorsen = (0.5 * (1j*k0)^2 +0.2808*(1j*k0) +0.01365)/((1j*k0)^2 + 0.3455*(1j*k0) +0.01365);

    % Arranging coefficients intTo matrices to solve following the
    % procedure from Demasi's book
    F =real(C_approx_Theodorsen);
    G = imag(C_approx_Theodorsen);    

    %The coefficients adapting all components from Theodorsen's loads with
    %strip theory have been identified, and will be named now.
    % R_11 = s*b^2 *1/3;
    Q11_re = R(1,1) * (2 * pi) * ((-(1/8 + a^2) * k0^2)  - 2 * (a^2 - 1/4) * G * k0  - 2 * (a + 1/2) * F);
    Q11_im = R(1,1) * (2 * pi) * (((1/2 - a) * k0) + 2 * (a^2 - 1/4) * F * k0  - 2 * (a + 1/2) * G );

    %Expressions from strip theory (coefficient for alpha, lift equation)
    %Multiply times b/4
      % R_31 = s*b/4;
      Q31_re = R(3,1) * (2 * pi) * ((a * k0^2 ) - 2 * (1/2 - a) * G * k0  + 2 * F ); 
      Q31_im = R(3,1) * (2 * pi) * (k0  + 2 * (1/2 - a) * F * k0  + 2 * G );
    
    %Expressions from strip theory for Q13 (h and moment equation)
    %Multiply times  b *1/4
    % R_13 = s*b/4;
    Q13_re = R(1,3)* (2 * pi) * (a * k0^2 + 2 * (a + 1/2) * G * k0 );
    Q13_im = R(1,3) * (2 * pi) * (-2 * (a + 1/2) * F * k0 );


%     Expressions from strip theory for Q33: h from lift equation
    % Multiply times 1/5
     % R_33 = s*1/5;
     Q33_re = R(3,3) * (2 * pi) * (-k0^2  - 2 * G * k0 );
     Q33_im = R(3,3) * (2 * pi) * (2 * F * k0 );

%Arranging the matrix and accounting for the presence of the span
% %Different order of the degrees of freedom now
    Q_re(:,:) = [Q11_re , Q13_re ;
                Q31_re,  Q33_re ];

    Q_im(:,:) = [Q11_im , Q13_im ;
                Q31_im , Q33_im];
%Matr
% % %Matrices are already dimensional, following derivation by NASA 2017
% 
    Q_re_dim = 0.5*rho*V^2*Q_re;
    Q_im_dim = 0.5*rho*V^2*Q_im;
end

function [lambdas_pos] = sol_flutter(b,M, k0, a, rho, V, A_struct, inv_A_struct, C_struct,s,R)

    lambdas_pos = zeros(1,length(A_struct));
    [Q_re_dim, Q_im_dim] =  matrices_aero(b,M, k0, a, rho, V,s,R);
    
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

function [R] = R_parameters(bending, torsion,y)
%This function obtains the R coefficients to multiply times the Q matrix
%and obtain the full flutter equations for the case of a finite wing with
%NO aileron applying strip theory

%Initial definition of analytical mode shapes to compare with the results
%later on
%s = 7.5 is the span
% s = 7.5;
% %Location of the aileron along the span of the wing
% y0_aileron = 0.3*s; % one third of the span for the first position of the aileron
% yf_aileron = 0.5*s; % half of the span for the last position of the aileron

%Number of points -1 along the span to integrate

% y_aileron = y0_aileron:d:yf_aileron;
% y = [0:d:y0_aileron-d, y_aileron, yf_aileron+d:d:s];
% 
% i0_aileron = find(y ==  y0_aileron);
% if_aileron = find( y == yf_aileron);



% %Aileron's mode shape
% aileron = ones(1,length(y_aileron)); %Just a constant in W and C

R = zeros(3);
%Integrals that we need to compute to find the coefficients in the
%aerodynamic loads matrix
int1 = torsion.^2;
R(1,1) = trapz(y, int1); %R11 in the matrix

% %This integral should only be computed where the aileron is present
% int2 = torsion(1,i0_aileron:if_aileron).*aileron;
% R(1,2) = trapz(y_aileron, int2);
% R(2,1) = R(1,2);

int3 = bending.*torsion;
R(1,3) = trapz(y,int3); %Both R13 and R31
R(3,1) = trapz(y,int3); 

% int4 = aileron.^2;
% R(2,2) = trapz(y_aileron,int4);
% 
% int5 = aileron.*bending(i0_aileron:if_aileron);
% R(2,3) = trapz(y_aileron,int5);

int6 = bending.^2;
R(3,3) = trapz(y, int6); %R33

end