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
%Mw --> Mass per unit area of the wing
%Ma --> Mass per unit area of the aileron

%Data from Wright and cooper to validate against it (value of a has been
%adapted to match that data, since it is given in terms of xf)
Struct = struct('b', 1, 's', 7.5, 'a', -0.2, 'EI', 4*10^7, 'GJ', 8*10^6, 'x_alpha', 0.04, 'Mw', 400, 'Ma', 400, 'k_beta', 5*10^3);
Struct.tol = max(Struct.EI,Struct.GJ)/1e8;
%Parametrizing the structure with N points to integrate along the span
N = 500;
dy = Struct.s/N;
dx = 2*Struct.b/N;

% %Location of the aileron along the span of the wing
y0_aileron = 0.1*Struct.s; % one third of the span for the first position of the aileron
yf_aileron = 0.9 *Struct.s; % half of the span for the last position of the aileron

%Location of the aileron along the chord of the wing
x0_aileron = 1.6*Struct.b; %Starting at 0.8 chord or 1.6 b
b0 = (x0_aileron-Struct.b)/Struct.b; %Non-dimensional distance between the midchord to the beginning of the aileron, necessary to compute the T paramters

%Number of points -1 along the span to integrate
y_aileron = y0_aileron:dy:yf_aileron;
y = [0:dy:y0_aileron-dy, y_aileron, yf_aileron+dy:dy:Struct.s];
i0_aileron_y = find(y ==  y0_aileron);
if_aileron_y = find(y == yf_aileron);

%Along the chord
x_aileron = x0_aileron:dx:2*Struct.b;
x = [0:dx:x0_aileron-dx, x_aileron]; % Part of wing without aileron and then the one with it
i0_aileronx = find(x == x0_aileron);

%Parametrizing the chord of the structure
x = 0:dx:2*Struct.b;

%Creating a grid for double numerical integration
[X,Y] = meshgrid(x,y);

% %Bending mode shape (grid for double integration, vector for single
% integration)
 bending_deflection = (y/Struct.s).^2;
 bending = (Y/Struct.s).^2;

% %Torsional mode shape
torsion = (Y/Struct.s).*(X-(Struct.b+Struct.b*Struct.a)); %We need the distance wrt elastic axis for the structure
torsion_deflection = (y/Struct.s); %For the aerodynamic loads we only need the deflection

%Aileron mode shape (equivalent to torsion but displaced)
aileron_deflection = zeros(1,length(y)); %Assuming the deflection is constant along y and 1 only where it is present
for i =i0_aileron_y:if_aileron_y
    aileron_deflection(i) = 1;
end
%Since we have a ceiling function defining the mode shape, we need to
%include an if statement to evaluate it
aileron = zeros(length(y),length(x));
for j = i0_aileron_y: if_aileron_y
    for i = i0_aileronx : length(x)
        aileron(j,i) = (x(i)-x0_aileron);
    end
end

R = R_parameters(bending_deflection,torsion_deflection,aileron_deflection,y);

%Aerodynamic data, density of the air
rho =  1.225; % in kg/m^3

%Computing structural matrices (assuming no structural damping)
[A_struct, inv_A_struct,  C_struct, omega_nat] = matrices_struct(Struct,bending,bending_deflection,torsion,torsion_deflection,aileron, aileron_deflection,y,x, i0_aileronx, i0_aileron_y, if_aileron_y); %frequency given in rad/s

%Establishing the velocity range
V=linspace(1,300,183);
%Preallocating the output variables
omega= zeros(length(V),length(A_struct));
damp = zeros(length(V),length(A_struct));

%Iterating along the velocities
for j = 1:length(V)
    if j ==61
        g=0;
    end
    [omegas,dampings] = pk_theodorsen(V(j), Struct.b, omega_nat,A_struct, inv_A_struct, C_struct, Struct.a, rho,R,b0,Struct.tol);
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
plot(V,omega(:,2)/(2*pi),'-ko','MarkerIndices',1:10:length(omega(:,2)));
hold on
plot(V,omega(:,3)/(2*pi),'-k^','MarkerIndices',1:10:length(omega(:,3)));
hold on
% plot(V,omega(:,3),'b');

% xlabel('$V/(\omega_\alpha \, b)$', 'Interpreter','latex')
xlabel('$V$ [m/s]', 'Interpreter','latex')
legend('Torsion','Aileron', 'Bending', 'Interpreter','latex');
xlim([0 V(end)])
% ylim([1 8])
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
plot(V, damp(:,2)*100, '-ko','MarkerIndices',1:10:length(damp(:,2)));
hold on
plot(V, damp(:,3)*100, '-k^','MarkerIndices',1:10:length(damp(:,3)));
hold on

% plot(V, damp(:,3), 'b');

ylim([-20, 10])
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

function [omegas, dampings] = pk_theodorsen(V,b,omega_prev,A_struct, inv_A_struct, C_struct,a, rho, R,b0,tol)
    omegas = zeros(1,length(A_struct));
    dampings = zeros(1,length(A_struct));

    k_guess = omega_prev*b/V;
    for l = 1:length(A_struct)
        switch l
            case 1
                inc = 1;
                tol = tol;
            case 2
                inc = 1.6;
                tol = tol;
            case 3
                inc = 1;
                tol = tol;
        end
        k0 = k_guess(l);
        %Iterative procedure on k
        
        min_diff = 200;
        num = 0;
        r = 1;
        k_lambdas(r) = tol*k0 +k0+10;

        %Evaluating the resulting errors from the limits of the interval
            %Maximum
            lambdas_pos_kmax = sol_flutter(b, k0, a, rho, V, A_struct, inv_A_struct, C_struct,R,b0);
            k_lambdas_kmax = imag(lambdas_pos_kmax(1:end))*b/V;     
            [min_diff_kmax, r_kmax] = min(abs(k_lambdas_kmax(:)-k_guess(l)),[], 'all');
            %Minimum
            lambdas_pos_kmin = sol_flutter(b, k0, a, rho, V, A_struct, inv_A_struct, C_struct,R,b0);
            k_lambdas_kmin = imag(lambdas_pos_kmin(1:end))*b/V;     
            [min_diff_kmin, r_kmin] = min(abs(k_lambdas_kmin(:)-k_guess(l)),[], 'all');
            
            min_diff_prev = min_diff;
        while num<200 && (abs(k_lambdas(r)-k_guess(l))> tol )
            if (abs(min_diff) <= abs(min_diff_prev) && k0 >0)  || k_lambdas(r) ==0
                k0 = k0 +0.618* inc;
            else
                k0 = k0 -0.6*inc;
            end
            % k0 = k0 +inc;
            lambdas_pos = sol_flutter(b, k0, a, rho, V, A_struct, inv_A_struct, C_struct,R,b0);
            % if abs(imag(lambdas_pos(1) - lambdas_pos(2)))<0.01 || abs(imag(lambdas_pos(2) - lambdas_pos(3)))<0.01 || abs(imag(lambdas_pos(1) - lambdas_pos(3)))<0.01
            %     break;
            % end
            k_lambdas = imag(lambdas_pos(1:end))*b/V;
                
            [~,r] = min(abs(k_lambdas(:)-k_guess(l)),[], 'all');
            min_diff_prev= min_diff;
            [min_diff] = (k_lambdas(r)-k_guess(l));
            
            %THis is different in Theodorsen and Demasi
            omegas(l) = imag(lambdas_pos(r));
            dampings(l) = real(lambdas_pos(r))/abs(omegas(l)); 

            %Tracking number of iterations
            num = num+1;
         end
         % disp(num)
         % if num == 200
         %     disp('not converged');
         %     disp( V);
         %     disp(l);
         % end
    end
end

function [A_struct, inv_A_struct, C_struct, omega_nat] = matrices_struct(Struct, bending,bending_deflection,torsion,torsion_deflection, aileron, aileron_deflection,y,x, i0_aileronx, i0_aileron_y, if_aileron_y)

%Returns the structural matrices of the system
    b = Struct.b;
    a = Struct.a;
    % x_beta = Struct.x_beta;
    x_alpha = Struct.x_alpha;
    s = Struct.s;
    EI = Struct.EI;
    GJ = Struct.GJ;
    Mw = Struct.Mw;
    Ma = Struct.Ma;
    k_beta = Struct.k_beta;

    %Computing the coefficients of the structure integrating the mode
    %shapes

    %From the kinetic energy with no aileron present (only in the region
    %with no aileron)
    A_struct = zeros(3);
    %WHEN USING MESHGRID Y GOES FIRST THEN GOES X COORDINATE
    %Contribution to the inertia matrix from bending diff in kinetic energy
    A_struct11w = Mw*trapz(y,trapz(x(1:i0_aileronx-1),torsion(:,(1:i0_aileronx-1)).^2,2)); %Double integration of wing where no aileron is present
    A_struct13w = Mw*trapz(y,trapz(x(1:i0_aileronx-1),bending(:,(1:i0_aileronx-1)).*torsion(:,(1:i0_aileronx-1)),2)); %Contribution to combination of bending and torsion
    A_struct31w = A_struct13w;
    A_struct33w = Mw* trapz(y,trapz(x(1:i0_aileronx-1), bending(:,(1:i0_aileronx-1)).^2,2) ); %Contribution to q_bdotdot

    %Performing the integrals of the trailing edge left
    %we need to split into 3 regions spanwise
    %1st region before the aileron
    A_struct11wtrl = Mw*trapz(y(1:i0_aileron_y),trapz(x(i0_aileronx:end), torsion(1:i0_aileron_y,i0_aileronx: end).^2,2));
    A_struct13wtrl = Mw*trapz(y(1:i0_aileron_y),trapz(x(i0_aileronx:end), bending(1:i0_aileron_y,i0_aileronx: end).*torsion(1:i0_aileron_y,i0_aileronx: end),2));
    A_struct31wtrl = A_struct13wtrl;
    A_struct33wtrl = Mw*trapz(y(1:i0_aileron_y),trapz(x(i0_aileronx:end), bending(1:i0_aileron_y,i0_aileronx: end).^2,2));

    %Aileron region
    A_struct11a = Ma * trapz(y(i0_aileron_y:if_aileron_y),trapz(x(i0_aileronx:end),torsion(i0_aileron_y:if_aileron_y, i0_aileronx:end).^2,2));
    A_struct12a = Ma * trapz(y(i0_aileron_y:if_aileron_y),trapz(x(i0_aileronx:end),torsion(i0_aileron_y:if_aileron_y, i0_aileronx:end).*aileron(i0_aileron_y:if_aileron_y, i0_aileronx:end),2));
    A_struct13a = Ma * trapz(y(i0_aileron_y:if_aileron_y),trapz(x(i0_aileronx:end),torsion(i0_aileron_y:if_aileron_y, i0_aileronx:end).*bending(i0_aileron_y:if_aileron_y, i0_aileronx:end),2));

    A_struct21a = A_struct12a;
    A_struct22a = Ma * trapz(y(i0_aileron_y:if_aileron_y),trapz(x(i0_aileronx:end),aileron(i0_aileron_y:if_aileron_y, i0_aileronx:end).^2,2));
    A_struct23a = Ma * trapz(y(i0_aileron_y:if_aileron_y),trapz(x(i0_aileronx:end),bending(i0_aileron_y:if_aileron_y, i0_aileronx:end).*aileron(i0_aileron_y:if_aileron_y, i0_aileronx:end),2));

    A_struct31a = A_struct13a;
    A_struct32a = A_struct23a;
    A_struct33a = Ma * trapz(y(i0_aileron_y:if_aileron_y),trapz(x(i0_aileronx:end),bending(i0_aileron_y:if_aileron_y, i0_aileronx:end).^2,2));
    
    %Region after the aileron, wing, trailing edge, right
    A_struct11wtrr = Mw*trapz(y(if_aileron_y:end),trapz(x(i0_aileronx:end), torsion(if_aileron_y:end,i0_aileronx: end).^2,2));
    A_struct13wtrr = Mw*trapz(y(if_aileron_y:end),trapz(x(i0_aileronx:end), bending(if_aileron_y:end,i0_aileronx: end).*torsion(if_aileron_y:end,i0_aileronx: end),2));
    A_struct31wtrr = A_struct13wtrr;
    A_struct33wtrr = Mw*trapz(y(if_aileron_y:end),trapz(x(i0_aileronx:end), bending(if_aileron_y:end,i0_aileronx: end).^2,2));
    
    %Structural inertia matrix
     % A_struct = [s*M/3*(8*b^3/3 - 4*b^2*(b+a*b) + (b+a*b)^2 *2*b), s*M/4*(4*b^2/2 - 2*b*(b+a*b)) ;
     %            s*M/4*(4*b^2/2 - 2*b*(b+a*b)),M*s*2*b/5 ];

     %Computing the coefficients of these matrices from the mode shapes
     A11 = A_struct11w + A_struct11wtrl + A_struct11a + A_struct11wtrr;
     A12 = A_struct12a;
     A13 = A_struct13w + A_struct13wtrl + A_struct13a + A_struct13wtrr;

     A21 = A12;
     A22 = A_struct22a;
     A23 = A_struct23a;

     A31 = A13;
     A32 = A23;
     A33 = A_struct33w + A_struct33wtrl + A_struct33a + A_struct33wtrr;

     %Arranging into matrix form
     A_struct = [A11, A12, A13;
                A21, A22, A23;
                A31, A32, A33];

    %Its inverse
    inv_A_struct = inv(A_struct);
     
    %No structural damping being considered

    %Computing the mode shapes from the potential energy for the stiffness
    %matrix
    torsion_strain = gradient(torsion_deflection,y); %Diff uses right hand differences, so using gradient instead
    bending_strain = gradient(gradient(bending_deflection,y),y);
    aileron_strain = gradient(aileron_deflection(i0_aileron_y:if_aileron_y),y(i0_aileron_y:if_aileron_y));
    C_33 = EI*trapz(y,(bending_strain).^2);
    C_22 = k_beta* trapz(y(i0_aileron_y:if_aileron_y), (aileron_deflection(i0_aileron_y:if_aileron_y)).^2);
    C_11 = GJ*trapz(y,(torsion_strain).^2);

    %Structural stiffness matrix
    % C_struct = [GJ/s, 0;
    %             0,  4*EI/s^3];
    C_struct = [C_11, 0, 0;
                0, C_22, 0;
                0, 0, C_33];

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

function [ Q_re_dim, Q_im_dim] = matrices_aero(b, k0, a, rho, V, R,b0)

% %Precomputing necessary quantities
    c = b0;
    d = sqrt(1 - c^2);
    p = -1/3 * d^3;
    Ts = T_parameters(c);
% 
    %All coefficients with Theodorsen's function
    C_Theodorsen =  besselh(1,2,k0)/(besselh(1,2,k0) + 1j*besselh(0,2,k0));
    % Arranging coefficients intTo matrices to solve following the
    % procedure from Demasi's book
    F =real(C_Theodorsen);
    G = imag(C_Theodorsen);
% 
%     %  %From 2017 NASA's paper
   
    %Calling the function computing the parameters to convert using strip
    %theory
    

    %The coefficients adapting all components from Theodorsen's loads with
    %strip theory have been identified, and will be named now.
    % R_11 = s*b^2 *1/3;
    Q11_re = R(1,1) * (2 * pi) * ((-(1/8 + a^2) * k0^2)  - 2 * (a^2 - 1/4) * G * k0  - 2 * (a + 1/2) * F);
    Q11_im = R(1,1) * (2 * pi) * (((1/2 - a) * k0) + 2 * (a^2 - 1/4) * F * k0  - 2 * (a + 1/2) * G );

    Q21_re  = R(2,1) * (2 * pi) *(1/pi)* ((Ts(7) + (c - a) * Ts(1)) * k0^2 - Ts(12) * (1/2 - a) * G * k0 + Ts(12) * F);
    Q21_im  = R(2,1) * (2) * ((p - Ts(1) - Ts(4) / 2) * k0 + Ts(12) * (1/2 - a) * F * k0 + Ts(12) * G);
% 

    %Expressions from strip theory (coefficient for alpha, lift equation)
    %Multiply times b/4
      % R_31 = s*b/4;
      Q31_re = R(3,1) * (2 * pi) * ((a * k0^2 ) - 2 * (1/2 - a) * G * k0  + 2 * F ); 
      Q31_im = R(3,1) * (2 * pi) * (k0  + 2 * (1/2 - a) * F * k0  + 2 * G );

      Q12_re  = R(1,2) * (2*pi / ( pi)) * ((Ts(7) + (c - a) * Ts(1)) * k0^2 + (Ts(4) + Ts(10)) + (a + 1/2) * Ts(11) * G * k0 - 2 * (a + 1/2) * Ts(10) * F);
      Q12_im  = R(1,2) * (2*pi / ( pi)) * (-(2*p + (1/2 - a) * Ts(4)) * k0 -(a + 1/2) * Ts(11) * F * k0 - 2 * (a + 1/2) * Ts(10) * G);

      Q22_re  = R(2,2) *  (2 * pi / ( pi^2)) *(Ts(3) * k0^2 + (Ts(5) - Ts(4) * Ts(10)) - (Ts(11) * Ts(12) / 2) * G * k0 + Ts(10) * Ts(12) * F);
      Q22_im  = R(2,2) * (2 * pi / (pi^2)) * ((-Ts(4) * Ts(11) / 2) * k0 + (Ts(11) * Ts(12) / 2) * F * k0 + Ts(10) * Ts(12) * G);

      Q32_re  = R(3,2) * (2 * pi / ( pi)) * (Ts(1) * k0^2 - Ts(11) * G * k0 + 2 * Ts(10) * F);
      Q32_im  = R(3,2) * (2 * pi / ( pi)) * (-Ts(4) * k0 + Ts(11) * F * k0 + 2 * Ts(10) * G);

% 
    
    %Expressions from strip theory for Q13 (h and moment equation)
    %Multiply times  b *1/4
    % R_13 = s*b/4;
    Q13_re = R(1,3)* (2 * pi) * (a * k0^2 + 2 * (a + 1/2) * G * k0 );
    Q13_im = R(1,3) * (2 * pi) * (-2 * (a + 1/2) * F * k0 );

    Q23_re  = R(2,3) * (2 * pi / ( pi)) * (Ts(1) * k0^2 - Ts(12) * G * k0);
    Q23_im  = R(2,3) * (2 * pi / ( pi)) * (Ts(12) * F * k0);

%     Expressions from strip theory for Q33: h from lift equation
    % Multiply times 1/5
     % R_33 = s*1/5;
     Q33_re = R(3,3) * (2 * pi) * (-k0^2  - 2 * G * k0 );
     Q33_im = R(3,3) * (2 * pi) * (2 * F * k0 );

    % %Arranging the matrix
    Q_re(:,:) = [Q11_re , Q12_re , Q13_re ;
                Q21_re , Q22_re , Q23_re ;
                Q31_re , Q32_re , Q33_re ];

    Q_im(:,:) = [Q11_im , Q12_im , Q13_im ;
               Q21_im , Q22_im , Q23_im ;
               Q31_im , Q32_im , Q33_im ];

%Arranging the matrix and accounting for the presence of the span
% %Different order of the degrees of freedom now
    % Q_re(:,:) = [Q11_re , Q13_re ;
    %             Q31_re,  Q33_re ];
    % 
    % Q_im(:,:) = [Q11_im , Q13_im ;
    %             Q31_im , Q33_im];
%Matr
% % %Matrices are already dimensional, following derivation by NASA 2017
% 
    Q_re_dim = 0.5*rho*V^2*Q_re;
    Q_im_dim = 0.5*rho*V^2*Q_im;
end

function [lambdas_pos] = sol_flutter(b, k0, a, rho, V, A_struct, inv_A_struct, C_struct,R,b0)

    lambdas_pos = zeros(1,length(A_struct));
    [Q_re_dim, Q_im_dim] =  matrices_aero(b, k0, a, rho, V,R,b0);
    
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

function [R] = R_parameters(bending_deflection, torsion_deflection,aileron_deflection,y)
%This function obtains the R coefficients to multiply times the Q matrix
%and obtain the full flutter equations for the case of a finite wing with
%Including an aileron and applying strip theory

R = zeros(3);
%Integrals that we need to compute to find the coefficients in the
%aerodynamic loads matrix
int1 = torsion_deflection.^2;
R(1,1) = trapz(y, int1); %R11 in the matrix

%This integral should only be computed where the aileron is present
int2 = torsion_deflection.*aileron_deflection;
R(1,2) = trapz(y, int2);
R(2,1) = R(1,2);

int3 = bending_deflection.*torsion_deflection;
R(1,3) = trapz(y,int3); %Both R13 and R31
R(3,1) = trapz(y,int3); 

int4 = aileron_deflection.^2;
R(2,2) = trapz(y,int4);

int5 = aileron_deflection.*bending_deflection;
R(2,3) = trapz(y,int5);

int6 = bending_deflection.^2;
R(3,3) = trapz(y, int6); %R33

end