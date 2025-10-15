function [F, A] = construct_F2andA2(solution, variableLocation)
%--- This function returnes the value of the ODE vector field at the given solution matrices. The ordering is in x1 direction

global V_s V_u1 V_u2 U kappa K1 K2  Q_opt R1 R2
global AlleeFlag J1 J2

segmentSize = numel(solution(:,1)); % number of (A)ij segments
numNonZero = 6^2 * segmentSize; %number of non-zero elements in A.
odeSystemSize = numel(solution);

N1 = solution(:,1);
Q1 = solution(:,2);
V1 = solution(:,3);
N2 = solution(:,4);
Q2 = solution(:,5);
V2 = solution(:,6);

%% Definign the additional growth term to model Allee effect with a qubic growth=================
if AlleeFlag
    % dB1_dN1 = 4 * K1(:) ./ ( (K1(:) - J1(:)).^2 );
    % dB2_dN2 = 4 * K2(:) ./ ( (K2(:) - J2(:)).^2 );
    % B1 = dB1_dN1 .* J1(:) .* ( N1 ./ J1(:) - 1);
    % B2 = dB2_dN2 .* J2(:) .* ( N2 ./ J2(:) - 1);
    
    expN1mJ1 = exp(-(N1 - J1) / 0.05);
    expN2mJ2 = exp(-(N2 - J2) / 0.05);
    B1 = 2.6 ./ ( 1 + expN1mJ1 ) - (2.6/2);
    B2 = 2.6 ./ ( 1 + expN2mJ2 ) - (2.6/2);
    dB1_dN1 = (2.6/0.05) * expN1mJ1 ./ ( 1 + expN1mJ1 ).^2;
    dB2_dN2 = (2.6/0.05) * expN2mJ2 ./ ( 1 + expN2mJ2 ).^2;

else %---Set Bi = 1 and dBi_dNi = 0 to remove the Allee effect
    B1 = 1;
    B2 = 1;
    dB1_dN1 = 0;
    dB2_dN2 = 0;
end

%% Calculating nonlinear terms =============================================
Vbar11 = V_u1; 
Vbar12 = (V_u1 + V_u2)/2;
Vbar21 = (V_u2 + V_u1)/2;
Vbar22 = V_u2;
Lambda11 = 1;
Lambda12 = sqrt(V_u1/Vbar12);
Lambda21 = sqrt(V_u2/Vbar21);
Lambda22 = 1;

Q1m2kVbar11 = Q1 - 2 * kappa * Vbar11;
Q1m2kVbar21 = Q1 - 2 * kappa * Vbar21;
Q2m2kVbar12 = Q2 - 2 * kappa * Vbar12;
Q2m2kVbar22 = Q2 - 2 * kappa * Vbar22;

Q1mQ2p2kVbar12 = Q1 - Q2m2kVbar12;
Q2mQ1p2kVbar21 = Q2 - Q1m2kVbar21;
Q1mQ2p2kVbar12_squared = Q1mQ2p2kVbar12.^2;
Q2mQ1p2kVbar21_squared = Q2mQ1p2kVbar21.^2;

V1_squared = V1.^2;
V2_squared = V2.^2;

V1p2Vbar11 = V1 + 2 * Vbar11;
V1p2Vbar21 = V1 + 2 * Vbar21;
V2p2Vbar12 = V2 + 2 * Vbar12;
V2p2Vbar22 = V2 + 2 * Vbar22;
V1p2Vbar11_squared = V1p2Vbar11.^2;
V1p2Vbar21_squared = V1p2Vbar21.^2;
V2p2Vbar12_squared = V2p2Vbar12.^2;
V2p2Vbar22_squared = V2p2Vbar22.^2;

V1pV1p2Vbar11 = V1 + V1p2Vbar11;
V1pV2p2Vbar12 = V1 + V2p2Vbar12;
V2pV1p2Vbar21 = V2 + V1p2Vbar21;
V2pV2p2Vbar22 = V2 + V2p2Vbar22;
V1pV1p2Vbar11_squared = V1pV1p2Vbar11.^2;
V1pV2p2Vbar12_squared = V1pV2p2Vbar12.^2;
V2pV1p2Vbar21_squared = V2pV1p2Vbar21.^2;
V2pV2p2Vbar22_squared = V2pV2p2Vbar22.^2;

Q_opt_squared = Q_opt(:).^2;
Q1_squared = Q1.^2;
Q2_squared = Q2.^2;
Q1_cubed = Q1_squared .* Q1;
Q2_cubed = Q2_squared .* Q2;
Q_optQ1 = Q_opt(:) .* Q1;
Q_optQ2 = Q_opt(:) .* Q2;
twoQ_optV1 = 2 * Q_opt(:) .* V1;
twoQ_optV2 = 2 * Q_opt(:) .* V2;
Q_optQ1_squared = Q_opt(:) .* Q1_squared;
Q_optQ2_squared = Q_opt(:) .* Q2_squared;
Q_opt_squaredQ1 = Q_opt_squared .* Q1;
Q_opt_squaredQ2 = Q_opt_squared .* Q2;
Q1mQ_opt = Q1 - Q_opt(:);
Q2mQ_opt = Q2 - Q_opt(:);
V1mQ1_squared = V1 - Q1_squared;
V2mQ2_squared = V2 - Q2_squared;

C11 = Lambda11 * sqrt(2*Vbar11) * exp( kappa^2 * Vbar11 ) ./ sqrt(V1pV1p2Vbar11);
C12 = Lambda12 * sqrt(2*Vbar12) * exp( kappa^2 * Vbar12 ) ./ sqrt(V1pV2p2Vbar12);
C21 = Lambda21 * sqrt(2*Vbar21) * exp( kappa^2 * Vbar21 ) ./ sqrt(V2pV1p2Vbar21);
C22 = Lambda22 * sqrt(2*Vbar22) * exp( kappa^2 * Vbar22 ) ./ sqrt(V2pV2p2Vbar22);
M11 = exp( -(2 * kappa * Vbar11)^2 ./ (2 * V1pV1p2Vbar11) );
M12 = exp( -Q1mQ2p2kVbar12_squared ./ (2 * V1pV2p2Vbar12) );
M21 = exp( -Q2mQ1p2kVbar21_squared ./ (2 * V2pV1p2Vbar21) );
M22 = exp( -(2 * kappa * Vbar22)^2 ./ (2 * V2pV2p2Vbar22) );
L11 = ( V1 .* Q1m2kVbar11 + V1p2Vbar11 .* Q1 ) ./ V1pV1p2Vbar11;
L12 = ( V1 .* Q2m2kVbar12 + V2p2Vbar12 .* Q1 ) ./ V1pV2p2Vbar12;
L21 = ( V2 .* Q1m2kVbar21 + V1p2Vbar21 .* Q2 ) ./ V2pV1p2Vbar21;
L22 = ( V2 .* Q2m2kVbar22 + V2p2Vbar22 .* Q2 ) ./ V2pV2p2Vbar22;
S11 = V1 .* V1p2Vbar11 ./ V1pV1p2Vbar11 + L11 .* (L11 - 2*Q1);
S12 = V1 .* V2p2Vbar12 ./ V1pV2p2Vbar12 + L12 .* (L12 - 2*Q1);
S21 = V2 .* V1p2Vbar21 ./ V2pV1p2Vbar21 + L21 .* (L21 - 2*Q2);
S22 = V2 .* V2p2Vbar22 ./ V2pV2p2Vbar22 + L22 .* (L22 - 2*Q2);
E1 = ( twoQ_optV1 + 2 * Q_optQ1_squared - Q_opt_squaredQ1 - 3 * V1 .* Q1 - Q1_cubed ) / (2*V_s);
E2 = ( twoQ_optV2 + 2 * Q_optQ2_squared - Q_opt_squaredQ2 - 3 * V2 .* Q2 - Q2_cubed ) / (2*V_s);
Y1 = U + ( twoQ_optV1 .* Q1 - 2 * Q_opt(:) .* Q1_cubed - Q_opt_squared .* V1mQ1_squared - 3 * V1_squared + Q1_cubed .* Q1 ) / (2*V_s);
Y2 = U + ( twoQ_optV2 .* Q2 - 2 * Q_opt(:) .* Q2_cubed - Q_opt_squared .* V2mQ2_squared - 3 * V2_squared + Q2_cubed .* Q2 ) / (2*V_s);

B1R1 = B1(:) .* R1(:);
B2R2 = B2(:) .* R2(:);
B1R1_K1 =  B1R1 ./ K1(:);
B2R2_K2 =  B2R2 ./ K2(:);
B1R1N1_K1 = B1R1_K1 .* N1;
B1R1N2_K1 = B1R1_K1 .* N2;
B2R2N1_K2 = B2R2_K2 .* N1;
B2R2N2_K2 = B2R2_K2 .* N2;
B1R1C11N1_K1 = B1R1N1_K1 .* C11;
B1R1C12N2_K1 = B1R1N2_K1 .* C12;
B2R2C21N1_K2 = B2R2N1_K2 .* C21;
B2R2C22N2_K2 = B2R2N2_K2 .* C22;
B1R1M11N1_K1 = B1R1N1_K1 .* M11;
B1R1M12N2_K1 = B1R1N2_K1 .* M12;
B2R2M21N1_K2 = B2R2N1_K2 .* M21;
B2R2M22N2_K2 = B2R2N2_K2 .* M22;
B1R1M11C11N1_K1 = B1R1C11N1_K1 .* M11;
B1R1M12C12N2_K1 = B1R1C12N2_K1 .* M12;
B2R2M21C21N1_K2 = B2R2C21N1_K2 .* M21;
B2R2M22C22N2_K2 = B2R2C22N2_K2 .* M22;
B1R1M11C11_K1 = B1R1_K1 .* M11 .* C11;
B1R1M12C12_K1 = B1R1_K1 .* M12 .* C12;
B2R2M21C21_K2 = B2R2_K2 .* M21 .* C21;
B2R2M22C22_K2 = B2R2_K2 .* M22 .* C22;

if AlleeFlag
    R1_K1 =  R1(:) ./ K1(:);
    R2_K2 =  R2(:) ./ K2(:);
    R1N1_K1 = R1_K1 .* N1;
    R1N2_K1 = R1_K1 .* N2;
    R2N1_K2 = R2_K2 .* N1;
    R2N2_K2 = R2_K2 .* N2;
    R1M11C11N1_K1 = R1N1_K1 .* M11 .* C11;
    R1M12C12N2_K1 = R1N2_K1 .* M12 .* C12;
    R2M21C21N1_K2 = R2N1_K2 .* M21 .* C21;
    R2M22C22N2_K2 = R2N2_K2 .* M22 .* C22;
else
    R1M11C11N1_K1 = B1R1M11C11N1_K1;
    R1M12C12N2_K1 = B1R1M12C12N2_K1;
    R2M21C21N1_K2 = B2R2M21C21N1_K2;
    R2M22C22N2_K2 = B2R2M22C22N2_K2;
end

G1 = B1R1 - B1R1M11C11N1_K1 - B1R1M12C12N2_K1 - ( Q1mQ_opt.^2 + V1 ) / (2 * V_s);
G2 = B2R2 - B2R2M21C21N1_K2 - B2R2M22C22N2_K2 - ( Q2mQ_opt.^2 + V2 ) / (2 * V_s);
B1R1mG1 =  B1R1 - G1;
B2R2mG2 =  B2R2 - G2;
H1 = B1R1mG1 .* Q1 - B1R1M11C11N1_K1 .* L11 - B1R1M12C12N2_K1 .* L12 + E1;
H2 = B2R2mG2 .* Q2 - B2R2M21C21N1_K2 .* L21 - B2R2M22C22N2_K2 .* L22 + E2;
W1 = B1R1mG1 .* V1mQ1_squared - B1R1M11C11N1_K1 .* S11 - B1R1M12C12N2_K1 .* S12 + Y1;
W2 = B2R2mG2 .* V2mQ2_squared - B2R2M21C21N1_K2 .* S21 - B2R2M22C22N2_K2 .* S22 + Y2;

%% Calculating F3 ==========================================================
F = zeros(odeSystemSize, 1);
F(variableLocation(1,:)) = G1 .* N1;
F(variableLocation(2,:)) = H1;
F(variableLocation(3,:)) = W1;
F(variableLocation(4,:)) = G2 .* N2;
F(variableLocation(5,:)) = H2;
F(variableLocation(6,:)) = W2;

%% Calculating A3 ==========================================================
dC11_dV1 = -C11 ./ V1pV1p2Vbar11;
dC12_dV1 = -C12 ./ (2 * V1pV2p2Vbar12);
dC12_dV2 = dC12_dV1;
dC21_dV1 = -C21 ./ (2 * V2pV1p2Vbar21);
dC21_dV2 = dC21_dV1;
dC22_dV2 = -C22 ./ V2pV2p2Vbar22;

dM12_dQ1 = -Q1mQ2p2kVbar12 ./ V1pV2p2Vbar12 .* M12;
dM12_dQ2 = -dM12_dQ1;
dM21_dQ1 = Q2mQ1p2kVbar21 ./ V2pV1p2Vbar21 .* M21;
dM21_dQ2 = -dM21_dQ1;

dM11_dV1 = (2 * kappa * Vbar11)^2 ./ V1pV1p2Vbar11_squared .* M11;
dM12_dV1 = Q1mQ2p2kVbar12_squared ./ (2 * V1pV2p2Vbar12_squared) .* M12;
dM12_dV2 = dM12_dV1;
dM21_dV1 = Q2mQ1p2kVbar21_squared ./ (2 * V2pV1p2Vbar21_squared) .* M21;
dM21_dV2 = dM21_dV1;
dM22_dV2 = (2 * kappa * Vbar22)^2 ./ V2pV2p2Vbar22_squared .* M22;

dL11_dQ1 = 1;
dL12_dQ1 = V2p2Vbar12 ./ V1pV2p2Vbar12;
dL12_dQ2 = V1 ./ V1pV2p2Vbar12;
dL21_dQ1 = V2 ./ V2pV1p2Vbar21;
dL21_dQ2 = V1p2Vbar21 ./ V2pV1p2Vbar21;
dL22_dQ2 = 1;

dL11_dV1 = -(4 * kappa * Vbar11^2) ./ V1pV1p2Vbar11_squared;
dL12_dV1 = - V2p2Vbar12 .* Q1mQ2p2kVbar12 ./ V1pV2p2Vbar12_squared;
dL12_dV2 = V1 .* Q1mQ2p2kVbar12 ./ V1pV2p2Vbar12_squared;
dL21_dV1 = V2 .* Q2mQ1p2kVbar21 ./ V2pV1p2Vbar21_squared;
dL21_dV2 = - V1p2Vbar21 .* Q2mQ1p2kVbar21 ./ V2pV1p2Vbar21_squared;
dL22_dV2 = -(4 * kappa * Vbar22^2) ./ V2pV2p2Vbar22_squared;

L11mQ1 = L11 - Q1;
L12mQ1 = L12 - Q1;
L21mQ2 = L21 - Q2;
L22mQ2 = L22 - Q2;

dS11_dQ1 = 2 * ( L11mQ1 .* dL11_dQ1 - L11 );
dS12_dQ1 = 2 * ( L12mQ1 .* dL12_dQ1 - L12 );
dS12_dQ2 = 2 * L12mQ1 .* dL12_dQ2;
dS21_dQ1 = 2 * L21mQ2 .* dL21_dQ1;
dS21_dQ2 = 2 * ( L21mQ2 .* dL21_dQ2 - L21 );
dS22_dQ2 = 2 * ( L22mQ2 .* dL22_dQ2 - L22 );

dS11_dV1 = (V1_squared + V1p2Vbar11_squared) ./ V1pV1p2Vbar11_squared + 2 * L11mQ1 .* dL11_dV1;
dS12_dV1 = V2p2Vbar12_squared ./ V1pV2p2Vbar12_squared + 2 * L12mQ1 .* dL12_dV1;
dS12_dV2 = V1_squared ./ V1pV2p2Vbar12_squared + 2 * L12mQ1 .* dL12_dV2;
dS21_dV1 = V2_squared ./ V2pV1p2Vbar21_squared + 2 * L21mQ2 .* dL21_dV1;
dS21_dV2 = V1p2Vbar21_squared ./ V2pV1p2Vbar21_squared + 2 * L21mQ2 .* dL21_dV2;
dS22_dV2 = (V2_squared + V2p2Vbar22_squared) ./ V2pV2p2Vbar22_squared + 2 * L22mQ2 .* dL22_dV2;

dE1_dQ1 = ( 4 * Q_optQ1 - Q_opt_squared - 3 * V1 - 3 * Q1_squared ) / (2 * V_s); 
dE2_dQ2 = ( 4 * Q_optQ2 - Q_opt_squared - 3 * V2 - 3 * Q2_squared ) / (2 * V_s);

dE1_dV1 = ( 2 * Q_opt(:) - 3 * Q1 ) / (2 * V_s); 
dE2_dV2 = ( 2 * Q_opt(:) - 3 * Q2 ) / (2 * V_s); 

dY1_dQ1 = ( twoQ_optV1 - 6 * Q_optQ1_squared + 2 * Q_opt_squaredQ1 + 4 * Q1_cubed ) / (2 * V_s); 
dY2_dQ2 = ( twoQ_optV2 - 6 * Q_optQ2_squared + 2 * Q_opt_squaredQ2 + 4 * Q2_cubed ) / (2 * V_s); 

dY1_dV1 = ( 2 * Q_optQ1 - Q_opt_squared - 6 * V1 ) / (2 * V_s); 
dY2_dV2 = ( 2 * Q_optQ2 - Q_opt_squared - 6 * V2 ) / (2 * V_s);  


dG1_dN1 = -B1R1M11C11_K1 + dB1_dN1 .* ( R1(:) - R1M11C11N1_K1 - R1M12C12N2_K1 );
dG1_dQ1 = -B1R1C12N2_K1 .* dM12_dQ1 - Q1mQ_opt / V_s;
dG1_dV1 = -B1R1C11N1_K1 .* dM11_dV1 - B1R1M11N1_K1 .* dC11_dV1 - B1R1C12N2_K1 .* dM12_dV1 - B1R1M12N2_K1 .* dC12_dV1 - 1 / (2 * V_s);
dG1_dN2 = -B1R1M12C12_K1;
dG1_dQ2 = -B1R1C12N2_K1 .* dM12_dQ2;
dG1_dV2 = -B1R1C12N2_K1 .* dM12_dV2 - B1R1M12N2_K1 .* dC12_dV2;

dG2_dN1 = -B2R2M21C21_K2;
dG2_dQ1 = -B2R2C21N1_K2 .* dM21_dQ1;
dG2_dV1 = -B2R2C21N1_K2 .* dM21_dV1 - B2R2M21N1_K2 .* dC21_dV1;
dG2_dN2 = -B2R2M22C22_K2 + dB2_dN2 .* ( R2(:) - R2M21C21N1_K2 - R2M22C22N2_K2 );
dG2_dQ2 = -B2R2C21N1_K2 .* dM21_dQ2 - Q2mQ_opt / V_s;
dG2_dV2 = -B2R2C21N1_K2 .* dM21_dV2 - B2R2M21N1_K2 .* dC21_dV2 - B2R2C22N2_K2 .* dM22_dV2 - B2R2M22N2_K2 .* dC22_dV2 - 1 / (2 * V_s);

B1R1N1_K1dM11C11_dV1 = B1R1C11N1_K1 .* dM11_dV1 + B1R1M11N1_K1 .* dC11_dV1;
B1R1N2_K1dM12C12_dV1 = B1R1C12N2_K1 .* dM12_dV1 + B1R1M12N2_K1 .* dC12_dV1;
B1R1N2_K1dM12C12_dV2 = B1R1C12N2_K1 .* dM12_dV2 + B1R1M12N2_K1 .* dC12_dV2;
B2R2N1_K2dM21C21_dV1 = B2R2C21N1_K2 .* dM21_dV1 + B2R2M21N1_K2 .* dC21_dV1;
B2R2N1_K2dM21C21_dV2 = B2R2C21N1_K2 .* dM21_dV2 + B2R2M21N1_K2 .* dC21_dV2;
B2R2N2_K2dM22C22_dV2 = B2R2C22N2_K2 .* dM22_dV2 + B2R2M22N2_K2 .* dC22_dV2;

dH1_dN1 = ( dB1_dN1 .* R1(:) - dG1_dN1 ) .* Q1 - B1R1M11C11_K1 .* L11 + dB1_dN1 .* ( - R1M11C11N1_K1 .* L11 - R1M12C12N2_K1 .* L12 );
dH1_dQ1 = -dG1_dQ1 .* Q1 + B1R1mG1 - B1R1M11C11N1_K1 .* dL11_dQ1 - B1R1C12N2_K1 .* ( dL12_dQ1 .* M12 + L12 .* dM12_dQ1 ) + dE1_dQ1;
dH1_dV1 = -dG1_dV1 .* Q1 - B1R1M11C11N1_K1 .* dL11_dV1 - L11 .* B1R1N1_K1dM11C11_dV1 - B1R1M12C12N2_K1 .* dL12_dV1 - L12 .* B1R1N2_K1dM12C12_dV1 + dE1_dV1;
dH1_dN2 = -dG1_dN2 .* Q1 - B1R1M12C12_K1 .* L12; 
dH1_dQ2 = -dG1_dQ2 .* Q1 - B1R1C12N2_K1 .* ( dL12_dQ2 .* M12 + L12 .* dM12_dQ2 );
dH1_dV2 = -dG1_dV2 .* Q1 - B1R1M12C12N2_K1 .* dL12_dV2 - L12 .* B1R1N2_K1dM12C12_dV2;

dH2_dN1 = -dG2_dN1 .* Q2 - B2R2M21C21_K2 .* L21;
dH2_dQ1 = -dG2_dQ1 .* Q2 - B2R2C21N1_K2 .* ( dL21_dQ1 .* M21 + L21 .* dM21_dQ1 );
dH2_dV1 = -dG2_dV1 .* Q2 - B2R2M21C21N1_K2 .* dL21_dV1 - L21 .* B2R2N1_K2dM21C21_dV1;
dH2_dN2 = ( dB2_dN2 .* R2(:) - dG2_dN2 ) .* Q2 - B2R2M22C22_K2 .* L22 + dB2_dN2 .* ( - R2M21C21N1_K2 .* L21 - R2M22C22N2_K2 .* L22 );
dH2_dQ2 = -dG2_dQ2 .* Q2 + B2R2mG2 - B2R2C21N1_K2 .* ( dL21_dQ2 .* M21 + L21 .* dM21_dQ2 ) - B2R2M22C22N2_K2 .* dL22_dQ2 + dE2_dQ2;
dH2_dV2 = -dG2_dV2 .* Q2 - B2R2M21C21N1_K2 .* dL21_dV2 - L21 .* B2R2N1_K2dM21C21_dV2 - B2R2M22C22N2_K2 .* dL22_dV2 - L22 .*B2R2N2_K2dM22C22_dV2 + dE2_dV2;

dW1_dN1 = ( dB1_dN1 .* R1(:) - dG1_dN1 ) .* V1mQ1_squared - B1R1M11C11_K1 .* S11 + dB1_dN1 .* ( - R1M11C11N1_K1 .* S11 - R1M12C12N2_K1 .* S12 );
dW1_dQ1 = -dG1_dQ1 .* V1mQ1_squared - 2 * Q1 .* B1R1mG1 - B1R1M11C11N1_K1 .* dS11_dQ1 - B1R1C12N2_K1 .* ( dS12_dQ1 .* M12 + S12 .* dM12_dQ1 ) + dY1_dQ1;
dW1_dV1 = -dG1_dV1 .* V1mQ1_squared + B1R1mG1 - B1R1M11C11N1_K1 .* dS11_dV1 - S11 .* B1R1N1_K1dM11C11_dV1 - B1R1M12C12N2_K1 .* dS12_dV1 - S12 .* B1R1N2_K1dM12C12_dV1 + dY1_dV1;
dW1_dN2 = -dG1_dN2 .* V1mQ1_squared - B1R1M12C12_K1 .* S12; 
dW1_dQ2 = -dG1_dQ2 .* V1mQ1_squared - B1R1C12N2_K1 .* ( dS12_dQ2 .* M12 + S12 .* dM12_dQ2 );
dW1_dV2 = -dG1_dV2 .* V1mQ1_squared - B1R1M12C12N2_K1 .* dS12_dV2 - S12 .* B1R1N2_K1dM12C12_dV2;

dW2_dN1 = -dG2_dN1 .* V2mQ2_squared - B2R2M21C21_K2 .* S21;
dW2_dQ1 = -dG2_dQ1 .* V2mQ2_squared - B2R2C21N1_K2 .* ( dS21_dQ1 .* M21 + S21 .* dM21_dQ1 );
dW2_dV1 = -dG2_dV1 .* V2mQ2_squared - B2R2M21C21N1_K2 .* dS21_dV1 - S21 .* B2R2N1_K2dM21C21_dV1;
dW2_dN2 = ( dB2_dN2 .* R2(:) - dG2_dN2 ) .* V2mQ2_squared - B2R2M22C22_K2 .* S22 + dB2_dN2 .* ( - R2M21C21N1_K2 .* S21 - R2M22C22N2_K2 .* S22 );
dW2_dQ2 = -dG2_dQ2 .* V2mQ2_squared - 2 * Q2 .* B2R2mG2 - B2R2C21N1_K2 .* ( dS21_dQ2 .* M21 + S21 .* dM21_dQ2 ) - B2R2M22C22N2_K2 .* dS22_dQ2 + dY2_dQ2;
dW2_dV2 = -dG2_dV2 .* V2mQ2_squared + B2R2mG2 - B2R2M21C21N1_K2 .* dS21_dV2 - S21 .* B2R2N1_K2dM21C21_dV2 - B2R2M22C22N2_K2 .* dS22_dV2 - S22 .* B2R2N2_K2dM22C22_dV2 + dY2_dV2;

rows = zeros(numNonZero, 1); % row indices of nonzero elements of A
columns = zeros(numNonZero, 1); % column indicis of nonzero elements of A
elements = zeros(numNonZero, 1); % nonzero elements of A
segmentCounter = 0;

%---constructing non-zero elements at 1st column of each block--------------------------------------
currentColumn = variableLocation(1,:);

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = G1 + N1 .* dG1_dN1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH1_dN1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW1_dN1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N2 .* dG2_dN1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH2_dN1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW2_dN1;
segmentCounter = segmentCounter + 1;

%---constructing non-zero elements at 2nd column of each block--------------------------------------
currentColumn = variableLocation(2,:);

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N1 .* dG1_dQ1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH1_dQ1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW1_dQ1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N2 .* dG2_dQ1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH2_dQ1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW2_dQ1;
segmentCounter = segmentCounter + 1;

%---constructing non-zero elements at 3rd column of each block--------------------------------------
currentColumn = variableLocation(3,:);

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N1 .* dG1_dV1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH1_dV1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW1_dV1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N2 .* dG2_dV1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH2_dV1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW2_dV1;
segmentCounter = segmentCounter + 1;

%---constructing non-zero elements at 4th column of each block--------------------------------------
currentColumn = variableLocation(4,:);

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N1 .* dG1_dN2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH1_dN2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW1_dN2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = G2 + N2 .* dG2_dN2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH2_dN2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW2_dN2;
segmentCounter = segmentCounter + 1;

%---constructing non-zero elements at 5th column of each block--------------------------------------
currentColumn = variableLocation(5,:);

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N1 .* dG1_dQ2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH1_dQ2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW1_dQ2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N2 .* dG2_dQ2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH2_dQ2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW2_dQ2;
segmentCounter = segmentCounter + 1;

%---constructing non-zero elements at 6th column of each block--------------------------------------
currentColumn = variableLocation(6,:);

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N1 .* dG1_dV2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH1_dV2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW1_dV2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N2 .* dG2_dV2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH2_dV2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW2_dV2;
segmentCounter = segmentCounter + 1;

%---constructing sparse A--------------------------------------------------------------------
A = sparse(rows, columns, elements, odeSystemSize, odeSystemSize);

end
