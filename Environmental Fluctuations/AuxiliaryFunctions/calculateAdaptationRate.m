function [random,  optimal, local] = calculateAdaptationRate(solution, variableLocation, modelParameters, discretizationParamaters )
%--- This function is the same as construct_F1 and Construct_F2andA2. It just returns the two components of F1 associated
%--- with random and optimal gene flow (F1 is the some of the two components), and the reaction (local) component F2 separaterly. 

Dx = discretizationParamaters.Dx;
D1 = modelParameters.D1;
D2 = modelParameters.D2;
V_u1 = modelParameters.V_u1;
V_u2 = modelParameters.V_u2;
Q_opt = modelParameters.Q_opt;
A_max1 = modelParameters.A_max1;
A_max2 = modelParameters.A_max2;
dQ_tilde = modelParameters.dQ_tilde;
epsilon = modelParameters.epsilon;
V_s = modelParameters.V_s;
U = modelParameters.U;
kappa = modelParameters.kappa;
K1 = modelParameters.K1;
K2 = modelParameters.K2;
R1 = modelParameters.R1;
R2 = modelParameters.R2;
AlleeFlag = modelParameters.AlleeFlag;
J1 = modelParameters.J1;
J2 = modelParameters.J2;

sigma1_squared = D1;
sigma2_squared = D2;

Ad1 = (A_max1 ./ V_u1) .* dQ_tilde;
AQd1 = Ad1 .* ( solution(:,2) - Q_opt );
AVd1 = Ad1 .* solution(:,3);

Ad2 = (A_max2 ./ V_u2) .* dQ_tilde;
AQd2 = Ad2 .* ( solution(:,5) - Q_opt );
AVd2 = Ad2 .* solution(:,6);

N1 = solution(:,1);
N2 = solution(:,4);

%--------------------------------------
solution_minusTwo = reflect(solution,[2 0]); % solution(i-2)  reflect also applies the reflecting boundary condition
solution_minusOne = reflect(solution,[1 0]); % solution(i-1)  
solution_plusOne = reflect(solution,[-1 0]); % solution(i+1)  
solution_plusTwo = reflect(solution,[-2 0]); % solution(i+2)

Ad1_minusTwo = reflect(Ad1,[2 0]); % Ad(i-2)  reflect also applies the reflecting boundary condition
Ad1_minusOne = reflect(Ad1,[1 0]); % Ad(i-1)  
Ad1_plusOne = reflect(Ad1,[-1 0]); % Ad(i+1)  
Ad1_plusTwo = reflect(Ad1,[-2 0]); % Ad(i+2)

Ad2_minusTwo = reflect(Ad2,[2 0]); % Ad(i-2)  reflect also applies the reflecting boundary condition
Ad2_minusOne = reflect(Ad2,[1 0]); % Ad(i-1)  
Ad2_plusOne = reflect(Ad2,[-1 0]); % Ad(i+1)  
Ad2_plusTwo = reflect(Ad2,[-2 0]); % Ad(i+2)

AQd1_minusTwo = reflect(AQd1,[2 0]); % AQd(i-2)  reflect also applies the reflecting boundary condition
AQd1_minusOne = reflect(AQd1,[1 0]); % AQd(i-1)  
AQd1_plusOne = reflect(AQd1,[-1 0]); % AQd(i+1)  
AQd1_plusTwo = reflect(AQd1,[-2 0]); % AQd(i+2)

AQd2_minusTwo = reflect(AQd2,[2 0]); % AQd(i-2)  reflect also applies the reflecting boundary condition
AQd2_minusOne = reflect(AQd2,[1 0]); % AQd(i-1)  
AQd2_plusOne = reflect(AQd2,[-1 0]); % AQd(i+1)  
AQd2_plusTwo = reflect(AQd2,[-2 0]); % AQd(i+2)
%----------------------------------------

dS = ( -solution_plusTwo + 8 * solution_plusOne - 8 * solution_minusOne + solution_minusTwo ) / (12 * Dx);
dN1_N1 = dS(:,1) ./ (N1 + epsilon);
dQ1 = dS(:,2);
dV1 = dS(:,3);
dN2_N2 = dS(:,4) ./ (N2 + epsilon);
dQ2 = dS(:,5);
dV2 = dS(:,6);

d2S = ( -solution_plusTwo + 16 * solution_plusOne - 30 * solution + 16 * solution_minusOne - solution_minusTwo ) / (12 * Dx^2);
d2S = d2S'; % reorders 2-dim array d2S so that linear indexing d2S(:) corresponds to current ordeing in U

AVd1_minusTwo = Ad1_minusTwo .* solution_minusTwo(:,3);
AVd1_minusOne = Ad1_minusOne .* solution_minusOne(:,3);
AVd1_plusOne = Ad1_plusOne .* solution_plusOne(:,3);
AVd1_plusTwo = Ad1_plusTwo .* solution_plusTwo(:,3);
dAVd1 = ( -AVd1_plusTwo + 8 * AVd1_plusOne - 8 * AVd1_minusOne + AVd1_minusTwo ) / (12 * Dx);

AVd2_minusTwo = Ad2_minusTwo .* solution_minusTwo(:,6);
AVd2_minusOne = Ad2_minusOne .* solution_minusOne(:,6);
AVd2_plusOne = Ad2_plusOne .* solution_plusOne(:,6);
AVd2_plusTwo = Ad2_plusTwo .* solution_plusTwo(:,6);
dAVd2 = ( -AVd2_plusTwo + 8 * AVd2_plusOne - 8 * AVd2_minusOne + AVd2_minusTwo ) / (12 * Dx);

ANQd1_minusTwo = AQd1_minusTwo .* solution_minusTwo(:,1);
ANQd1_minusOne = AQd1_minusOne .* solution_minusOne(:,1);
ANQd1_plusOne = AQd1_plusOne .* solution_plusOne(:,1);
ANQd1_plusTwo = AQd1_plusTwo .* solution_plusTwo(:,1);
dANQd1 = ( -ANQd1_plusTwo + 8 * ANQd1_plusOne - 8 * ANQd1_minusOne + ANQd1_minusTwo ) / (12 * Dx);

ANQd2_minusTwo = AQd2_minusTwo .* solution_minusTwo(:,4);
ANQd2_minusOne = AQd2_minusOne .* solution_minusOne(:,4);
ANQd2_plusOne = AQd2_plusOne .* solution_plusOne(:,4);
ANQd2_plusTwo = AQd2_plusTwo .* solution_plusTwo(:,4);
dANQd2 = ( -ANQd2_plusTwo + 8 * ANQd2_plusOne - 8 * ANQd2_minusOne + ANQd2_minusTwo ) / (12 * Dx);

random = d2S(:);
random(variableLocation(1,:)) = random(variableLocation(1,:)) * sigma1_squared;
random(variableLocation(2,:)) = ( random(variableLocation(2,:)) + 2 * dN1_N1 .* dQ1 ) * sigma1_squared;
random(variableLocation(3,:)) = ( random(variableLocation(3,:)) + 2 * dN1_N1 .* dV1 + 2 * dQ1.^2 ) * sigma1_squared;
random(variableLocation(4,:)) = random(variableLocation(4,:)) * sigma2_squared;
random(variableLocation(5,:)) = ( random(variableLocation(5,:)) + 2 * dN2_N2 .* dQ2 ) * sigma2_squared;
random(variableLocation(6,:)) = ( random(variableLocation(6,:)) + 2 * dN2_N2 .* dV2 + 2 * dQ2.^2 ) * sigma2_squared;

optimal = zeros(size(random));
optimal(variableLocation(1,:)) = - dANQd1;
optimal(variableLocation(2,:)) = - dAVd1 - dN1_N1 .* AVd1 - dQ1 .* AQd1;
optimal(variableLocation(3,:)) = - 2 * dQ1 .* AVd1 - dV1 .* AQd1;
optimal(variableLocation(4,:)) = - dANQd2;
optimal(variableLocation(5,:)) = - dAVd2 - dN2_N2 .* AVd2 - dQ2 .* AQd2;
optimal(variableLocation(6,:)) = - 2 * dQ2 .* AVd2 - dV2 .* AQd2;



odeSystemSize = numel(solution);

Q1 = solution(:,2);
V1 = solution(:,3);
Q2 = solution(:,5);
V2 = solution(:,6);

%% Definign the additional growth term to model Allee effect=================
if AlleeFlag  
    expN1mJ1 = exp(-(N1 - J1) / 0.05);
    expN2mJ2 = exp(-(N2 - J2) / 0.05);
    B1 = 2.6 ./ ( 1 + expN1mJ1 ) - (2.6/2);
    B2 = 2.6 ./ ( 1 + expN2mJ2 ) - (2.6/2);

else %---Set Bi = 1 and dBi_dNi = 0 to remove the Allee effect
    B1 = 1;
    B2 = 1;
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

V1pV1p2Vbar11 = V1 + V1p2Vbar11;
V1pV2p2Vbar12 = V1 + V2p2Vbar12;
V2pV1p2Vbar21 = V2 + V1p2Vbar21;
V2pV2p2Vbar22 = V2 + V2p2Vbar22;

Q_opt_squared = Q_opt(:).^2;
Q1_squared = Q1.^2;
Q2_squared = Q2.^2;
Q1_cubed = Q1_squared .* Q1;
Q2_cubed = Q2_squared .* Q2;
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
B1R1M11C11N1_K1 = B1R1C11N1_K1 .* M11;
B1R1M12C12N2_K1 = B1R1C12N2_K1 .* M12;
B2R2M21C21N1_K2 = B2R2C21N1_K2 .* M21;
B2R2M22C22N2_K2 = B2R2C22N2_K2 .* M22;

G1 = B1R1 - B1R1M11C11N1_K1 - B1R1M12C12N2_K1 - ( Q1mQ_opt.^2 + V1 ) / (2 * V_s);
G2 = B2R2 - B2R2M21C21N1_K2 - B2R2M22C22N2_K2 - ( Q2mQ_opt.^2 + V2 ) / (2 * V_s);
B1R1mG1 =  B1R1 - G1;
B2R2mG2 =  B2R2 - G2;
H1 = B1R1mG1 .* Q1 - B1R1M11C11N1_K1 .* L11 - B1R1M12C12N2_K1 .* L12 + E1;
H2 = B2R2mG2 .* Q2 - B2R2M21C21N1_K2 .* L21 - B2R2M22C22N2_K2 .* L22 + E2;
W1 = B1R1mG1 .* V1mQ1_squared - B1R1M11C11N1_K1 .* S11 - B1R1M12C12N2_K1 .* S12 + Y1;
W2 = B2R2mG2 .* V2mQ2_squared - B2R2M21C21N1_K2 .* S21 - B2R2M22C22N2_K2 .* S22 + Y2;

%% Calculating F3 (local)==========================================================
local = zeros(odeSystemSize, 1);
local(variableLocation(1,:)) = G1 .* N1;
local(variableLocation(2,:)) = H1;
local(variableLocation(3,:)) = W1;
local(variableLocation(4,:)) = G2 .* N2;
local(variableLocation(5,:)) = H2;
local(variableLocation(6,:)) = W2;

end