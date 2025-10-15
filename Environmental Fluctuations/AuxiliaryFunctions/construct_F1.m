function F = construct_F1(solution, variableLocation, Dx)
%--- This function returnes the value of the F1 component of ODE vector field. The ordering is in x1 direction.

global D1 D2 V_u1 V_u2 Q_opt A_max1 A_max2 dQ_tilde  epsilon
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

F = d2S(:);
F(variableLocation(1,:)) = F(variableLocation(1,:)) * sigma1_squared   - dANQd1;
F(variableLocation(2,:)) = ( F(variableLocation(2,:)) + 2 * dN1_N1 .* dQ1 ) * sigma1_squared ...
    - dAVd1 - dN1_N1 .* AVd1 - dQ1 .* AQd1;
F(variableLocation(3,:)) = ( F(variableLocation(3,:)) + 2 * dN1_N1 .* dV1 + 2 * dQ1.^2 ) * sigma1_squared ...
    - 2 * dQ1 .* AVd1 - dV1 .* AQd1;
F(variableLocation(4,:)) = F(variableLocation(4,:)) * sigma2_squared   - dANQd2;
F(variableLocation(5,:)) = ( F(variableLocation(5,:)) + 2 * dN2_N2 .* dQ2 ) * sigma2_squared ...
    - dAVd2 - dN2_N2 .* AVd2 - dQ2 .* AQd2;
F(variableLocation(6,:)) = ( F(variableLocation(6,:)) + 2 * dN2_N2 .* dV2 + 2 * dQ2.^2 ) * sigma2_squared ...
    - 2 * dQ2 .* AVd2 - dV2 .* AQd2;

end