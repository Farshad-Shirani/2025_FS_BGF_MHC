function A = construct_A1(solution, variableLocation, Dx)
%--- This function gives the sparse matrix resulted from discritization in x1 direction & linearization

global D1 D2 V_u1 V_u2 Q_opt A_max1 A_max2 dQ_tilde epsilon
sigma1_squared = D1;
sigma2_squared = D2;

segmentSize = numel(solution(:,1)); % number of (A)ij segments
numNonZero = 2 * (9 + 15 + 15) * segmentSize; %number of non-zero elements in A.
odeSystemSize = numel(solution);

%-----------------------------------------------
solution_minusTwo = reflect(solution,[2 0]); % solution(i-2)  reflect also applies the reflecting boundary condition
solution_minusOne = reflect(solution,[1 0]); % solution(i-1)  
solution_plusOne = reflect(solution,[-1 0]); % solution(i+1)  
solution_plusTwo = reflect(solution,[-2 0]); % solution(i+2)  
%----------------------------------------------

dS = ( -solution_plusTwo + 8*solution_plusOne - 8*solution_minusOne + solution_minusTwo ) / (12 * Dx);

%---creating linear arrays
N1 = solution(:,1);
N2 = solution(:,4);
Q1 = solution(:,2);
Q2 = solution(:,5);
V1 = solution(:,3);
V2 = solution(:,6);

dN1 = dS(:,1);
dN2 = dS(:,4);
dQ1 = dS(:,2);
dQ2 = dS(:,5);
dV1 = dS(:,3);
dV2 = dS(:,6);

dN1_N1 = dN1 ./ (N1 + epsilon);
dN2_N2 = dN2 ./ (N2 + epsilon);
dQ1_N1 = dQ1 ./ (N1 + epsilon);
dQ2_N2 = dQ2 ./ (N2 + epsilon);
dV1_N1 = dV1 ./ (N1 + epsilon);
dV2_N2 = dV2 ./ (N2 + epsilon);

Ad1 = (A_max1 ./ V_u1) .* dQ_tilde;
ANd1 = Ad1 .* N1;
AVd1 = Ad1 .* V1;
AQd1 = Ad1 .* (Q1 - Q_opt);

Ad2 = (A_max2 ./ V_u2) .* dQ_tilde;
ANd2 = Ad2 .* N2;
AVd2 = Ad2 .* V2;
AQd2 = Ad2 .* (Q2 - Q_opt);

AVd1_N1 = AVd1 ./ (N1 + epsilon);
AddN1_N1 = Ad1 .* dN1_N1;
AddQ1 = Ad1 .* dQ1;
AddV1 = Ad1 .* dV1;

AVd2_N2 = AVd2 ./ (N2 + epsilon);
AddN2_N2 = Ad2 .* dN2_N2;
AddQ2 = Ad2 .* dQ2;
AddV2 = Ad2 .* dV2;

rows = zeros(numNonZero, 1); % row indices of nonzero elements of A
columns = zeros(numNonZero, 1); % column indicis of nonzero elements of A
elements = zeros(numNonZero, 1); % nonzero elements of A

segmentCounter = 0;
%% diagonal blocks of A (center block of (A)ij)===========================================================
firstColumn = variableLocation(1,:); % indices of the first (left) columns of the center block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma1_squared;
segmentCounter = segmentCounter + 1;

%------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -2 * dN1_N1 .* dQ1_N1 * sigma1_squared   + dN1_N1 .* AVd1_N1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma1_squared   - AddQ1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = - AddN1_N1;
segmentCounter = segmentCounter + 1;

%-------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -2 * dN1_N1 .* dV1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = - AddV1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma1_squared   - 2 * AddQ1;
segmentCounter = segmentCounter + 1;

%------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma2_squared;
segmentCounter = segmentCounter + 1;

%------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) =  -2 * dN2_N2 .* dQ2_N2 * sigma2_squared  + dN2_N2 .* AVd2_N2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma2_squared   - AddQ2 ;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = - AddN2_N2;
segmentCounter = segmentCounter + 1;

%-------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = -2 * dN2_N2 .* dV2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = - AddV2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma2_squared   - 2 * AddQ2;
segmentCounter = segmentCounter + 1;


%% 1st lower diagonal of A (first left block of (A)ij)====================================================
firstColumn = reflect(variableLocation(1,:), 1); % indices of the first (left) columns of the first left block of (A)ij. reflect applies the reflecting boundary condition

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma1_squared   + 8 / (12*Dx) * AQd1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 8 / (12*Dx) * ANd1;
segmentCounter = segmentCounter + 1;

%----------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dQ1_N1 * sigma1_squared   + 8 / (12*Dx) .* AVd1_N1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN1_N1 ) * sigma1_squared ...
    + 8 / (12*Dx) * AQd1 - AddQ1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = 8 / (12*Dx) * Ad1 - AddN1_N1;
segmentCounter = segmentCounter + 1;

%-----------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dV1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * (-8) / (12*Dx) * dQ1 * sigma1_squared ...
    + 2 * 8 / (12*Dx) * AVd1 - AddV1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN1_N1 ) * sigma1_squared ...
    -2 * AddQ1 + 8 / (12*Dx) * AQd1;
segmentCounter = segmentCounter + 1;

%------------------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma2_squared    + 8 / (12*Dx) * AQd2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = 8 / (12*Dx) * ANd2;
segmentCounter = segmentCounter + 1;

%-----------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dQ2_N2 * sigma2_squared + 8 / (12*Dx) .* AVd2_N2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN2_N2 ) * sigma2_squared ...
    + 8 / (12*Dx) * AQd2 - AddQ2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = 8 / (12*Dx) * Ad2 - AddN2_N2;
segmentCounter = segmentCounter + 1;

%------------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dV2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = 4 * (-8) / (12*Dx) * dQ2 * sigma2_squared ...
    + 2 * 8 / (12*Dx) * AVd2 - AddV2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN2_N2 ) * sigma2_squared ...
    - 2 * AddQ2 + 8 / (12*Dx) * AQd2;
segmentCounter = segmentCounter + 1;


%% 2nd lower diagonal of A (secend left block of (A)ij)======================================================
firstColumn = reflect(variableLocation(1,:), 2); % indices of the first (left) columns of the second left block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma1_squared   - 1 / (12*Dx) * AQd1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx) * ANd1;
segmentCounter = segmentCounter + 1;

%---------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dQ1_N1 * sigma1_squared   - 1 / (12*Dx) * AVd1_N1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN1_N1 ) * sigma1_squared ...
    - 1 / (12*Dx) * AQd1 - AddQ1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = - 1 / (12*Dx) * Ad1 - AddN1_N1;
segmentCounter = segmentCounter + 1;

%---------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dV1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 / (12*Dx) * dQ1 * sigma1_squared   - 2 / (12*Dx)* AVd1 - AddV1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN1_N1 ) * sigma1_squared ...
    - 2 * AddQ1 - 1 / (12*Dx) * AQd1;
segmentCounter = segmentCounter + 1;

%---------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma2_squared    - 1 / (12*Dx) * AQd2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx) * ANd2;
segmentCounter = segmentCounter + 1;

%---------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dQ2_N2 * sigma2_squared   - 1 / (12*Dx) * AVd2_N2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN2_N2 ) * sigma2_squared ...
    - 1 / (12*Dx) * AQd2 - AddQ2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = - 1 / (12*Dx) * Ad2 - AddN2_N2;
segmentCounter = segmentCounter + 1;

%------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dV2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = 4 / (12*Dx) * dQ2 * sigma2_squared  - 2 / (12*Dx)* AVd2 - AddV2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN2_N2 ) * sigma2_squared ...
    - 2 * AddQ2 - 1 / (12*Dx) * AQd2;
segmentCounter = segmentCounter + 1;

%% 1st upper diagonal of A (first right block of (A)ij)======================================================
firstColumn = reflect(variableLocation(1,:), -1); % indices of the first (left) columns of the first right block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma1_squared    - 8 / (12*Dx) * AQd1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = -8 / (12*Dx) * ANd1;
segmentCounter = segmentCounter + 1;

%---------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dQ1_N1 * sigma1_squared   - 8 / (12*Dx) .* AVd1_N1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN1_N1 ) * sigma1_squared ...
    - 8 / (12*Dx) * AQd1 - AddQ1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = -8 / (12*Dx) * Ad1 - AddN1_N1;
segmentCounter = segmentCounter + 1;

%----------------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dV1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * 8 / (12*Dx) * dQ1 * sigma1_squared ...
    - 2 * 8 / (12*Dx) * AVd1 - AddV1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN1_N1 ) * sigma1_squared ...
    -2 * AddQ1 - 8 / (12*Dx) * AQd1;
segmentCounter = segmentCounter + 1;

%--------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma2_squared    - 8 / (12*Dx) * AQd2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = -8 / (12*Dx) * ANd2;
segmentCounter = segmentCounter + 1;

%---------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dQ2_N2 * sigma2_squared   - 8 / (12*Dx) .* AVd2_N2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN2_N2 ) * sigma2_squared ...
    - 8 / (12*Dx) * AQd2 - AddQ2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = -8 / (12*Dx) * Ad2 - AddN2_N2;
segmentCounter = segmentCounter + 1;

%---------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dV2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = 4 * 8 / (12*Dx) * dQ2 * sigma2_squared ...
    - 2 * 8 / (12*Dx) * AVd2 - AddV2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN2_N2 ) * sigma2_squared ...
    -2 * AddQ2 - 8 / (12*Dx) * AQd2;
segmentCounter = segmentCounter + 1;


%% 2nd upper diagonal of A (secend right block of (A)ij)=====================================================
firstColumn = reflect(variableLocation(1,:), -2); % indices of the first (left) columns of the second right block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma1_squared   + 1 / (12*Dx) * AQd1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 1 / (12*Dx) * ANd1;
segmentCounter = segmentCounter + 1;

%---------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dQ1_N1 * sigma1_squared   + 1 / (12*Dx) * AVd1_N1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 * (-1) / (12*Dx) * dN1_N1 ) * sigma1_squared ...
    + 1 / (12*Dx) * AQd1 - AddQ1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = 1 / (12*Dx) * Ad1 - AddN1_N1;
segmentCounter = segmentCounter + 1;

%-----------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dV1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * (-1) / (12*Dx) * dQ1 * sigma1_squared   + 2 / (12*Dx)* AVd1 - AddV1;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 * (-1) / (12*Dx) * dN1_N1 ) * sigma1_squared ...
    - 2 * AddQ1 + 1 / (12*Dx) * AQd1;
segmentCounter = segmentCounter + 1;

%----------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma2_squared    + 1 / (12*Dx) * AQd2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = 1 / (12*Dx) * ANd2;
segmentCounter = segmentCounter + 1;

%-----------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dQ2_N2 * sigma2_squared   + 1 / (12*Dx) * AVd2_N2;
segmentCounter = segmentCounter + 1;


segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 * (-1) / (12*Dx) * dN2_N2 ) * sigma2_squared ...
    + 1 / (12*Dx) * AQd2 - AddQ2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = 1 / (12*Dx) * Ad2 - AddN2_N2;
segmentCounter = segmentCounter + 1;

%-----------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dV2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = 4 * (-1) / (12*Dx) * dQ2 * sigma2_squared   + 2 / (12*Dx)* AVd2 - AddV2;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2  * (-1) / (12*Dx) * dN2_N2 ) * sigma2_squared ...
    - 2 * AddQ2 + 1 / (12*Dx) * AQd2;
segmentCounter = segmentCounter + 1;


%% Constructing sparse A ==================================================================================

A = sparse(rows, columns, elements, odeSystemSize, odeSystemSize);

end

