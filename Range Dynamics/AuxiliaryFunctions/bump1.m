function y = bump1(x)
% bump function with radius 1 in one dimensional space (peak normalized to 1)

x ( abs(x) >= 1) = NaN; % x values outside the bump's support
%y = exp(- 1 ./ (1 - x.^2) ); %peak not normalized to 1
y = exp(1 - 1 ./ (1 - x.^2) ); 

y(isnan(y)) = 0;

end