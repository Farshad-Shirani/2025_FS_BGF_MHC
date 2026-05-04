function y = smoothStep(x, x_begin, x_end, smoothness)
% a step function with smooth transition
% y changes from 0 to 1, smoothly, over the interval [x_begin, x_end]
% smoothness paramter, typically 0.5, determines how sharp are the corners (angels) at x_begin and x_end 

y = log( 1 + exp(smoothness * (x - x_begin)) ) - log( 1 + exp(smoothness * (x - x_end)) );

y = y / ( (x_end - x_begin) * smoothness );

end

