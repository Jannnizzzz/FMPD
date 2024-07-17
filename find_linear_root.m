
function [retval] = find_linear_root(x, y)

% function has duplicate x values, no root
if x(1) == x(2)
    
    retval = 0;
    
elseif y(1) == y(2)
    
    if y(1)*y(2) == 0.0
        retval = x(1);
        % duplicate y values in root function
    else
        retval = 0.5*(x(1)+x(2));
    end
    
else
    retval = -y(1)*(x(2)-x(1))/(y(2)-y(1)) + x(1);
end

return