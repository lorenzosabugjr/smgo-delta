function f = fn_styblinski( x_e )
    % Styblinski-Tang function
    % https://en.wikipedia.org/wiki/Test_functions_for_optimization
    % f(x*) = -39.16616 * n at x* = ( -2.903534, -2.903534, ... )    

    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end

    f = sum( x .^ 4 - 16 * x .^ 2 + 5 * x ) / 2 + 2*(rand(1) - 0.5);
    
end