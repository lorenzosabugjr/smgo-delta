function f = fn_rosenbrock( x_e )
    % Rosenbrock function
    % https://en.wikipedia.org/wiki/Test_functions_for_optimization
    % f(x*) = 0 at x* = ( 1, 1, ... )
    
    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end
    
    f = sum( 100 * ( x( 2:end ) - x( 1:end-1 ) .^ 2 ) .^ 2 + ( 1 - x( 1:end-1 ) ) .^ 2 );

end