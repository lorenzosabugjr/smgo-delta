function f = fn_deb1( x_e )
    % Deb function #1
    % https://www.al-roomi.org/benchmarks/unconstrained/n-dimensions/231-deb-s-function-no-01
    % f(x*) = -1 * n at x* = <many>    

    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end
    
    f = -1 / length(x) * sum( ( sin( 5*pi*x ) ).^ 6 );
    
end