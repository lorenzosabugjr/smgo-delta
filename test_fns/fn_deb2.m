function f = fn_deb2( x_e )
    % Deb function #2
    % https://al-roomi.org/benchmarks/unconstrained/n-dimensions/232-deb-s-function-no-02
    % f(x*) = -1 * n at x* = <many>  
  
    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end

    f = -1 / length(x) * sum( ( sin( 5*pi*(x .^ 0.75 - 0.05 ) ) ).^ 6 );

end