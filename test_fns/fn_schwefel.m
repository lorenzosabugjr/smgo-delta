function f = fn_schwefel( x_e )
    % Generalized Schwefel function
    % https://www.al-roomi.org/benchmarks/unconstrained/n-dimensions/176-generalized-schwefel-s-problem-2-26
    % f(x*) = ?418.982887272433799807913601398 * n at x* = ( 420.968746, 420.968746, ... )
    
    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end
    
    f = - sum( x .* sin( sqrt( abs( x ) ) ) );

end