function f = fn_parabplane( x_e )
    % simple parabolic function 

    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end

    f = max( 4 - sum( ( x + 3.0 ) .^ 2 ), sum( ( x - [ -5 -5 ]' ) .^ 2 ) - 64 ) + 0.5*(rand(1) - 0.5);
end