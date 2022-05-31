function f = fn_g23( x )
    f = -9 * x( 5 ) - 15 * x( 8 ) + 6 * x( 1 ) + 16 * x( 2 ) + 10 * ( x( 6 ) + x( 7 ) );
    g1 = x( 9 ) * x( 3 ) + 0.02 * x( 6 ) - 0.025 * x( 5 );
    g2 = x( 9 ) * x( 4 ) + 0.02 * x( 7 ) - 0.015 * x( 8 );
    f = [ f -g1 -g2 ];
end