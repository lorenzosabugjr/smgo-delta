function f = fn_g12( x )
    f = -( 100 - ( x( 1 ) - 5 ) ^ 2 - ( x( 2 ) - 5 ) ^ 2 - ( x( 3 ) - 5 ) ^ 2 ) / 100;
    g1 = ( x( 1 ) - s2g( x( 1 ) ) ) ^ 2 + ( x( 2 ) - s2g( x( 2 ) ) ) ^ 2 + ( x( 3 ) - s2g( x( 3 ) ) ) ^ 2 - 0.0625;
    f = [ f -g1 ];
end

function sn = s2g( x )
    if round( x ) == 10
        sn = 9;
    elseif round( x ) == 0
        sn = 1;
    else
        sn = round( x );
    end
end