function f = real2normd( x, bounds )
    f = ( x -  + bounds( : , 1 ) ) ./ ( bounds( : , 2 ) - bounds( : , 1 ) );
end