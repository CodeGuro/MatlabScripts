function steady_vecx = findSteady_lin( n, matK, vecb )
    steady_vecx = ( eye( n ) - matK ) \ vecb;
end

