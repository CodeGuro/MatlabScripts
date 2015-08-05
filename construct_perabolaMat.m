function peraMat = construct_perabolaMat( matdeltaX, n, p )
    
    rowCount = n*(n+1)/2;
    columnCount = n;
    peraMat = nan( rowCount, columnCount );

    for segment = 1:n
        i_para = rowCount - (n-(segment-1))*(n-(segment-2))/2 + 1;
        for row = segment:n
            for column = 1:n
                j_para = column;
                peraMat( i_para, j_para ) = matdeltaX( segment, column ) * matdeltaX( row, column );
            end
            i_para = i_para + 1;
        end
    end
    
end

