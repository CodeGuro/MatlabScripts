function result = epsilonfunc( i, j, vecX, matK, matN )
    result = power( vecX( j ) / matK( i, j ), matN( i, j ) );
end

