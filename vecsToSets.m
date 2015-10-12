function result = vecsToSets( v )
    result = {};
    for i = 1 : size( v, 1 )
        result{ i } = v( i, : );
    end
end