function result =  makeAlphas( vec )
    for i = 1 : size( vec, 1 )
        result( i ) = rand();
    end
    result = result';
end