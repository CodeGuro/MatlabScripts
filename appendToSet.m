function set = appendToSet( set, data )
   
   setSize = size( set, 2 );
   for i = 1 : size( data, 2 )
       set{ setSize + i } = data{ i };
   end
end