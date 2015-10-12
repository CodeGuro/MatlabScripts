function result = insert( vector, vals, idx )
    c=false(1,length(vector)+length(idx));
    c(idx)=true;
    result=nan(size(c));
    result(~c)=vector;
    result(c)=vals;
end

