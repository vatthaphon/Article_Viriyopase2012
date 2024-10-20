function val=isEq(val1,val2)
    if (((val2-eps)<val1) && (val1<(val2+eps)))
        val=1;
    else
        val=0;
    end
end