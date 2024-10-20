function val=getEarliestESet(e_List_From, e_List_To, e_List_Int)
    [s_r,s_i]=sort(e_List_Int(:,1));
    
    rows=size(s_i,1);
    min_Val=e_List_Int(s_i(1,1),1);
    
    iSp_List=0;iSp_Arr_List=0;
    if (isEq(e_List_From(s_i(1,1),1),e_List_To(s_i(1,1),1))==1)
        iSp_List=iSp_List+1;
        tmpSp_List(iSp_List,1)=s_i(1,1);        
    else
        iSp_Arr_List=iSp_Arr_List+1;
        tmpSp_Arr_List(iSp_Arr_List,1)=s_i(1,1);
    end
    
    index=2;
    while ((index<=rows) && (isEq(min_Val,e_List_Int(s_i(index,1),1))==1))        
        if (isEq(e_List_From(s_i(index,1),1),e_List_To(s_i(index,1),1))==1)
            iSp_List=iSp_List+1;
            tmpSp_List(iSp_List,1)=s_i(index,1);
        else
            iSp_Arr_List=iSp_Arr_List+1;
            tmpSp_Arr_List(iSp_Arr_List,1)=s_i(index,1);
        end                
        index=index+1;
    end
    
    % Concat. the spike and the spike-arr. events %
    if (exist('tmpSp_List')==1)
        if (exist('tmpSp_Arr_List')==1)
            val=[tmpSp_List;tmpSp_Arr_List];
        else
            val=tmpSp_List;
        end
    else
        val=tmpSp_Arr_List;
    end        
end