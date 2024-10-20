function val=findEvent(e_List_From,e_List_To,neuronFrom,neuronTo)
    val=find((e_List_From==neuronFrom).*(e_List_To==neuronTo));
end