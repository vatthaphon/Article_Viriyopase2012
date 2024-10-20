function addEvent(fromNeuron, toNeuron, int)
global i_e_List
global e_List_From e_List_To e_List_Int

if (fromNeuron~=toNeuron)
    i_e_List=i_e_List+1;
    e_List_From(i_e_List,1)=fromNeuron;
    e_List_To(i_e_List,1)=toNeuron;
    e_List_Int(i_e_List,1)=int;
else
    % Update the natual-spike 1st neuron in the list %
    i_neuron=findEvent(e_List_From,e_List_To,fromNeuron,toNeuron);
    cnt=size(i_neuron,1);
    
    % There exist something to be updated %
    if (cnt==1)
        e_List_Int(i_neuron)=int;
    elseif (cnt==0)
        i_e_List=i_e_List+1;
        e_List_From=[fromNeuron;e_List_From];
        e_List_To=[toNeuron;e_List_To];
        e_List_Int=[int;e_List_Int];
    else
        disp('addEvent:Strange!');
    end
end

end