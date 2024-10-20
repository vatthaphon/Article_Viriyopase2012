function w=check(sp_Arr, n_Spike, w)
T_0=25;

eta=10^(-3);
w_out=-1;
w_in=1;
tau_syn=6/T_0;
tau_p=1/T_0;
tau_m=20/T_0;
A_p=1.5;
A_m=-1.5;

sp_Arr_i=1;
n_Spike_i=1;

is_sp_Arr_NaN=isnan(sp_Arr(sp_Arr_i, 1));
is_n_Spike_NaN=isnan(n_Spike(n_Spike_i, 1));

while ~((is_sp_Arr_NaN==1) && (is_n_Spike_NaN==1))
    if ~(is_sp_Arr_NaN==1)
        n_Spike_before_ft=(n_Spike<=sp_Arr(sp_Arr_i, 1));
        
        delta_t=sp_Arr(sp_Arr_i, 1)-n_Spike(n_Spike_before_ft, 1);
        
        w=w+eta*w_in+eta*Sum_F_PRE_EVENT(delta_t, tau_p, tau_m, A_p, A_m);
        
        sp_Arr_i=sp_Arr_i+1;
        is_sp_Arr_NaN=isnan(sp_Arr(sp_Arr_i, 1));
    end
    
    if ~(is_n_Spike_NaN==1)
        sp_Arr_before_ft=(sp_Arr<=n_Spike(n_Spike_i, 1));  
        
        delta_t=n_Spike(n_Spike_i, 1)-sp_Arr(sp_Arr_before_ft, 1);
        
        w=w+eta*w_out+eta*Sum_F_POST_EVENT(-delta_t, tau_syn, tau_p, tau_m, A_p, A_m);        
        
        n_Spike_i=n_Spike_i+1;
        is_n_Spike_NaN=isnan(n_Spike(n_Spike_i, 1));
    end    
end

end

function val=Sum_F_PRE_EVENT(delta_t, tau_p, tau_m, A_p, A_m)

    val=sum(A_p.*exp(-delta_t./tau_p)+A_m.*exp(-delta_t./tau_m));            
end

function val=Sum_F_POST_EVENT(delta_t, tau_syn, tau_p, tau_m, A_p, A_m)

    tau_p_tilde=tau_syn*tau_p/(tau_syn+tau_p);
    tau_m_tilde=tau_syn*tau_m/(tau_syn+tau_m);
    
    val=sum((exp(delta_t./tau_syn).*(A_p.*(1-delta_t./tau_p_tilde)+A_m.*(1-delta_t./tau_m_tilde))));      
end
