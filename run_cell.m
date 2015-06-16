% 1 - Rib; 2 - Nuk; 3 - Pol
clear;
seed = 12346;
p_routing = [0.2 0.2 0.4 0.2];
p_RNA_routing = [0.7 0.3]; % RNA - asRNA
p_asRNA_start_end = [1 0]; % end+  start(die)
a = 15;
n_initial_pol = a;
n_initial_rib = a;
p_initial_pol = rand (n_initial_pol,1);
n_rib = n_initial_rib; %Ribosoma
n_rnk_nuk = a; %RNK nukleaza
n_rnk_pol = n_initial_pol; %RNK polimeraza
n_asrna = 3; %asRNA
l_pol_die_rate = 0.01;
l_nuk_die_rate = 0.009;
l_rib_die_rate = 0.004;
l_asrna_die_rate = 0.01;
l_eaten_rate = 0.001;
l_rib_serv_rate = 0.01;
g_addasrna_rate = 0.0;
l_addasrna_rate_1 = g_addasrna_rate;%Ribosome
l_addasrna_rate_2 = g_addasrna_rate; %RNA nukl
l_addasrna_rate_3 = 0.01;%0.005 %RNA poly
l_addasrna_rate_4 = g_addasrna_rate;
l_queue_time = 1;

split = 1000-1; %1000%
l_steps_x = 4; %4
l_steps_y = 4; %4
l_step_seed = 3; %4
rnk_pol_discr = zeros(l_steps_y,l_steps_x,1:split);
rib_discr = zeros(l_steps_y,l_steps_x,1:split);
result_frequency = zeros(500,400,1);
for s = 1:l_step_seed
for k = 1:l_steps_y

    for i = 1:l_steps_x
        strcat('seed=',int2str(s), ' y=',int2str(k),' x=',int2str(i))
        open('Cell_v1');
        set_param('Cell_v1', 'StopTime', int2str(split+1));
        sim('Cell_v1');
        l_size_pol = out_rnk_pol.length;
        l_size_rib = out_rib.length;
        for j = 1:split
            rnk_pol_discr(k,i,j)=GetTimeValue(j,l_size_pol,out_rnk_pol.time,out_rnk_pol.data);
            rib_discr(k,i,j)=GetTimeValue(j,l_size_rib,out_rib.time,out_rib.data);
            if rib_discr(k,i,j) > 0 && rnk_pol_discr(k,i,j) > 0 
                
            result_frequency(rib_discr(k,i,j),rnk_pol_discr(k,i,j),1) = result_frequency(rib_discr(k,i,j),rnk_pol_discr(k,i,j),1)+1;
            end
        end
        n_rib = n_rib + 10;
    end
    n_rib = n_initial_rib;
    n_rnk_pol = n_rnk_pol + 10;
    p_initial_pol = rand (n_rnk_pol,1);
    
end
seed = (seed-321);
n_rnk_pol = n_initial_pol;
n_rib = n_initial_rib;
p_initial_pol = rand (n_rnk_pol,1);
end
normalization = sum(sum(result_frequency));
result_frequency = result_frequency/normalization;
surface (result_frequency(1:200, 1:100));


%plot( rib_discr(1,1,1:split),rnk_pol_discr(1,1,1:split), rib_discr(5,5,1:split),rnk_pol_discr(5,5,1:split));
 xlabel('RNA poly');
 ylabel('Ribosome');
%  hleg1 = legend('ribosome','nuklease','polymerase','queue');