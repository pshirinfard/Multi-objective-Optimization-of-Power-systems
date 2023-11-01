clc
%% Separating unfeasible pop from ParetoFront
%% Constants Variables
p_pv=[0 0 0 0 0 0 0.873262749000000 10.4705915900000 23.8948660100000 35.1873519600000 43.6959561700000 48.2092406000000 50 49.5876259200000 45.7649610500000 39.0085899400000 29.7993778700000 17.5194771800000 5.43648868500000 0.519391570000000 0 0 0 0];
p_wt=[20.9150326800000 19.6078431400000 13.5947712400000 9.41176470600000 11.7647058800000 5.22875817000000 2.61437908500000 3.13725490200000 13.3333333300000 21.6993464100000 19.3464052300000 19.0849673200000 23.0065359500000 29.2810457500000 37.9084967300000 38.6928104600000 40 37.3856209200000 19.3464052300000 7.58169934600000 5.39869281000000 1.04575163400000 0.784313725000000 0.672580000000000];
p_load=[17.25126706 17.75509512 10.00595705 7.199076723 10.39090325 6.459291425 4.04721553 12.63862444 22.64929401 51.75750825 65.34071726 75.75428218 77.31428483 77.74352497 80.42926362 72.84998377 60.69659525 49.75108386 22.11604468 7.973288849 7.151659507 6.330828628 2.085264091 1.2563]';

iteration = 24    ;
        p_pv   =  p_pv(iteration) ;
        p_wt   =  p_wt(iteration) ;
        p_load =  p_load(iteration);
        soc= 0.5 ;
p=1;
for k1=1:pop_size% bedune takhatiharo joda mikone
         if p_load>=p_pv+p_wt % Discharge
             
           ps_wt    = x(k1,1);
           p_dis    = x(k1,2);

            
           ps_pv(k1)=p_load-(ps_wt+p_dis);
           
           co(k1,1) = ps_wt-p_wt;
           co(k1,2) = ps_pv(k1)-p_pv;
           co(k1,3) = p_dis-soc;
           co(k1,4) = -ps_pv(k1);

         else % Charge
             
           ps_wt    = x(k1,1);
           ps_wt2es = x(k1,2);
           ps_pv2es = x(k1,3);
           
           ps_pv(k1)=p_load-(ps_wt);
           
           co(k1,1) = ps_wt2es+ps_wt-p_wt; 
           co(k1,2) = ps_pv2es+ps_pv(k1)-p_pv;   
           p_charge(k1) = ps_wt2es+ps_pv2es;
           co(k1,3) = soc+p_charge(k1)-10;
           co(k1,4) = -(ps_wt2es+ps_wt-p_wt);
           co(k1,5) = -(ps_pv2es+ps_pv(k1)-p_pv);   
           co(k1,6) = -ps_pv(k1);
         end
         
 


%      tol=1;
%         if ( (co(k1,1)<= tol) &&(co(k1,2)<= tol) && (co(k1,3)<= tol) &&...
%                 (co(k1,4)<= tol) && (co(k1,5)<= tol))
%             j(p)=k1;
%             p=p+1;
%         end
    
end

%% Use Utopia Point to find Best solution
for i = 1:size(co, 1)
    max_violation(i) = max(co(i, :));
end

[best_sol,best_pop]=min(max_violation);

% size_j = numel(j);
% 
% f_utopia(1,1:size_j)=0;% preallocation
% for o = 1:size_j% fasele taa noghte ideal o hesab mikone
%     f_utopia(o)=sqrt((ff(j(o),1)+10^3)^2+(ff(j(o),2))^2);
% end
% 
% [best_sol,best_pop]=min(f_utopia);
% Min ro joda mikone ham faselasho mige ham shomarasho tu j
% L=j(best_pop);% L shomareye behtarin ozve jamiate

     if p_load>=p_pv+p_wt % Discharge
             
           ps_wt    = x(best_pop,1);
           p_dis   = x(best_pop,2);

           ps_wt2es = 0;
           ps_pv2es = 0;
           soc=soc-p_dis+ps_wt2es+ps_pv2es;
vertcat([ps_wt ps_pv(best_pop) p_dis ps_wt2es ps_pv2es soc],[ff(best_pop,1:2) 0 0 0 0],[co(best_pop,1:4) 0 0]);

     else % Charge
             
           ps_wt    = x(best_pop,1);
           ps_wt2es = x(best_pop,2);
           ps_pv2es = x(best_pop,3);
           p_dis    = 0;
           soc=soc-p_dis+ps_wt2es+ps_pv2es;

vertcat([ps_wt ps_pv(best_pop) p_dis ps_wt2es ps_pv2es soc],[ff(best_pop,1:2) 0 0 0 0],[co(best_pop,1:6)]);

     end

% Final_Solutions=[ps_wt pwt_es ppv_es ps_wt2es ps_pv2es soc;ff(best_pop) 0 0 0 0;0 co(best_pop)];


