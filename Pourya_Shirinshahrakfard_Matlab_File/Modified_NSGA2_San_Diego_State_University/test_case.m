function [fit, err]=test_case(x)
global p V
%% Unconstrained Test functions (for p=1 to p=9)
if p==1     % Test case problem 1
    f1=(4*x(1)^2)+(4*x(2)^2);                       
    f2=((x(1)-5)^2)+((x(2)-5)^2);               
end

if p==2     % ZDT1 from Deb paper NSGA2
    cons=[0];
    f1 = x(1);
       g=1+(9*sum(x(2:V),2)/(V-1));            
       f2 = g*(1-sqrt(x(1)/g));                  
end

if p==3     % kUR from Deb
    f1=(-10*exp(-0.2*(sqrt(x(1)^2+x(2)^2))))+(-10*exp(-0.2*(sqrt(x(2)^2+x(3)^2))));
    f2=((abs(x(1))^0.8) + (5*sin(x(1))^3))+((abs(x(2))^0.8) + (5*sin(x(2))^3))+((abs(x(3))^0.8) + (5*sin(x(3))^3));
end


if p==4    % SCH frm Deb paper
    f1=x.*x;
    f2=(x-2).^2;
end

if p==5     % ZDT2
    f1 = (x(1));
    g=1+(9*sum(x(2:V),2)/(V-1));             
    f2 =((1-(x(1)/g)^2));                
end   

if p==6     % Test case problem 2
    f1=1-exp(-sum((x-1/sqrt(V)).^2,2));
    f2=1-exp(-sum((x+1/sqrt(V)).^2,2));
end

if p==7     % ZDT3
    f1 = x(1);                                 
    g=1+(9*sum(x(2:V),2)/(V-1));               
    f2 = (1-(sqrt(x(1)/g)) - ((x(1)/g)*sin(10*pi*x(1))));               
end  

if p==8     % ZDT4       
    f1 = x(1);  temp=0;
    for ii = 2: V
        temp=temp+((x(ii))^2)-(10*cos(4*pi*x(ii)));
    end
    g= 1 + (10*(V-1)) + temp;           
    f2 = (1-sqrt(x(1)/g));                 
end  

if p==9     % ZDT6       
    f1 = 1-(exp(-4*x(1)))*(sin(6*pi*x(1)))^6; 
    g=1+(9*(sum(x(2:V),2)/(V-1))^0.25);        
    f2 = (1-(f1/g)^2);                     
end  
err= zeros(1,1);

%% Constrained Test functions (for p=10 to p=14)

if p==10     %BNH 
%% Constants Variables
p_pv=[0 0 0 0 0 0 0.873262749000000 10.4705915900000 23.8948660100000 35.1873519600000 43.6959561700000 48.2092406000000 50 49.5876259200000 45.7649610500000 39.0085899400000 29.7993778700000 17.5194771800000 5.43648868500000 0.519391570000000 0 0 0 0];
p_wt=[20.9150326800000 19.6078431400000 13.5947712400000 9.41176470600000 11.7647058800000 5.22875817000000 2.61437908500000 3.13725490200000 13.3333333300000 21.6993464100000 19.3464052300000 19.0849673200000 23.0065359500000 29.2810457500000 37.9084967300000 38.6928104600000 40 37.3856209200000 19.3464052300000 7.58169934600000 5.39869281000000 1.04575163400000 0.784313725000000 0.672580000000000];
p_load=[17.25126706 17.75509512 10.00595705 7.199076723 10.39090325 6.459291425 4.04721553 12.63862444 22.64929401 51.75750825 65.34071726 75.75428218 77.31428483 77.74352497 80.42926362 72.84998377 60.69659525 49.75108386 22.11604468 7.973288849 7.151659507 6.330828628 2.085264091 1.2563]';

iteration = 12    ;
        p_pv   =  p_pv(iteration) ;
        p_wt   =  p_wt(iteration) ;
        p_load =  p_load(iteration);
        soc    =  0.5 ;
%% Constraints and OBJ functions
         
         if p_load>=p_pv+p_wt % Discharge
             
           ps_wt    = x(1);
           p_dis   = x(2);

            
           ps_pv=p_load-(ps_wt+p_dis);
             
           f1 =20*ps_wt+25*ps_pv+30*p_dis; % Min
           f2 =-(20*ps_wt+25*ps_pv+30*p_dis);% Max
           
           c(1,1) = (ps_wt-p_wt)*((ps_wt-p_wt)>=0)+((-1)*((ps_wt-p_wt)<0));
           c(1,2) = (ps_pv-p_pv)*((ps_pv-p_pv)>=0)+((-1)*((ps_pv-p_pv)<0));
           c(1,3) = ((p_dis-soc))*(((p_dis-soc))>=0)+((-1)*((p_dis-soc)<0));
           c(1,4) = (-100*(ps_pv))*((-100*(ps_pv))>=0)+((-1)*((-100*(ps_pv))<0));
           
         else % Charge
             
           ps_wt    = x(1);
           ps_wt2es = x(2);
           ps_pv2es = x(3);
           
           ps_pv=p_load-(ps_wt);
           
           f1 =20*ps_wt+25*ps_pv+10*ps_wt2es+12.5*ps_pv2es; % Min
           f2 =-(20*ps_wt+25*ps_pv+10*ps_wt2es+12.5*ps_pv2es);% Max
           c(1,1) = (ps_wt2es+ps_wt-p_wt)*((ps_wt2es+ps_wt-p_wt)>=0)+((-1)*((ps_wt2es+ps_wt-p_wt)<0)); 
           c(1,2) = (ps_pv2es+ps_pv-p_pv)*((ps_pv2es+ps_pv-p_pv)>=0)+((-1)*((ps_pv2es+ps_pv-p_pv)<0));   
           p_charge = ps_wt2es+ps_pv2es;
           c(1,3) = 100*(soc+p_charge-10)*((soc+p_charge-10)>=0)+((-1)*((soc+p_charge-10)<0));
           c(1,4) = (-(ps_wt2es+ps_wt-p_wt)*(-(ps_wt2es+ps_wt-p_wt)>=0)+((-1)*(-(ps_wt2es+ps_wt-p_wt)<0)));
           c(1,5) = (-(ps_pv2es+ps_pv-p_pv)*(-(ps_pv2es+ps_pv-p_pv)>=0)+((-1)*(-(ps_pv2es+ps_pv-p_pv)<0)));   
           c(1,6) = (-100*(ps_pv))*((-100*(ps_pv))>=0)+((-1)*((-100*(ps_pv))<0));

         end

           
    err=(c>0).*c;
end
if p==11     %SRN  
    f1=(x(1)-2)^2+(x(2)-1)^2+2;
    f2=9*x(1)-(x(2)-1)^2;
    c(1,1)=x(1)^2+x(2)^2-225;
    c(1,2)=x(1)-(3*x(2))+10;
    err=(c>0).*c;
end
if p==12     %TNK
    f1=x(1);
    f2=x(2);
    c(1,1)=-x(1)^2-x(2)^2+1+(0.1*cos(16*atan((x(1)/x(2))))); 
    c(1,2)=(x(1)-0.5)^2+(x(2)-0.5)^2-0.5;
    err=(c>0).*c;
end

if p==13     % OSY 
    f1=-((25*(x(1)-2)^2)+((x(2)-2)^2)+((x(3)-1)^2)+((x(4)-4)^2)+((x(5)-1)^2));
    f2=(x(1)^2)+(x(2)^2)+(x(3)^2)+(x(4)^2)+(x(5)^2)+(x(6)^2);
    c(1,1)=-x(1)-x(2)+2;
    c(1,2)=-6+x(1)+x(2);
    c(1,3)=-2+x(2)-x(1);
    c(1,4)=-2+x(1)-3*x(2);
    c(1,5)=-4+((x(3)-3)^2)+x(4);
    c(1,6)=-((x(5)-3)^2)-x(6)+4;
    err=(c>0).*c;
end

if p==14    % CONSTR
    f1=x(1);
    f2=(1+x(2))/(x(1));
    c(1,1)=-x(2)-(9*x(1))+6;
    c(1,2)=+x(2)-9*x(1)+1;
    err=(c>0).*c;
end
fit=[f1 f2];

