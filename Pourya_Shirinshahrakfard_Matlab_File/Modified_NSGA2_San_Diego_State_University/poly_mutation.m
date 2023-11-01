function mutated_child = poly_mutation(y)
global V xl xu  pm etam
% choose the type of mutation
modified=3;
if modified==0
mu=pm;                    % Mutation Rate
sigma=0.1*(xl-xu);        % Mutation Step Size

    nVar=V;
    
    nMu=ceil(mu*nVar);

    j=randsample(nVar,nMu);
    
    mutated_child=y;
    
    mutated_child(j)=y(j)+mean(sigma)*randn(size(j));

elseif modified==1
%% Polynomial mutation including boundary constraint
del=min((y-xl),(xu-y))./(xu-xl);%Always Less than 0.5 depend on location of zhen
t=rand(1,V);
loc_mut=t<pm;        
u=rand(1,V);
delq=(u<=0.5).*((((2*u)+((1-2*u).*((1-del).^(etam+1)))).^(1/(etam+1)))-1)+(u>0.5).*(1-((2*(1-u))+(2*(u-0.5).*((1-del).^(etam+1)))).^(1/(etam+1)));
c=y+delq.*loc_mut.*(xu-xl);
mutated_child=c;
elseif modified==2
   %% Description
% 1. Input is the crossovered child of size (1,V) in the vector 'y' from 'genetic_operator.m'.
% 2. Output is in the vector 'mutated_child' of size (1,V).
%% Polynomial mutation including boundary constraint
del1=(y-xl)/(xu-xl);
del2=(xu-y)/(xu-xl);

t=rand(1,V);
loc_mut=t<pm;        
u=rand(1,V);
delq=(u<=0.5).*((((2*u)+((1-2*u).*((1-del1).^(etam+1)))).^(1/(etam+1)))-1)+(u>0.5).*(1-((2*(1-u))+(2*(u-0.5).*((1-del2).^(etam+1)))).^(1/(etam+1)));
c=y+delq.*loc_mut.*(xu-xl);
mutated_child=c;    


elseif modified==3
t=rand(1,V);
loc_mut=t<pm;

delta1=(y-xl)/(xu-xl);
delta2=(xu-y)/(xu-xl);

        
u=rand(1,V);
if u>pm
  delta=min(delta1,delta2);
else
  delta=(u<0.5)*delta1+(u>=0.5)*delta2;  
end
  delq=(u<=0.5).*((((2*u)+((1-2*u).*((1-delta).^(etam+1)))).^(1/(etam+1)))-1)+(u>0.5).*(1-((2*(1-u))+(2*(u-0.5).*((1-delta).^(etam+1)))).^(1/(etam+1)));

c=y+delq.*loc_mut.*(xu-xl);
mutated_child=c;    
end
end