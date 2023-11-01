clear all
clc
global V M xl xu etac etam p pop_size pm


%% code starts
M=2;
p=input('Test problem index  :');
pop_size=400;           % Population size
no_runs=1;              % Number of runs
gen_max=200;            % MAx number of generations - stopping criteria
fname='test_case';      % Objective function and constraint evaluation

if p==13,  % OSY
    pop_size=100; 
    no_runs=10;
end;                   
if (p==2 || p==5 || p==7), gen_max=1000; end;

if p<=9     % Unconstrained test functions
tV=[2;30;3;1;30;4;30;10;10];
V=tV(p);
txl=[-5*ones(1,V);zeros(1,V);-5*ones(1,V);-1000*ones(1,V);zeros(1,V);-1/sqrt(V)*ones(1,V);zeros(1,V); 0 -5*ones(1,V-1);zeros(1,V)]; 
txu=[10*ones(1,V); ones(1,V);5*ones(1,V);1000*ones(1,V);ones(1,V);1/sqrt(V) *ones(1,V);ones(1,V);1 5*ones(1,V-1);ones(1,V)];
xl=(txl(p,1:V));            % lower bound vector
xu=(txu(p,1:V));            % upper bound vectorfor 
etac = 20;                  % distribution index for crossover
etam = 20;                  % distribution index for mutation / mutation constant
else        % Constrained test functions
p1=p-9;
tV=[2;2;2;6;2];
% V=tV(p1);

txl=[0 0 0 0 0 0;-20 -20 0 0 0 0;0 0 0 0 0 0;0 0 1 0 1 0;0.1 0 0 0 0 0]; 
txu=[5 3 0 0 0 0;20 20 0 0 0 0;pi pi 0 0 0 0;10 10 5 6 5 10;1 5 0 0 0 0];
% xl=(txl(p1,1:V));           % lower bound vector
p_pv=[0 0 0 0 0 0 0.873262749000000 10.4705915900000 23.8948660100000 35.1873519600000 43.6959561700000 48.2092406000000 50 49.5876259200000 45.7649610500000 39.0085899400000 29.7993778700000 17.5194771800000 5.43648868500000 0.519391570000000 0 0 0 0];
p_wt=[20.9150326800000 19.6078431400000 13.5947712400000 9.41176470600000 11.7647058800000 5.22875817000000 2.61437908500000 3.13725490200000 13.3333333300000 21.6993464100000 19.3464052300000 19.0849673200000 23.0065359500000 29.2810457500000 37.9084967300000 38.6928104600000 40 37.3856209200000 19.3464052300000 7.58169934600000 5.39869281000000 1.04575163400000 0.784313725000000 0.672580000000000];
p_load=[17.25126706 17.75509512 10.00595705 7.199076723 10.39090325 6.459291425 4.04721553 12.63862444 22.64929401 51.75750825 65.34071726 75.75428218 77.31428483 77.74352497 80.42926362 72.84998377 60.69659525 49.75108386 22.11604468 7.973288849 7.151659507 6.330828628 2.085264091 1.2563]';

iteration = 12    ;
        p_pv   =  p_pv(iteration) ;
        p_wt   =  p_wt(iteration) ;
        p_load =  p_load(iteration);
% Constants Variables
if p_load>=p_pv+p_wt
        V = 2;
        b(1)=p_wt;
        b(2)=10;
else
        V = 3;
        b(1)=p_wt;
        b(2)=10;
        b(3)=10;        
end

xl=zeros(1,V);
xu=b;

etac = 20;                  % distribution index for crossover
etam = 100;                 % distribution index for mutation / mutation constant
end
pm=1/V;                     % Mutation Probability

Q=[];
for run = 1:no_runs 
%% Initial population 
xl_temp=repmat(xl, pop_size,1);
xu_temp=repmat(xu, pop_size,1);
x = xl_temp+((xu_temp-xl_temp).*rand(pop_size,V));

%% Evaluate objective function
for i =1:pop_size
[ff(i,:) err(i,:)] =feval(fname, x(i,:));           % Objective function evaulation 
end
error_norm=normalisation(err);                      % Normalisation of the constraint violation
population_init=[x ff error_norm];
[population front]=NDS_CD_cons(population_init);    % Non domination Sorting on initial population
    
%% Generation Starts
for gen_count=1:gen_max
% selection (Parent Pt of 'N' pop size)
parent_selected=tour_selection(population);                     % 10 Tournament selection
%% Reproduction (Offspring Qt of 'N' pop size)
child_offspring  = genetic_operator(parent_selected(:,1:V));    % SBX crossover and polynomial mutation

for ii = 1:pop_size
[fff(ii,:) err(ii,:)]=feval(fname, child_offspring(ii,:));      % objective function evaluation for offspring
end

error_norm=normalisation(err);                                  
child_offspring=[child_offspring fff error_norm];

%% INtermediate population (Rt= Pt U Qt of 2N size)
population_inter=[population(:,1:V+M+1) ; child_offspring(:,1:V+M+1)];
[population_inter_sorted front]=NDS_CD_cons(population_inter);              % Non domination Sorting on offspring
%% Replacement - N
new_pop=replacement(population_inter_sorted, front);
population=new_pop;
end
new_pop=sortrows(new_pop,V+1);
paretoset(run).trial=new_pop(:,1:V+M+1);
Q = [Q; paretoset(run).trial];                      % Combining Pareto solutions obtained in each run
end

%% Result and Pareto plot
if run==1
plot(new_pop(:,V+1),new_pop(:,V+2),'*')
else                                        
[pareto_filter front]=NDS_CD_cons(Q);               % Applying non domination sorting on the combined Pareto solution set
rank1_index=find(pareto_filter(:,V+M+2)==1);        % Filtering the best solutions of rank 1 Pareto
pareto_rank1=pareto_filter(rank1_index,1:V+M)
plot(pareto_rank1(:,V+1),pareto_rank1(:,V+2),'*')   % Final Pareto plot
end
xlabel('objective function 1')
ylabel('objective function 2')
if p==1
    title(' 1 - Test case 1')
elseif p==2
    title(' 2  - ZDT1')
elseif p==3 
    title(' 3  - KUR')
elseif p==4
    title(' 4  - SCH')
elseif p==5
    title(' 5  - ZDT2')
elseif p==6
    title(' 6  - Test case 3')
elseif p==7
    title(' 7  - ZDT3')
elseif p==8
    title(' 8  - ZDT4')
elseif p==9
    title(' 9  - ZDT6')
elseif p==10
    title(' 10 - BNH')
elseif p==11
    title(' 11 - SRN')
elseif p==12
    title(' 12 - TNK')
elseif p==13
    title(' 13 - OSY')
elseif p==14
    title(' 14 - CONSTR')
end