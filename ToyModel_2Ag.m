clear all;
%close all;
clc;

%%case with 2 Ags
n_Ag = 2;
n_founders = 3;
rep = 9;
n_max_Bcells = n_founders*2^rep;
n_cycle_max = 100;
n_trial_max = 1000;
a_act = 10;
a_threshold = 20;
a_min = -6;
conc = 1.5;

p_mut = 0.2; % per division.
p_CDR = 0.3;
p_FR_lethal = 0.9;
p_recycle = 0.85;
t_cell_selection = 0.6;

founders = rand(1,n_founders);
founder_b_cells = zeros(n_Ag, n_founders);

b_cells = zeros(n_trial_max, n_Ag, n_max_Bcells);
exit_cells = zeros(n_trial_max, n_cycle_max, n_Ag, floor(n_max_Bcells/4));
number_recycled_b_cells = zeros(n_trial_max, n_cycle_max);
number_exit_cells = zeros(n_trial_max, n_cycle_max);

tic;

%%INITIALIZATION + proliferation: 
%% each founder needs to meet a_act for one Ag
cycle_number = 1;
for f = 1:n_founders
    rand_f = rand;
    %in case of 2 Ags if <0.5 then meets activation for Ag_1
        if founders(1,f) < 0.5
            founder_b_cells(1,f) = a_act;
            founder_b_cells(2,f) = a_min + (a_act - a_min)*rand_f;
        else
            founder_b_cells(2,f) = a_act;
            founder_b_cells(1,f) = a_min + (a_act - a_min)*rand_f;
        end
end   
 
number_recycled_b_cells(:,cycle_number) = size(founder_b_cells,2);
number_exit_cells(:,cycle_number) = 0;

cycle_number = cycle_number +1;

for f = 1:n_founders
    f_start = (f-1)*2^rep+1;   
    %in case of 2 Ags if <0.5 then meets activation for Ag_1
    for b = f_start:f_start+2^rep-1
        for n = 1:n_trial_max
            for i=1:n_Ag
                b_cells(n,i,b) = founder_b_cells(i,f);
            end
        end
    end
end
number_recycled_b_cells(:,cycle_number) = size(b_cells,3);
number_exit_cells(:,cycle_number) = 0;

cycle_number = cycle_number +1;

%%%% The GC starts at cycle 3

%% Stochastic process: reproduce the GC reaction many times.
% b_cells_trial = zeros(1, n_Ag, n_max_Bcells);
% exit_cells_trial = zeros(1, n_cycle_max, n_Ag, floor(n_max_Bcells/4));
% number_recycled_b_cells_trial = zeros(1, n_cycle_max);
% number_exit_cells_trial = zeros(1, n_cycle_max);

trial_number = 1;
while trial_number <= n_trial_max
    disp(['TRIAL NUMBER ' num2str(trial_number)]);
    b_cells_trial = b_cells(trial_number,:,:);
    number_recycled_b_cells_trial = number_recycled_b_cells(trial_number,:);
    exit_cells_trial = exit_cells(trial_number, :, :,:);
    number_exit_cells_trial = number_exit_cells(trial_number,:);
    
    [ b_cells_trial, number_recycled_b_cells_trial, exit_cells_trial, number_exit_cells_trial, final_cycle ] = runTrial( b_cells_trial, exit_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, conc, a_act, a_threshold, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, n_max_Bcells, n_cycle_max );    
 
    for j = 1:size(b_cells_trial,2) %n_Ag
        for k = 1:size(b_cells_trial,3) % n_max_bcells
            b_cells(trial_number, j,k) = b_cells_trial(1,j,k);
        end
    end   

    for i = 1:final_cycle
        number_recycled_b_cells(trial_number,i) = number_recycled_b_cells_trial(i);
        number_exit_cells(trial_number,i) = number_exit_cells_trial(i);
        
        for j = 1:size(b_cells_trial,2) %n_Ag
            for k = 1:size(exit_cells_trial,4)     
                exit_cells(trial_number,i,j,k) = exit_cells_trial(1,i,j,k);
            end
        end
    end
    
    trial_number = trial_number +1;
end

toc; 
%% Analyze trials
%pop_over_time = analysis(number_recycled_b_cells, n_trial_max, p_mut, p_recycle, t_cell_selection);
%[pop_time, breadth, number_exit_cells_mean, neutralized] = analysis( number_recycled_b_cells, number_exit_cells, exit_cells, n_trial_max, a_act, n_cycle_max, p_mut, p_recycle, t_cell_selection );
[ pop_time, total_exit_cells, neutralized, breadth ] = analysis( number_recycled_b_cells, number_exit_cells, exit_cells, n_trial_max, a_act, n_cycle_max, p_mut, p_recycle, t_cell_selection, conc);



