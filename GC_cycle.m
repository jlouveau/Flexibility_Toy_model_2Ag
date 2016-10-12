function [ new_exit_cells, b_cells_trial ] = GC_cycle( b_cells_trial, conc, a_act, a_threshold, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection)
%GC_cycle will modify the b_cells matrix (remove those that have lethal
%mutation etc...) and it will create the output of new_plasma_cells

%%DARK ZONE: mutation
%% for each B cell, determine if there is mutation, whether it's in the CDR or FR and the type of mutation. Then change the affinities (or delete) accordingly.
%disp(['GC_cycle line 7 ' num2str(size(b_cells_trial))]);
%disp(['number of b cells beginning of cycle ' num2str(length(b_cells))]);
% 1st division + SHM
daughters1 = division_and_mutation(b_cells_trial, a_act, a_threshold, p_mut, p_CDR, p_FR_lethal);
%disp(['number of b cells after division 1 ' num2str(length(b_cells))]);
%2nd division + SHM
b_cells_trial = division_and_mutation(daughters1, a_act, a_threshold, p_mut, p_CDR, p_FR_lethal);
%disp(['number of b cells after division 2 ' num2str(length(b_cells))]);

%%LIGHT ZONE: selection
%% B cells that have affinity at least higher than a threshold and are in the top portion remain. 

b_cells_trial = selection(b_cells_trial, conc, a_act, t_cell_selection);

%%RECYCLE
%% randomly pick exit_cells for the selected b_cells.
n_selected = floor(t_cell_selection*size(b_cells_trial,3));

n_exit = floor((1 - p_recycle)*n_selected);
%disp(['n exit ' num2str(n_exit)]);

new_exit_cells = zeros(1, size(b_cells_trial,2), n_exit);

for k = 1:n_exit
    ind = randi(size(b_cells_trial,3));
        new_exit_cells(:, :,k) = b_cells_trial(:,:,ind); %add the chosen b cells to the exit b cells
        b_cells_trial(:,:,ind) = []; % remove that b_cells from the list of GC b cells.
end


