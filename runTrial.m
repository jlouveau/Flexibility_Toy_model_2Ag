function [ b_cells_trial, number_recycled_b_cells_trial, exit_cells_trial, number_exit_cells_trial, final_cycle ] = runTrial( b_cells_trial, exit_cells_trial, number_recycled_b_cells_trial, number_exit_cells_trial, conc, a_act, a_threshold, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection, n_max_Bcells, n_cycle_max, cycle_number )

cycle_number = 2;

%%GC REACTION:
%% the founders seed the GC and undergo AM. GC_cycle modifies b_cells and adds the new plasma cells. In this toy model, exit_cells designed the cells that exit the GC at the end of a cycle (memory + plasma cells).

%disp('RunTrial l12');
while 1
    cycle_number = cycle_number +1;
    %disp(['Cycle Number ' num2str(cycle_number)]);

    [new_exit_cells, b_cells_trial] = GC_cycle(b_cells_trial, conc, a_act, a_threshold, p_mut, p_CDR, p_FR_lethal, p_recycle, t_cell_selection);
    number_recycled_b_cells_trial(cycle_number) = size(b_cells_trial,3);
    number_exit_cells_trial(cycle_number) = size(new_exit_cells, 3);
    
    for l = 1:number_exit_cells_trial(cycle_number)
        for m = 1:size(new_exit_cells,2)
            exit_cells_trial(1, cycle_number, m, l) = new_exit_cells(1,m,l);
        end
    end
    
    if cycle_number > n_cycle_max -1 || number_recycled_b_cells_trial(cycle_number) > n_max_Bcells || isempty(b_cells_trial)
        break
    end
end

final_cycle = cycle_number - 1;
end
