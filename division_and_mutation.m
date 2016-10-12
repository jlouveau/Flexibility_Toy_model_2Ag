function [ daughters ] = division_and_mutation( b_cells_trial, a_act, a_threshold, p_mut, p_CDR, p_FR_lethal )
%disp(['div and mut line 2 ' num2str(size(b_cells_trial))]);
%division
daughters = cat(3,b_cells_trial, b_cells_trial);
%round of SHM
%disp(['div and mut line 6 ' num2str(size(daughters))]);

for n = 1:size(daughters,3)
    rand_mut = rand;   
    if rand_mut < p_mut %if there is mutation, then call function mutation and modify the daughter cell
        %disp (['number of daughters before mutation ' num2str(size(daughters,3))]);
        mutant = mutation(daughters(1,:,n), a_act, a_threshold, p_CDR, p_FR_lethal);
        if ~isempty(mutant)
            daughters(1,:,n) = mutant(:);
        else            
            daughters(1,:,n) = [NaN NaN];
        end
    end
end

%disp(['div and mut line 20 ' num2str(size(daughters))]);
daughters = daughters(:, :, isnan(daughters(1,1,:)) < 1);
%disp(['div and mut line 22 ' num2str(size(daughters))]);
%disp (['number of daughters after mutation ' num2str(size(daughters,2))]);

