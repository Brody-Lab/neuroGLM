function cell_info = make_cell_info(Cells)
    % makes a reduced Cells structure with no spike times but just info
    % about each cell which is useful to have alongside, for instance, GLM
    % fits
    fields = {'waveform','rec','ap_meta','nTrials','penetration','last_modified','rat','sess_date','sessid','bank','electrode','unitCount','unitISIRatio','unitLRatio','unitIsoDist','unitVppRaw','clusterNotes','hemisphere','probe_serial','distance_from_tip','DV','AP','ML','regions'};
    for f=1:length(fields)
        if isfield(Cells,fields{f})
            cell_info.(fields{f})=Cells.(fields{f});
            switch fields{f}
                case 'waveform'
                    if isfield(cell_info.waveform,'meanWfGlobalRaw')
                        cell_info.waveform=rmfield(cell_info.waveform,'meanWfGlobalRaw');
                    end
                    if isfield(cell_info.waveform,'mean_uV')
                        cell_info.waveform=rmfield(cell_info.waveform,'mean_uV');
                    end                
                case 'rec'
                    cell_info.rec = rmfield(cell_info.rec,'sync');
            end
        end
    end
end

