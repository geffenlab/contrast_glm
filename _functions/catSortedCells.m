function [spikeData, sessionData] = catSortedCells(pathList,ops)

% load template
t = load(ops.template_dir);

sess_cnt = 0;
unit_cnt = 0;
sessionData = [];
spikeData = [];
for i = 1:size(pathList,2)
    for j = 1:size(pathList{i},1)

        % get the root folder
        parts = split(pathList{i}{j},'/');
        root = ['/' fullfile(parts{1:9})];
        
        fprintf('_________________________________________\n')
        fprintf('MOUSE %s - SESSION %s\n',parts{6},parts{7})
        
        % separate blocks
        [block, cellInfo, labels,~,events] = ...
            splitSpikesByTemplates_allev(root,...
                                         t.template, ...
                                         ops.drift_correct);
        
        
        % session info
        rm_fields = {'spikes','clust','clustID','clustLabel'};
        sess_cnt = sess_cnt + 1;
        sessionData(sess_cnt).sessID = sprintf('%s_%s',...
                                               cellInfo{1,1}, ...
                                               cellInfo{1,2});
        sessionData(sess_cnt).block = rmfield(block,rm_fields);
        sessionData(sess_cnt).events = events;
        
        % load waveforms for this session
        wave = extractWaveforms(pathList{i}{j});
        
        % spike info
        for c = 1:size(cellInfo,1)
            
            unit_cnt = unit_cnt+1;
            
            % cell ID
            if strcmp(cellInfo{c,5},'good')
                lbl = 'su';
            else
                lbl = 'mu';
            end
            spikeData(unit_cnt).cellID = sprintf('%s_%03d_%03d_%s',...
                                                 sessionData(sess_cnt).sessID,...
                                                 cellInfo{c,4}, ...
                                                 cellInfo{c,3},...
                                                 lbl);
            spikeData(unit_cnt).sessionID = sessionData(sess_cnt).sessID;
            
            % spikes per block
            for b = 1:length(block)
                spikeData(unit_cnt).spikes{b} = ...
                    block(b).spikes(block(b).clust == block(b).clustID(c));
            end
            
            % waveform info
            spikeData(unit_cnt).spkcount = wave.spkcount(c);
            spikeData(unit_cnt).mfr = cellInfo{c,7};
            spikeData(unit_cnt).peakchan = wave.peakchan(c);
            spikeData(unit_cnt).waveform = wave.peakwaveform(c,:);
            spikeData(unit_cnt).trough_peak = wave.trough_peak(c);
            spikeData(unit_cnt).peak_inflect = wave.peak_inflect(c);
            spikeData(unit_cnt).FWHM = wave.FWHM(c);
            
        end
        
    end
    
end