function create_bads_reports(base_folder_path, subs, params)
% Import necessary packages
import mlreportgen.report.*
import mlreportgen.dom.*

% Loop over each subject in the list
for subNumber = subs
    % Convert subNumber to a two-digit string
    subStr = sprintf('%02d', subNumber);

    % Create the report name
    reportName = ['sub_', subStr, '_sensorlevel'];

    % Define the subject folder path
    subjectFolderPath = fullfile(base_folder_path, ['sub_', subStr]);

    % Create a new report
    rpt = Report(fullfile(subjectFolderPath, reportName), 'pdf');

    % Add a title page
    titlepg = TitlePage;
    titlepg.Title = ['Subject ', subStr, ' Preprocessing'];
    add(rpt, titlepg);

    % Add the table of contents
    toc = TableOfContents();
    add(rpt, toc);

    % Define sections
    sections = {'opm', 'squid', 'opmeeg', 'squideeg'};
    sections2 = {'opm', 'squidmag', 'opmeeg', 'squideeg'};

    %%
    chapter = Chapter('Parameters');
    processStruct(params, 'params', chapter);
    add(rpt, chapter);
    
    
    %% Bad Channels chapter
    chapter = Chapter('Bad Channels');
    chapter.Numbered = false; % Remove chapter numbering
    for i_section = 1:length(sections)
        filePath = fullfile(subjectFolderPath, ['sub_', subStr, '_' sections{i_section} '_badchs.mat']);
        if isfile(filePath)
            data = load(filePath);
            section = Section(sections(i_section));
            section.Numbered = false; % Remove section numbering
            fields = fieldnames(data);
            totalLength = 0;
            for j = 1:length(fields)
                varData = data.(fields{j});
                if isempty(varData)
                    varLength = 0;
                    varList = '';
                elseif iscell(varData)
                    varLength = length(varData);
                    varList = strjoin(varData, ', ');
                else
                    varLength = max(size(varData));
                    varList = num2str(varData(:)');
                    varList = regexprep(varList, '\s+', ', ');
                end
                totalLength = totalLength + varLength;
                para = Paragraph([fields{j},' (n=', num2str(varLength), '): ', varList]);
                add(section, para);
            end
            para = Paragraph(['n_{total}: ', num2str(totalLength)]);
            add(section, para);
            add(chapter, section);
        end
    end
    add(rpt, chapter);

    %% Bad Trials chapter
    chapter = Chapter('Bad Trials');
    chapter.Numbered = false; % Remove chapter numbering
    for i_section = 1:length(sections)
        filePath = fullfile(subjectFolderPath, ['sub_', subStr, '_' sections{i_section} '_badtrls.mat']);
        if isfile(filePath)
            data = load(filePath);
            section = Section(sections(i_section));
            section.Numbered = false; % Remove section numbering
            fields = fieldnames(data);
            totalLength = 0;
            for j = 1:length(fields)
                varData = data.(fields{j});
                if isempty(varData)
                    varLength = 0;
                    varList = '';
                elseif iscell(varData)
                    varLength = length(varData);
                    varList = strjoin(varData, ', ');
                else
                    varLength = max(size(varData));
                    varList = num2str(varData(:)');
                    varList = regexprep(varList, '\s+', ', ');
                end
                totalLength = totalLength + varLength;
                para = Paragraph([fields{j},' (n=', num2str(varLength), '): ', varList]);
                add(section, para);
            end
            para = Paragraph(['n_{total}: ', num2str(totalLength)]);
            add(section, para);
            add(chapter, section);
        end
    end
    add(rpt, chapter);

    %% ICA chapter
    chapter = Chapter('ICA');
    chapter.Numbered = false; % Remove chapter numbering
    for i_section = 1:length(sections)
        filePath = fullfile(subjectFolderPath,['sub_' subStr '_' sections{i_section} '_ica_comp.mat']);
        if isfile(filePath)
            data = load(filePath);
            section = Section(sections(i_section));
            section.Numbered = false; % Remove section numbering
            fields = {'ecg_comp_idx', 'eog1_comp_idx', 'eog2_comp_idx'};
            totalLength = 0;
            for j = 1:length(fields)
                varData = data.(fields{j});
                if isempty(varData)
                    varLength = 0;
                    varList = '';
                elseif iscell(varData)
                    varLength = length(varData);
                    varList = strjoin(varData, ', ');
                else
                    varLength = max(size(varData));
                    varList = num2str(varData(:)');
                    varList = strjoin(arrayfun(@num2str, varData(:)', 'UniformOutput', false), ', ');
                end
                totalLength = totalLength + varLength;
                para = Paragraph([fields{j}, ': ', varList]);
                add(section, para);
            end

            % Define the folder and the pattern
            pattern = ['sub_' subStr '_' sections{i_section} '_ica_rejected_comps*.jpg']; % Replace 'your_string' with the starting string
            files = dir(fullfile(subjectFolderPath, 'figs', pattern));

            if ~isempty(files)
                for i_file = 1:length(files)
                    img = Image(fullfile(subjectFolderPath,'figs',files(i_file).name));
                    img.Style = {ScaleToFit};
                    add(section, img);
                end
            end
            add(chapter, section);
        end
    end
    add(rpt, chapter);

    %% Butterfly plots chapter
    chapter = Chapter('Timelocked');
    chapter.Numbered = false; % Remove chapter numbering
       
    for i_trigger = 1:length(params.trigger_labels)
        % Add rows and cells to the table and insert the images
        section = Section(['Trigger: ' params.trigger_labels(i_trigger)]);
        section.Numbered = false; % Remove section numbering
        
        tbl = Table();
        tbl.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
        for i = 1:2
            row = TableRow();
            for j = 1:2
                imgIndex = (j-1)*2 + i;
                img = Image(fullfile(subjectFolderPath,'figs',['sub_' subStr '_' sections2{imgIndex} '_butterfly_audodd-' params.trigger_labels{i_trigger} '.jpg']));
                img.Style = {Width('8cm'), ScaleToFit};
                entry = TableEntry();
                append(entry, img);
                append(row, entry);
            end
            append(tbl, row);
        end
    
        add(section, tbl);
        add(chapter, section);
        add(chapter, PageBreak());
    end
    add(rpt, chapter);

%     %% M100 chapter
%     chapter = Chapter('M100');
%     chapter.Numbered = false; % Remove chapter numbering
%     for i_section = 1:length(sections)
%         section = Section(sections(i_section));
%         section.Numbered = false; % Remove section numbering
% 
%         if contains(sections(i_section),'eeg')
%             fieldMultiplier = 1e9;
%             fieldUnit = '[nV]';
%         else
%             fieldMultiplier = 1e15;
%             fieldUnit = '[fT]';
%         end

%         % Load the .mat file
%         data = load(fullfile(subjectFolderPath,['sub_' subStr '_' sections2{i_section} '_M100.mat']));
%         M100 = data.M100;
%         
%         % Create the table with the required data
%         num_trigger = length(params.trigger_labels);
%         T = table('Size', [8, num_trigger], 'VariableTypes', repmat({'double'}, 1, num_trigger), 'VariableNames', params.trigger_labels);
%         T.Properties.RowNames = {['peak_amplitude ' fieldUnit], ['max_amplitude ' fieldUnit], ['min_amplitude ' fieldUnit], 'peak_latency [ms]', 'SNR_prestim', 'SNR_stderr', ['std_prestim ' fieldUnit], ['stderr ' fieldUnit]};
%         
%         for i_trigger = 1:num_trigger
%             trigger_data = M100{i_trigger};
%             T{['peak_amplitude ' fieldUnit], params.trigger_labels{i_trigger}} = fieldMultiplier*trigger_data.peak_amplitude;
%             T{['max_amplitude ' fieldUnit], params.trigger_labels{i_trigger}} = fieldMultiplier*trigger_data.max_amplitude;
%             T{['min_amplitude ' fieldUnit], params.trigger_labels{i_trigger}} = fieldMultiplier*trigger_data.min_amplitude;
%             T{'peak_latency [ms]', params.trigger_labels{i_trigger}} = 1e3*trigger_data.peak_latency;
%             T{'SNR_prestim', params.trigger_labels{i_trigger}} = trigger_data.peak_amplitude / trigger_data.prestim_std;
%             T{'SNR_stderr', params.trigger_labels{i_trigger}} = trigger_data.peak_amplitude / trigger_data.std_error;
%             T{['std_prestim ' fieldUnit], params.trigger_labels{i_trigger}} = fieldMultiplier*trigger_data.prestim_std;
%             T{['stderr ' fieldUnit], params.trigger_labels{i_trigger}} = fieldMultiplier*trigger_data.std_error;
%         end
%         
%         % Convert the MATLAB table to a DOM table
%         domTable = Table();
%         domTable.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
%         
%         % Add header row
%         headerRow = TableRow();
%         append(headerRow, TableEntry(' '));
%         for i_trigger = 1:num_trigger
%             headerEntry = TableEntry(params.trigger_labels{i_trigger});
%             headerEntry.Style = {HAlign('center'), Bold()};
%             append(headerRow, headerEntry);
%         end
%         append(domTable, headerRow);
%         
%         formatString = {'%.1f','%.1f','%.1f','%.f','%.1f','%.2f','%.1f','%.1f'};
% 
%         % Add data rows
%         for i_row = 1:height(T)
%             row = TableRow();
%             rowNameEntry = TableEntry(T.Properties.RowNames{i_row});
%             rowNameEntry.Style = {Bold()};
%             append(row, rowNameEntry);
%             for j = 1:num_trigger
%                 entry = TableEntry(num2str(T{i_row, j},formatString{i_row}));
%                 entry.Style = {HAlign('right')}; 
%                 append(row, entry);
%             end
%             append(domTable, row);
%         end
%         add(section,domTable);
%         add(chapter,section)
%     end
% 
%     % Max channel plots
%     for i_trigger = 1:length(params.trigger_labels)
%         % Add rows and cells to the table and insert the images
%         section = Section(['Trigger: ' params.trigger_labels(i_trigger)]);
%         section.Numbered = false; % Remove section numbering
%         
%         tbl = Table();
%         tbl.Style = {Border('solid'), Width('100%'), RowSep('solid'), ColSep('solid')};
%         for i = 1:2
%             row = TableRow();
%             for j = 1:2
%                 imgIndex = (j-1)*2 + i;
%                 img = Image(fullfile(subjectFolderPath,'figs',['sub_' subStr '_' sections2{imgIndex} '_evoked_peakchannel' num2str(i_trigger) '.jpg']));
%                 img.Style = {Width('8cm'), ScaleToFit};
%                 entry = TableEntry();
%                 append(entry, img);
%                 append(row, entry);
%             end
%             append(tbl, row);
%         end
%     
%         add(section, tbl);
%         add(chapter, section);
%         add(chapter, PageBreak());
%     end
% 
%     add(rpt, chapter);

    %% Close the report
    close(rpt);
end

function processStruct(s, parentName, chapter)
    localFields = fieldnames(s); % Make fields local
    for localI = 1:numel(localFields) % Make i local
        fieldName = localFields{localI};
        fieldValue = s.(fieldName);
        fullName = strcat(parentName, '.', fieldName);
        
        if isstruct(fieldValue)
            % If the field is a struct, process it recursively
            processStruct(fieldValue, fullName, chapter);
        elseif ismatrix(fieldValue) && isnumeric(fieldValue)
            % If the field is a 1-D array, convert it to a comma-separated string
            rowStrs = arrayfun(@(i) strjoin(arrayfun(@num2str, fieldValue(i,:), 'UniformOutput', false), ', '), 1:size(fieldValue,1), 'UniformOutput',false);
            valueStr = strjoin(rowStrs, '; ');
            append(chapter, Paragraph([fullName, ': ', valueStr]));
        elseif iscell(fieldValue) && all(cellfun(@ischar, fieldValue))
            % If the field is a cell array of strings, convert it to a comma-separated string
            valueStr = strjoin(fieldValue, ', ');
            append(chapter, Paragraph([fullName, ': ', valueStr]));
        elseif isempty(fieldValue)
            % Otherwise, convert the value to a string and append it
            append(chapter, Paragraph([fullName, ': []']));
        else
            % Otherwise, convert the value to a string and append it
            if size(fieldValue,1) > size(fieldValue,2)
                fieldValue = fieldValue';
            end
            valueStr = num2str(fieldValue);
            append(chapter, Paragraph([fullName, ': ', valueStr]));
        end
    end
end
end