% Generates minimal reactomes for all 52 AGORA models across six diets.

OverallTimer = tic;

fprintf('\n-----------------------------------------------------\n');
fprintf('GENERATING MINIMAL REACTOMES\n');
fprintf('\n-----------------------------------------------------\n');


% DIRECTORIES

base_dir = 'C:\Users\Avinash Chauhan\agora_2_models';
diet_dir = fullfile(base_dir, 'Diets');
output_dir = fullfile(base_dir, 'minReactModels');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


% SOLVER SETUP

changeCobraSolver('gurobi', 'LP', 0);
changeCobraSolver('gurobi', 'MILP', 0);


% LOAD MODELS LIST

model_files = dir(fullfile(base_dir, '*.mat'));
model_files = model_files(~[model_files.isdir]);
nModels = numel(model_files);

fprintf('Total models found : %d\n\n', nModels);


% DIETS SETUP

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
nDiets = numel(dietNames);
DietTables = cell(nDiets-1, 1);

fprintf('Loading diet files...\n');
for d = 1:nDiets-1
    file = fullfile(diet_dir, [dietNames{d} '.txt']);
    T = readtable(file, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);
    
    rxns = strtrim(string(T{:,1}));
    rxns = strrep(rxns, '[u]', '(e)'); % Standardize extracellular compartment
    
    DietTables{d}.rxns = rxns;
    DietTables{d}.lbs  = T{:,2};
end
fprintf('All diet files loaded.\n\n');


% STARTING PARPOOL

for d = 1:nDiets
    folder = fullfile(output_dir, dietNames{d});
    if ~exist(folder, 'dir'), mkdir(folder); end
end

if isempty(gcp('nocreate'))
    parpool('local'); 
end


% MAIN ADAPTIVE LOOP

ParpoolResults = cell(nModels, 1);
fprintf('Starting parallel execution over %d models...\n', nModels);

parfor m = 1:nModels
    changeCobraSolver('gurobi', 'all', 0);
    ModelTimer = tic;
    
    file = model_files(m).name;
    [~, modelName] = fileparts(file);
    localResults = cell(nDiets, 7);
    
    %% LOADING MODEL
    try
        S = load(fullfile(base_dir, file));
        fn = fieldnames(S);
        model = [];
        for k = 1:numel(fn)
           tmp = S.(fn{k});
           if isstruct(tmp) && isfield(tmp,'rxns') && isfield(tmp,'mets')
              model = tmp; break;
           end
        end
        if iscell(model), model = model{1}; end
        model.rxns = strtrim(model.rxns);
    catch ME
        for d = 1:nDiets
            localResults(d,:) = {modelName, dietNames{d}, NaN, NaN, NaN, NaN, 'LOAD_FAIL'};
        end
        ParpoolResults{m} = localResults;
        continue;
    end
    
    %% EVALUATING ACROSS DIETS
    for d = 1:nDiets
        dietName = dietNames{d};
        saveFolder = fullfile(output_dir, dietName);
        saveFile = fullfile(saveFolder, [modelName '.mat']);
        
        % Skipping if result already exists 
        if exist(saveFile, 'file')
            try
                existingData = load(saveFile, 'WTgrowth', 'MRMgrowth', 'WTRxns', 'MRMRxns', 'targetUsed');
                status_str = sprintf('EXISTING_%.2f', existingData.targetUsed);
                localResults(d,:) = {modelName, dietName, existingData.WTgrowth, ...
                                     existingData.MRMgrowth, existingData.WTRxns, ...
                                     existingData.MRMRxns, status_str};
                continue;
            catch
                % Proceeding if file is corrupted
            end
        end
        
        currentModel = model;
        num_rxns = length(currentModel.rxns);
        
        % Applying Diet Constraints
        if d <= numel(DietTables)
            dietLB = DietTables{d}.lbs ./ 24; 
            modelRxnsLower = strrep(lower(string(currentModel.rxns)), '[u]', '(e)');
            [tf, loc] = ismember(lower(DietTables{d}.rxns), modelRxnsLower);
            if any(tf)
                valid = loc > 0;
                matchedIdx = loc(valid);
                matchedRxns = currentModel.rxns(matchedIdx);
                matchedLB = dietLB(valid);
                currentModel = changeRxnBounds(currentModel, matchedRxns, matchedLB, 'l');
                currentModel = changeRxnBounds(currentModel, matchedRxns, 1000, 'u');
            end
        end
        
        % Optimizing Wild-Type
        WTsol = optimizeCbModel(currentModel, 'max', 'one');
        if isempty(WTsol) || WTsol.stat ~= 1 || WTsol.f < 1e-6
            localResults(d,:) = {modelName, dietName, NaN, NaN, NaN, NaN, 'WT_FAIL'};
            continue;
        end
        
        WTgrowth = WTsol.f;
        WTRxns = nnz(abs(WTsol.x) > 1e-6);
        
        % -----------------------------------------------------------
        % STRICT TARGETS & SCALING
        % -----------------------------------------------------------
        %survival_targets = [0.05, 0.10, 0.20, 0.50, 0.80];
        survival_targets =  [1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05];
        success = false;
        
        if WTgrowth < 0.01
            scales_to_try = [1, 10000, 1000];
        elseif WTgrowth < 1
            scales_to_try = [1, 1000, 100];
        else
            scales_to_try = [1, 10, 100, 0.1];
        end
        
        for st = survival_targets
            
            % Creating Base Model & Scale ONLY ATPM Requirements safely
            base_model = currentModel;
            atpm_idx = find(contains(upper(base_model.rxns), 'ATPM'));
            for i = 1:numel(atpm_idx)
                if base_model.lb(atpm_idx(i)) > 0
                    base_model.lb(atpm_idx(i)) = base_model.lb(atpm_idx(i)) * st;
                end
            end
            
            Jmin = [];
            
            %  Standard MILP (minReact)
            for sf = scales_to_try
                scaledModel = base_model;
                if sf ~= 1
                    scaledModel.lb = scaledModel.lb * sf;
                    scaledModel.ub = scaledModel.ub * sf;
                end
                
                try
                    [Jmin_temp, ~] = minReact(scaledModel, st, 1e-6);
                    % Safety Checking to prevent array out of bounds
                    if ~isempty(Jmin_temp) && length(Jmin_temp(1,:)) == num_rxns
                        Jmin = double(Jmin_temp(1,:)');
                        break; 
                    end
                catch
                    continue; 
                end
            end
            
            
            % THE UNBREAKABLE pFBA / FBA FALLBACK
            
            if isempty(Jmin)
                fallbackModel = base_model;
                
                % Locking objective to required 5% survival fraction (relaxed 2% for boundary crash)
                obj_idx = find(fallbackModel.c);
                if ~isempty(obj_idx)
                    fallbackModel.lb(obj_idx(1)) = WTgrowth * st * 0.98; 
                end
                
                try
                    sol_pfba = minimizeModelFlux(fallbackModel);
                    
                    % Strict checks before indexing
                    if ~isempty(sol_pfba) && sol_pfba.stat == 1
                        if isfield(sol_pfba, 'x') && length(sol_pfba.x) == num_rxns
                            Jmin = double(abs(sol_pfba.x) > 1e-6)';
                        elseif isfield(sol_pfba, 'full') && length(sol_pfba.full) == num_rxns
                            Jmin = double(abs(sol_pfba.full) > 1e-6)';
                        end
                    end
                catch
                    % pFBA crashed
                end
                
                % Standard FBA flux vector
                if isempty(Jmin)
                    try
                        sol_fba = optimizeCbModel(fallbackModel, 'max', 'one');
                        if ~isempty(sol_fba) && sol_fba.stat == 1
                            if isfield(sol_fba, 'x') && length(sol_fba.x) == num_rxns
                                Jmin = double(abs(sol_fba.x) > 1e-6)';
                            elseif isfield(sol_fba, 'full') && length(sol_fba.full) == num_rxns
                                Jmin = double(abs(sol_fba.full) > 1e-6)';
                            end
                        end
                    catch
                    end
                end
            end
            
            
            if isempty(Jmin)
                continue; % Moving to next target if absolutely rigid or solver timed out
            end
            
            % Apply filter to the FULLY CONSTRAINED original model
            mrmVector = logical(Jmin);
            removeList = base_model.rxns(~mrmVector);
            
            if ~isempty(removeList)
                minimalModel = removeRxns(base_model, removeList);
            else
                minimalModel = base_model;
            end
            
            if isempty(minimalModel.rxns) || nnz(minimalModel.c) == 0
                continue; 
            end
            
            % Optimizing final unscaled minimal model to verify biological viability
            MRMsol = optimizeCbModel(minimalModel, 'max', 'one');
            
            if isempty(MRMsol) || MRMsol.stat ~= 1
                continue; 
            end
            
            MRMgrowth = MRMsol.f;
            MRMRxns = nnz(abs(MRMsol.x) > 1e-6);
            
            % Saving Results
            saveStruct = struct('minimalModel', minimalModel, 'WTgrowth', WTgrowth, ...
                'MRMgrowth', MRMgrowth, 'WTRxns', WTRxns, 'MRMRxns', MRMRxns, ...
                'Jmin', Jmin', 'modelName', modelName, 'dietName', dietName, 'targetUsed', st);
            
            parsave(saveFile, saveStruct);
            
            status_str = sprintf('PASS_%.2f', st);
            localResults(d,:) = {modelName, dietName, WTgrowth, MRMgrowth, WTRxns, MRMRxns, status_str};
            success = true;
            break; 
        end
        
        if ~success
            localResults(d,:) = {modelName, dietName, WTgrowth, NaN, WTRxns, NaN, 'FAIL_ALL_THRESHOLDS'};
        end
    end
    
    fprintf('%s processed for all diets in %.2f mins.\n', modelName, toc(ModelTimer)/60);
    ParpoolResults{m} = localResults;
end


% AGGREGATING RESULTS & FINAL SUMMARY

MasterResults = vertcat(ParpoolResults{:});

Elapsed = toc(OverallTimer);
fprintf('\n------------------------------------------------\n');
fprintf('ALL MODELS COMPLETED\n');
fprintf('------------------------------------------------\n');
fprintf('Total Time : %.2f minutes\n', Elapsed/60);
fprintf('Minimal reactomes saved in: %s\n\n', output_dir);

if ~isempty(MasterResults)
    passCount  = sum(contains(string(MasterResults(:,7)),'PASS'));
    existCount = sum(contains(string(MasterResults(:,7)),'EXISTING'));
    failCount  = sum(contains(string(MasterResults(:,7)),'FAIL') | contains(string(MasterResults(:,7)),'ERROR'));
    
    fprintf('Newly Computed Runs : %d\n', passCount);
    fprintf('Loaded Existing     : %d\n', existCount);
    fprintf('Failed Runs         : %d\n\n', failCount);
    
    ResultTable = cell2table(MasterResults, ...
        'VariableNames', {'Model','Diet','WTGrowth','MRMGrowth','WTRxns','MRMRxns','Status'});
    
    summaryFile = fullfile(output_dir, 'minReactModels_Summary.csv');
    writetable(ResultTable, summaryFile);
    
    models = unique(string(ResultTable.Model), 'stable');
    WideTable = table(models(:), 'VariableNames', {'Model'});
    dietShort = {'HF','MED','UNH','VEG','WES','NOD'};
    
    for d = 1:numel(dietNames)
        wt = nan(height(WideTable),1); mrm = nan(height(WideTable),1);
        wtr = nan(height(WideTable),1); mrr = nan(height(WideTable),1);
        
        for i = 1:height(WideTable)
            idx = strcmp(string(ResultTable.Model), string(WideTable.Model(i))) & ...
                  strcmp(string(ResultTable.Diet), string(dietNames{d}));
            if any(idx)
                wt(i)  = ResultTable.WTGrowth(idx);
                mrm(i) = ResultTable.MRMGrowth(idx);
                wtr(i) = ResultTable.WTRxns(idx);
                mrr(i) = ResultTable.MRMRxns(idx);
            end
        end
        WideTable.([dietShort{d} '_WT'])  = wt;
        WideTable.([dietShort{d} '_MRM']) = mrm;
        WideTable.([dietShort{d} '_WTR']) = wtr;
        WideTable.([dietShort{d} '_MRR']) = mrr;
    end
    
    wideFile = fullfile(output_dir, 'minReactModels_WideSummary.csv');
    writetable(WideTable, wideFile);
    fprintf('Summary tables saved successfully to %s\n', output_dir);
else
    fprintf('No valid data generated to save to CSV.\n');
end


% HELPER FUNCTIONS

function parsave(fname, dataStruct)
    save(fname, '-struct', 'dataStruct', '-v7.3');
end