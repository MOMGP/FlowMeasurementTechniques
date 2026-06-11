clear;
clc;
close all;

cfg = defaultConfig();
ensureOutputFolders(cfg);

fprintf('Exporting main-branch calculations for Python plotting...\n');
fprintf('Repository root: %s\n', cfg.repoRoot);
fprintf('Output folder:   %s\n\n', cfg.outputDir);

[hwaProfiles, calibrationPoints, calibrationCurve] = processHwaLikeMain(cfg);
writetable(hwaProfiles, fullfile(cfg.tablesDir, 'hwa_profiles.csv'));
writetable(calibrationPoints, fullfile(cfg.tablesDir, 'hwa_calibration_points.csv'));
writetable(calibrationCurve, fullfile(cfg.tablesDir, 'hwa_calibration_curve.csv'));

[pivMeanFields, pivRmsFields, pivProfiles, pivProcessingSummary] = processPivLikeMain(cfg);
writetable(pivMeanFields, fullfile(cfg.tablesDir, 'piv_mean_fields.csv'));
writetable(pivRmsFields, fullfile(cfg.tablesDir, 'piv_rms_fields.csv'));
writetable(pivProfiles, fullfile(cfg.tablesDir, 'piv_profile_xc12.csv'));
writetable(pivProcessingSummary, fullfile(cfg.tablesDir, 'piv_processing_summary.csv'));

[comparisonProfiles, comparisonAtHwaPoints] = compareHwaAndPiv(cfg, hwaProfiles, pivProfiles);
writetable(comparisonProfiles, fullfile(cfg.tablesDir, 'method_comparison_profiles.csv'));
writetable(comparisonAtHwaPoints, fullfile(cfg.tablesDir, 'method_comparison_at_hwa_points.csv'));

fprintf('Done. CSV tables are in:\n%s\n', cfg.tablesDir);
fprintf('Next step: python analysis\\plot_main_outputs.py\n');

function cfg = defaultConfig()
    scriptDir = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(scriptDir);
    if ~isfolder(fullfile(repoRoot, 'data'))
        repoRoot = scriptDir;
    end

    cfg.repoRoot = repoRoot;
    cfg.dataDir = fullfile(repoRoot, 'data');
    cfg.outputDir = fullfile(repoRoot, 'outputs');
    cfg.tablesDir = fullfile(cfg.outputDir, 'tables');
    cfg.figuresDir = fullfile(cfg.outputDir, 'figures');

    cfg.hwaDir = fullfile(cfg.dataDir, 'HWA', 'Group15');
    cfg.pivProcessedDir = fullfile(cfg.dataDir, 'PIV processed');

    cfg.trailingEdgeY_mm = 55;

    % In the DaVis coordinate system the airfoil wake is downstream toward
    % negative x. The HWA traverse is at x/c = 1.2, i.e. about 20 mm
    % downstream of the trailing edge. The trailing edge is near x = -75 mm
    % in the exported PIV fields, so the comparison station is near -95 mm.
    cfg.pivProfileX_mm = -95;

    % DaVis Vx is negative for the tunnel streamwise direction in this data.
    cfg.pivStreamwiseSign = -1;

    cfg.pivCases = struct( ...
        'aoa_deg', {0, 5, 15}, ...
        'case_id', {'aoa0', 'aoa5', 'aoa15'}, ...
        'processedFolder', {'ProcessedAoa0', 'Processed5', 'Processed15'});
end

function ensureOutputFolders(cfg)
    folders = {cfg.outputDir, cfg.tablesDir, cfg.figuresDir};
    for iFolder = 1:numel(folders)
        if ~isfolder(folders{iFolder})
            mkdir(folders{iFolder});
        end
    end
end

function [profiles, calibrationPoints, calibrationCurve] = processHwaLikeMain(cfg)
    calibrationPath = fullfile(cfg.hwaDir, 'calibration.txt');
    calibrationData = readmatrix(calibrationPath);
    calibrationData = calibrationData(:, 1:2);
    calibrationData = calibrationData(all(isfinite(calibrationData), 2), :);

    uCal = calibrationData(:, 1);
    signalCal = calibrationData(:, 2);

    % Same calibration strategy as HW_Post.m on main.
    [signalUnique, ~, groupIdx] = unique(signalCal, 'sorted');
    uUnique = accumarray(groupIdx, uCal, [], @mean);

    calibrationPoints = table(uCal, signalCal, ...
        'VariableNames', {'velocity_ms', 'signal_v'});

    signalGrid = linspace(min(signalUnique), max(signalUnique), 500).';
    calibrationCurve = table(signalGrid, ...
        interp1(signalUnique, uUnique, signalGrid, 'pchip', 'extrap'), ...
        'VariableNames', {'signal_v', 'velocity_ms'});

    dataFiles = dir(fullfile(cfg.hwaDir, 'aoa*_10.txt'));
    records = struct('aoa_deg', {}, 'position_mm', {}, 'y_relative_mm', {}, ...
        'mean_velocity_ms', {}, 'rms_velocity_ms', {}, 'sample_rate_hz', {}, ...
        'n_samples', {}, 'file_name', {});

    for iFile = 1:numel(dataFiles)
        fileName = dataFiles(iFile).name;
        tokens = regexp(fileName, '^aoa(-?\d+)_([0-9]+)_10\.txt$', 'tokens', 'once');
        if isempty(tokens)
            continue;
        end

        aoa = str2double(tokens{1});
        position_mm = str2double(tokens{2});
        filePath = fullfile(dataFiles(iFile).folder, fileName);

        rawData = readmatrix(filePath);
        time = rawData(:, 1);
        probeSignal = rawData(:, 2);
        valid = isfinite(time) & isfinite(probeSignal);
        time = time(valid);
        probeSignal = probeSignal(valid);

        velocity = interp1(signalUnique, uUnique, probeSignal, 'pchip', 'extrap');
        meanVelocity = mean(velocity, 'omitnan');
        rmsVelocity = sqrt(mean((velocity - meanVelocity).^2, 'omitnan'));
        sampleRate_hz = 1 / mean(diff(time));

        records(end + 1).aoa_deg = aoa; %#ok<AGROW>
        records(end).position_mm = position_mm;
        records(end).y_relative_mm = position_mm - cfg.trailingEdgeY_mm;
        records(end).mean_velocity_ms = meanVelocity;
        records(end).rms_velocity_ms = rmsVelocity;
        records(end).sample_rate_hz = sampleRate_hz;
        records(end).n_samples = numel(velocity);
        records(end).file_name = string(fileName);
    end

    if isempty(records)
        error('No valid HWA records could be processed.');
    end

    profiles = sortrows(struct2table(records), {'aoa_deg', 'position_mm'});
end

function [meanFields, rmsFields, profiles, processingSummary] = processPivLikeMain(cfg)
    meanFields = table();
    rmsFields = table();
    profiles = table();
    processingSummary = table();

    for iCase = 1:numel(cfg.pivCases)
        pivCase = cfg.pivCases(iCase);
        caseDir = fullfile(cfg.pivProcessedDir, pivCase.processedFolder);

        [avgFile, stdevFile] = findAverageAndStdevFiles(fullfile(caseDir, 'Overlap50MP3AvgStDev'));
        avgField = readDavisField(avgFile, true);
        stdevField = readDavisField(stdevFile, true);

        meanFields = [meanFields; fieldToTable(avgField, pivCase, cfg, 'Overlap50MP3')]; %#ok<AGROW>
        rmsFields = [rmsFields; rmsFieldToTable(stdevField, pivCase, cfg, 'Overlap50MP3')]; %#ok<AGROW>
        profiles = [profiles; extractPivProfile(cfg, pivCase, avgField, stdevField, 'Overlap50MP3')]; %#ok<AGROW>
        processingSummary = [processingSummary; summarizeProcessingFolders(cfg, pivCase)]; %#ok<AGROW>
    end
end

function [avgFile, stdevFile] = findAverageAndStdevFiles(folderPath)
    avgFile = findFirstExisting(folderPath, {'avg.dat', 'Avg.dat', 'B00001.dat'});
    stdevFile = findFirstExisting(folderPath, {'stdev.dat', 'RMS.dat', 'B00002.dat'});
    if avgFile == "" || stdevFile == ""
        error('Could not find average/stdev files in %s.', folderPath);
    end
end

function path = findFirstExisting(folderPath, names)
    path = "";
    for iName = 1:numel(names)
        candidate = fullfile(folderPath, names{iName});
        if isfile(candidate)
            path = string(candidate);
            return;
        end
    end
end

function field = readDavisField(filePath, trimEdge)
    filePath = char(filePath);
    [nCols, nRows] = extractDavisGridSize(filePath);
    data = importdata(filePath, ' ', 3);
    M = data.data;

    field.X = reshape(M(:, 1), nCols, nRows).';
    field.Y = reshape(M(:, 2), nCols, nRows).';
    field.Vx = reshape(M(:, 3), nCols, nRows).';
    field.Vy = reshape(M(:, 4), nCols, nRows).';
    field.valid = reshape(M(:, 5), nCols, nRows).' > 0;

    if trimEdge
        field.X = field.X(1:end - 1, 1:end - 1);
        field.Y = field.Y(1:end - 1, 1:end - 1);
        field.Vx = field.Vx(1:end - 1, 1:end - 1);
        field.Vy = field.Vy(1:end - 1, 1:end - 1);
        field.valid = field.valid(1:end - 1, 1:end - 1);
    end
end

function [nCols, nRows] = extractDavisGridSize(filePath)
    fid = fopen(filePath, 'r');
    if fid < 0
        error('Cannot open %s.', filePath);
    end
    cleanup = onCleanup(@() fclose(fid));
    fgetl(fid);
    fgetl(fid);
    thirdLine = fgetl(fid);
    tokens = regexp(thirdLine, 'I=(\d+),\s*J=(\d+)', 'tokens', 'once');
    if isempty(tokens)
        error('Could not parse DaVis grid size from %s.', filePath);
    end
    nCols = str2double(tokens{1});
    nRows = str2double(tokens{2});
end

function result = fieldToTable(field, pivCase, cfg, processing)
    u = cfg.pivStreamwiseSign * field.Vx;
    v = field.Vy;
    result = table( ...
        repmat(string(pivCase.case_id), numel(field.X), 1), ...
        repmat(pivCase.aoa_deg, numel(field.X), 1), ...
        repmat(string(processing), numel(field.X), 1), ...
        field.X(:), field.Y(:), field.Y(:) - cfg.trailingEdgeY_mm, ...
        field.Vx(:), field.Vy(:), u(:), v(:), hypot(u(:), v(:)), field.valid(:), ...
        'VariableNames', {'case_id', 'aoa_deg', 'processing', 'x_mm', 'y_mm', ...
        'y_relative_mm', 'davis_vx_ms', 'davis_vy_ms', 'u_streamwise_ms', ...
        'v_normal_ms', 'speed_ms', 'valid'});
end

function result = rmsFieldToTable(field, pivCase, cfg, processing)
    result = table( ...
        repmat(string(pivCase.case_id), numel(field.X), 1), ...
        repmat(pivCase.aoa_deg, numel(field.X), 1), ...
        repmat(string(processing), numel(field.X), 1), ...
        field.X(:), field.Y(:), field.Y(:) - cfg.trailingEdgeY_mm, ...
        abs(field.Vx(:)), abs(field.Vy(:)), hypot(field.Vx(:), field.Vy(:)), field.valid(:), ...
        'VariableNames', {'case_id', 'aoa_deg', 'processing', 'x_mm', 'y_mm', ...
        'y_relative_mm', 'u_rms_ms', 'v_rms_ms', 'speed_rms_ms', 'valid'});
end

function profile = extractPivProfile(cfg, pivCase, avgField, stdevField, processing)
    xColumns = mean(avgField.X, 1, 'omitnan');
    [~, selectedColumn] = min(abs(xColumns - cfg.pivProfileX_mm));

    valid = avgField.valid(:, selectedColumn) & stdevField.valid(:, selectedColumn);
    y = avgField.Y(:, selectedColumn);
    u = cfg.pivStreamwiseSign * avgField.Vx(:, selectedColumn);
    uRms = abs(stdevField.Vx(:, selectedColumn));

    profile = table( ...
        repmat(string(pivCase.case_id), numel(y), 1), ...
        repmat(pivCase.aoa_deg, numel(y), 1), ...
        repmat(string(processing), numel(y), 1), ...
        repmat(xColumns(selectedColumn), numel(y), 1), ...
        y, y - cfg.trailingEdgeY_mm, u, uRms, valid, ...
        'VariableNames', {'case_id', 'aoa_deg', 'processing', 'x_profile_mm', ...
        'y_mm', 'y_relative_mm', 'mean_velocity_ms', 'rms_velocity_ms', 'valid'});
end

function summary = summarizeProcessingFolders(cfg, pivCase)
    caseDir = fullfile(cfg.pivProcessedDir, pivCase.processedFolder);
    folders = dir(fullfile(caseDir, '*AvgStDev'));
    summary = table();

    for iFolder = 1:numel(folders)
        folderPath = fullfile(folders(iFolder).folder, folders(iFolder).name);
        try
            [avgFile, stdevFile] = findAverageAndStdevFiles(folderPath);
            avgField = readDavisField(avgFile, true);
            stdevField = readDavisField(stdevFile, true);
        catch
            continue;
        end

        valid = avgField.valid;
        u = cfg.pivStreamwiseSign * avgField.Vx;
        uRms = abs(stdevField.Vx);

        row = table( ...
            string(pivCase.case_id), pivCase.aoa_deg, ...
            erase(string(folders(iFolder).name), 'AvgStDev'), ...
            nnz(valid), numel(valid), nnz(valid) / numel(valid), ...
            mean(u(valid), 'omitnan'), std(u(valid), 'omitnan'), ...
            mean(uRms(valid), 'omitnan'), ...
            'VariableNames', {'case_id', 'aoa_deg', 'processing', ...
            'n_valid_vectors', 'n_total_vectors', 'valid_fraction', ...
            'mean_u_streamwise_ms', 'spatial_std_u_streamwise_ms', ...
            'mean_u_rms_ms'});
        summary = [summary; row]; %#ok<AGROW>
    end
end

function [comparisonProfiles, comparisonAtHwaPoints] = compareHwaAndPiv(cfg, hwaProfiles, pivProfiles)
    hwa = table( ...
        repmat("HWA", height(hwaProfiles), 1), ...
        string("hwa_" + string(hwaProfiles.aoa_deg)), hwaProfiles.aoa_deg, ...
        hwaProfiles.y_relative_mm, hwaProfiles.mean_velocity_ms, ...
        hwaProfiles.rms_velocity_ms, true(height(hwaProfiles), 1), ...
        'VariableNames', {'method', 'case_id', 'aoa_deg', 'y_relative_mm', ...
        'mean_velocity_ms', 'rms_velocity_ms', 'valid'});

    piv = table( ...
        repmat("PIV", height(pivProfiles), 1), ...
        pivProfiles.case_id, pivProfiles.aoa_deg, pivProfiles.y_relative_mm, ...
        pivProfiles.mean_velocity_ms, pivProfiles.rms_velocity_ms, pivProfiles.valid, ...
        'VariableNames', hwa.Properties.VariableNames);

    comparisonProfiles = [hwa; piv];

    rows = table();
    aoaValues = intersect(unique(hwaProfiles.aoa_deg), unique(pivProfiles.aoa_deg));
    for iAoa = 1:numel(aoaValues)
        aoa = aoaValues(iAoa);
        h = hwaProfiles(hwaProfiles.aoa_deg == aoa, :);
        p = pivProfiles(pivProfiles.aoa_deg == aoa & pivProfiles.valid, :);
        if height(p) < 2
            continue;
        end

        [pY, order] = sort(p.y_relative_mm);
        pMean = p.mean_velocity_ms(order);
        pRms = p.rms_velocity_ms(order);
        inRange = h.y_relative_mm >= min(pY) & h.y_relative_mm <= max(pY);
        h = h(inRange, :);

        pivMeanAtHwa = interp1(pY, pMean, h.y_relative_mm, 'linear');
        pivRmsAtHwa = interp1(pY, pRms, h.y_relative_mm, 'linear');

        thisRows = table( ...
            repmat(aoa, height(h), 1), h.y_relative_mm, ...
            h.mean_velocity_ms, pivMeanAtHwa, pivMeanAtHwa - h.mean_velocity_ms, ...
            h.rms_velocity_ms, pivRmsAtHwa, pivRmsAtHwa - h.rms_velocity_ms, ...
            'VariableNames', {'aoa_deg', 'y_relative_mm', ...
            'hwa_mean_velocity_ms', 'piv_mean_velocity_ms', 'piv_minus_hwa_mean_ms', ...
            'hwa_rms_velocity_ms', 'piv_rms_velocity_ms', 'piv_minus_hwa_rms_ms'});
        rows = [rows; thisRows]; %#ok<AGROW>
    end

    comparisonAtHwaPoints = rows;
end
