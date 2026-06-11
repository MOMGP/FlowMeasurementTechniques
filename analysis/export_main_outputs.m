clear;
clc;
close all;

cfg = defaultConfig();
ensureOutputFolders(cfg);

fprintf('Exporting main-branch calculations for Python plotting...\n');
fprintf('Repository root: %s\n', cfg.repoRoot);
fprintf('Output folder:   %s\n\n', cfg.outputDir);

[hwaProfiles, calibrationPoints, calibrationCurve, calibration] = processHwaLikeMain(cfg);
writetable(hwaProfiles, fullfile(cfg.tablesDir, 'hwa_profiles.csv'));
writetable(calibrationPoints, fullfile(cfg.tablesDir, 'hwa_calibration_points.csv'));
writetable(calibrationCurve, fullfile(cfg.tablesDir, 'hwa_calibration_curve.csv'));

[hwaAutocorrelation, hwaConvergence] = processHwaAutocorrelation(cfg, calibration);
writetable(hwaAutocorrelation, fullfile(cfg.tablesDir, 'hwa_autocorrelation.csv'));
writetable(hwaConvergence, fullfile(cfg.tablesDir, 'hwa_convergence.csv'));

[hwaSpectra, hwaSpectralPeaks] = processHwaSpectra(cfg, calibration, hwaProfiles);
writetable(hwaSpectra, fullfile(cfg.tablesDir, 'hwa_spectra.csv'));
writetable(hwaSpectralPeaks, fullfile(cfg.tablesDir, 'hwa_spectral_peaks.csv'));

[pivMeanFields, pivRmsFields, pivInstantaneousFields, pivProfiles, pivProcessingSummary] = processPivLikeMain(cfg);
writetable(pivMeanFields, fullfile(cfg.tablesDir, 'piv_mean_fields.csv'));
writetable(pivRmsFields, fullfile(cfg.tablesDir, 'piv_rms_fields.csv'));
writetable(pivInstantaneousFields, fullfile(cfg.tablesDir, 'piv_instantaneous_fields.csv'));
writetable(pivProfiles, fullfile(cfg.tablesDir, 'piv_profile_xc12.csv'));
writetable(pivProcessingSummary, fullfile(cfg.tablesDir, 'piv_processing_summary.csv'));

[pivProcessingFields, pivWindowFields, pivWindowSummary, pivEnsembleFields, ...
    pivDeltaTFields, selfVsDavisFields, selfVsDavisSummary] = processPivParameterExports(cfg);
writetable(pivProcessingFields, fullfile(cfg.tablesDir, 'piv_processing_fields.csv'));
writetable(pivWindowFields, fullfile(cfg.tablesDir, 'piv_window_size_fields.csv'));
writetable(pivWindowSummary, fullfile(cfg.tablesDir, 'piv_window_size_summary.csv'));
writetable(pivEnsembleFields, fullfile(cfg.tablesDir, 'piv_ensemble_fields.csv'));
writetable(pivDeltaTFields, fullfile(cfg.tablesDir, 'piv_delta_t_fields.csv'));
writetable(selfVsDavisFields, fullfile(cfg.tablesDir, 'piv_self_vs_davis_fields.csv'));
writetable(selfVsDavisSummary, fullfile(cfg.tablesDir, 'piv_self_vs_davis_summary.csv'));

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

    cfg.hwaTrailingEdgeY_mm = 55;
    cfg.hwaAutocorrMaxLag_s = 1;
    cfg.hwaTargetUncertainty = 0.01;
    cfg.hwaCoverageFactor = 3;
    cfg.spectraWindowSamples = 16384;
    cfg.spectraOverlapFraction = 0.5;
    cfg.spectraPlotMaxFrequency_Hz = 5000;
    cfg.minimumPeakFrequency_Hz = 10;
    cfg.minimumPeakSeparation_Hz = 20;
    cfg.nDominantPeaks = 3;

    % PIV and HWA do not use the same vertical coordinate convention.
    % HWA y is reported relative to the zero-AoA trailing-edge traverse
    % reference. In the DaVis exports, the same zero-AoA reference is near
    % Y = 62 mm, and the positive vertical direction is opposite to the HWA
    % traverse sign. Use this only for profile comparisons, not for raw field
    % display.
    cfg.pivReferenceY_mm = 62;
    cfg.pivYSignForHwaComparison = -1;

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

    cfg.pivShortDtCase = struct( ...
        'aoa_deg', 15, ...
        'case_id', 'aoa15_short_dt', ...
        'processedFolder', 'Processed15dt');

    cfg.selfPiv.windowSizes_px = [16, 32, 64];
    cfg.selfPiv.deltaT_s = 74e-6;
    cfg.selfPiv.rawImagePath = fullfile(cfg.dataDir, 'PIV', 'FMT Results', ...
        'aoa15_final', 'B00001.tif');
    cfg.selfPiv.calibrationImagePath = fullfile(cfg.dataDir, 'PIV', 'FMT Results', ...
        'cal_final', 'B00001.tif');
end

function ensureOutputFolders(cfg)
    folders = {cfg.outputDir, cfg.tablesDir, cfg.figuresDir};
    for iFolder = 1:numel(folders)
        if ~isfolder(folders{iFolder})
            mkdir(folders{iFolder});
        end
    end
end

function [profiles, calibrationPoints, calibrationCurve, calibration] = processHwaLikeMain(cfg)
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

    calibration.signalUnique = signalUnique;
    calibration.uUnique = uUnique;

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
        records(end).y_relative_mm = position_mm - cfg.hwaTrailingEdgeY_mm;
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

function velocity = hwaSignalToVelocity(calibration, signal)
    velocity = interp1(calibration.signalUnique, calibration.uUnique, signal, 'pchip', 'extrap');
end

function [time, velocity] = readHwaVelocity(filePath, calibration)
    rawData = readmatrix(filePath);
    time = rawData(:, 1);
    signal = rawData(:, 2);
    valid = isfinite(time) & isfinite(signal);
    time = time(valid);
    signal = signal(valid);
    velocity = hwaSignalToVelocity(calibration, signal);
end

function [autocorrelation, convergence] = processHwaAutocorrelation(cfg, calibration)
    filePath = fullfile(cfg.hwaDir, 'correlationtest.txt');
    [time, velocity] = readHwaVelocity(filePath, calibration);
    sampleRate_hz = 1 / mean(diff(time));

    fluctuation = velocity - mean(velocity, 'omitnan');
    maxLagSamples = min(numel(fluctuation) - 1, round(cfg.hwaAutocorrMaxLag_s * sampleRate_hz));
    [rhoFull, lags] = xcorr(fluctuation, maxLagSamples, 'coeff');
    positive = lags >= 0;
    lag_s = lags(positive);
    lag_s = lag_s(:) / sampleRate_hz;
    rho = rhoFull(positive);
    rho = rho(:);

    firstZeroIndex = find(rho <= 0 & lag_s > 0, 1, 'first');
    if isempty(firstZeroIndex)
        firstZeroIndex = numel(rho);
    end

    integralTime_s = trapz(lag_s(1:firstZeroIndex), rho(1:firstZeroIndex));
    meanVelocity = mean(velocity, 'omitnan');
    rmsVelocity = sqrt(mean((velocity - meanVelocity).^2, 'omitnan'));
    requiredUncorrelatedSamples = (cfg.hwaCoverageFactor * rmsVelocity / ...
        (cfg.hwaTargetUncertainty * abs(meanVelocity)))^2;
    requiredSamplingTime_s = 2 * integralTime_s * requiredUncorrelatedSamples;

    autocorrelation = table(lag_s, rho, ...
        'VariableNames', {'lag_s', 'rho'});
    convergence = table(meanVelocity, rmsVelocity, sampleRate_hz, ...
        lag_s(firstZeroIndex), integralTime_s, cfg.hwaCoverageFactor, ...
        cfg.hwaTargetUncertainty, requiredUncorrelatedSamples, requiredSamplingTime_s, ...
        'VariableNames', {'mean_velocity_ms', 'rms_velocity_ms', 'sample_rate_hz', ...
        'first_zero_crossing_s', 'integral_time_s', 'coverage_factor', ...
        'target_uncertainty_fraction', 'uncorrelated_samples_required', ...
        'sampling_time_required_s'});
end

function [spectra, spectralPeaks] = processHwaSpectra(cfg, calibration, hwaProfiles)
    aoa15Profiles = hwaProfiles(hwaProfiles.aoa_deg == 15, :);
    [~, shearIndex] = max(aoa15Profiles.rms_velocity_ms);
    shearLayerPosition_mm = aoa15Profiles.position_mm(shearIndex);

    requests = table( ...
        ["AoA 0 trailing edge"; "AoA 5 trailing edge"; "AoA 15 shear layer"], ...
        [0; 5; 15], ...
        [cfg.hwaTrailingEdgeY_mm; cfg.hwaTrailingEdgeY_mm; shearLayerPosition_mm], ...
        'VariableNames', {'case_label', 'aoa_deg', 'position_mm'});

    spectra = table();
    spectralPeaks = table();

    for iRequest = 1:height(requests)
        aoa = requests.aoa_deg(iRequest);
        position_mm = requests.position_mm(iRequest);
        filePath = fullfile(cfg.hwaDir, sprintf('aoa%d_%d_10.txt', aoa, position_mm));
        [time, velocity] = readHwaVelocity(filePath, calibration);
        [frequency_hz, phi_uu, df_hz] = computeWelchSpectrum(time, velocity, ...
            cfg.spectraWindowSamples, cfg.spectraOverlapFraction);

        thisSpectrum = table( ...
            repmat(requests.case_label(iRequest), numel(frequency_hz), 1), ...
            repmat(aoa, numel(frequency_hz), 1), ...
            repmat(position_mm, numel(frequency_hz), 1), ...
            frequency_hz(:), phi_uu(:), repmat(df_hz, numel(frequency_hz), 1), ...
            'VariableNames', {'case_label', 'aoa_deg', 'position_mm', ...
            'frequency_hz', 'phi_uu', 'frequency_resolution_hz'});
        spectra = [spectra; thisSpectrum]; %#ok<AGROW>

        peaks = dominantFrequencies(frequency_hz, phi_uu, cfg.minimumPeakFrequency_Hz, ...
            cfg.spectraPlotMaxFrequency_Hz, cfg.minimumPeakSeparation_Hz, cfg.nDominantPeaks);
        peaks(end + 1:cfg.nDominantPeaks) = NaN;
        peakTable = table( ...
            repmat(requests.case_label(iRequest), cfg.nDominantPeaks, 1), ...
            repmat(aoa, cfg.nDominantPeaks, 1), ...
            repmat(position_mm, cfg.nDominantPeaks, 1), ...
            (1:cfg.nDominantPeaks).', peaks(:), ...
            'VariableNames', {'case_label', 'aoa_deg', 'position_mm', ...
            'peak_number', 'frequency_hz'});
        spectralPeaks = [spectralPeaks; peakTable]; %#ok<AGROW>
    end
end

function [frequency_hz, phi_uu, frequencyResolution_hz] = computeWelchSpectrum( ...
        time, velocity, windowSamples, overlapFraction)
    time = time(:);
    velocity = velocity(:);
    velocity = velocity - mean(velocity, 'omitnan');
    sampleRate_hz = 1 / mean(diff(time));
    windowSamples = min(windowSamples, numel(velocity));
    overlapSamples = floor(overlapFraction * windowSamples);
    window = hann(windowSamples, 'periodic');
    [phi_uu, frequency_hz] = pwelch(velocity, window, overlapSamples, ...
        windowSamples, sampleRate_hz);
    frequencyResolution_hz = sampleRate_hz / windowSamples;
end

function peakFrequencies_hz = dominantFrequencies(frequency_hz, phi_uu, ...
        minFrequency_hz, maxFrequency_hz, minSeparation_hz, nPeaks)
    valid = frequency_hz >= minFrequency_hz & frequency_hz <= maxFrequency_hz;
    candidate = false(size(valid));
    candidate(2:end - 1) = valid(2:end - 1) & ...
        phi_uu(2:end - 1) > phi_uu(1:end - 2) & ...
        phi_uu(2:end - 1) >= phi_uu(3:end);

    peakFrequencies = frequency_hz(candidate);
    peakEnergies = phi_uu(candidate);
    [~, order] = sort(peakEnergies, 'descend');
    peakFrequencies = peakFrequencies(order);

    peakFrequencies_hz = [];
    for iPeak = 1:numel(peakFrequencies)
        farEnough = isempty(peakFrequencies_hz) || ...
            all(abs(peakFrequencies(iPeak) - peakFrequencies_hz) >= minSeparation_hz);
        if farEnough
            peakFrequencies_hz(end + 1) = peakFrequencies(iPeak); %#ok<AGROW>
        end
        if numel(peakFrequencies_hz) == nPeaks
            break;
        end
    end
end

function [meanFields, rmsFields, instantaneousFields, profiles, processingSummary] = processPivLikeMain(cfg)
    meanFields = table();
    rmsFields = table();
    instantaneousFields = table();
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

        if pivCase.aoa_deg == 15
            instantPath = fullfile(caseDir, 'Overlap50MP3', 'B00001.dat');
            instantField = readDavisField(instantPath, true);
            instantaneousFields = fieldToTable(instantField, pivCase, cfg, 'instantaneous');
        end
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
        field.X(:), field.Y(:), pivYRelativeToHwa(cfg, field.Y(:)), ...
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
        field.X(:), field.Y(:), pivYRelativeToHwa(cfg, field.Y(:)), ...
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
        y, pivYRelativeToHwa(cfg, y), u, uRms, valid, ...
        'VariableNames', {'case_id', 'aoa_deg', 'processing', 'x_profile_mm', ...
        'y_mm', 'y_relative_mm', 'mean_velocity_ms', 'rms_velocity_ms', 'valid'});
end

function yRelative = pivYRelativeToHwa(cfg, y_mm)
    yRelative = cfg.pivYSignForHwaComparison * (y_mm - cfg.pivReferenceY_mm);
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

function [processingFields, windowFields, windowSummary, ensembleFields, deltaTFields, ...
        selfVsDavisFields, selfVsDavisSummary] = processPivParameterExports(cfg)
    aoa15Case = cfg.pivCases([cfg.pivCases.aoa_deg] == 15);
    processingSpecs = {
        '0% overlap, single pass', 'Overlap0SinglePassAvgStDev';
        '50% overlap, single pass', 'Overlap50SinglePassAvgStDev';
        '50% overlap, 3-pass', 'Overlap50MP3AvgStDev';
    };

    processingFields = table();
    for iSpec = 1:size(processingSpecs, 1)
        label = string(processingSpecs{iSpec, 1});
        folderName = processingSpecs{iSpec, 2};
        folderPath = fullfile(cfg.pivProcessedDir, aoa15Case.processedFolder, folderName);
        [avgFile, ~] = findAverageAndStdevFiles(folderPath);
        field = readDavisField(avgFile, true);
        thisTable = fieldToTable(field, aoa15Case, cfg, label);
        processingFields = [processingFields; thisTable]; %#ok<AGROW>
    end

    [windowFields, windowSummary] = computeWindowSizeStudy(cfg);
    ensembleFields = processPivEnsembleStudy(cfg, aoa15Case);
    deltaTFields = processPivDeltaTStudy(cfg, aoa15Case);
    [selfVsDavisFields, selfVsDavisSummary] = processSelfVsDavis(cfg);
end

function ensembleFields = processPivEnsembleStudy(cfg, aoa15Case)
    specs = {
        '10 sequential samples', 'Overlap50MP3Img10Inc1AvgStDev';
        '10 separated samples', 'Overlap50MP3Img10Inc10AvgStDev';
        '100 samples', 'Overlap50MP3AvgStDev';
    };
    ensembleFields = table();
    for iSpec = 1:size(specs, 1)
        label = string(specs{iSpec, 1});
        folderPath = fullfile(cfg.pivProcessedDir, aoa15Case.processedFolder, specs{iSpec, 2});
        [avgFile, ~] = findAverageAndStdevFiles(folderPath);
        field = readDavisField(avgFile, true);
        thisTable = fieldToTable(field, aoa15Case, cfg, label);
        ensembleFields = [ensembleFields; thisTable]; %#ok<AGROW>
    end
end

function deltaTFields = processPivDeltaTStudy(cfg, aoa15Case)
    specs = {
        'Original Delta t', aoa15Case, 'Overlap50MP3AvgStDev';
        'Short Delta t', cfg.pivShortDtCase, 'Overlap50MP3AvgStDev';
    };
    deltaTFields = table();
    for iSpec = 1:size(specs, 1)
        label = string(specs{iSpec, 1});
        pivCase = specs{iSpec, 2};
        folderPath = fullfile(cfg.pivProcessedDir, pivCase.processedFolder, specs{iSpec, 3});
        [avgFile, ~] = findAverageAndStdevFiles(folderPath);
        field = readDavisField(avgFile, true);
        thisTable = fieldToTable(field, pivCase, cfg, label);
        deltaTFields = [deltaTFields; thisTable]; %#ok<AGROW>
    end
end

function [selfVsDavisFields, selfVsDavisSummary] = processSelfVsDavis(cfg)
    pixelToMm = calibratePivImageLikeMain(cfg.selfPiv.calibrationImagePath);
    referenceField = readDavisField(fullfile(cfg.pivProcessedDir, ...
        'Processed15', 'Overlap0SinglePass', 'B00001.dat'), true);

    rawImage = imread(cfg.selfPiv.rawImagePath);
    halfImageHeight_px = floor(size(rawImage, 1) / 2);
    windowSize = 32;
    xOffset_mm = min(referenceField.X(:)) - 0.5 * windowSize * pixelToMm;
    yOffset_mm = max(referenceField.Y(:)) - ...
        (floor(halfImageHeight_px / windowSize) - 0.5) * windowSize * pixelToMm;

    selfField = computeFixedWindowPivLikeMain(cfg.selfPiv.rawImagePath, windowSize, ...
        cfg.selfPiv.deltaT_s, pixelToMm, cfg.pivStreamwiseSign);
    selfField = applyDavisMaskToWindowField(selfField, referenceField, xOffset_mm, yOffset_mm);

    davisU = cfg.pivStreamwiseSign * referenceField.Vx;
    davisV = referenceField.Vy;
    selfU = reshape(selfField.u_streamwise_ms, size(referenceField.X));
    selfV = reshape(selfField.v_normal_ms, size(referenceField.X));
    valid = referenceField.valid & reshape(selfField.valid, size(referenceField.X));
    du = selfU - davisU;
    dv = selfV - davisV;

    selfVsDavisFields = table( ...
        referenceField.X(:), referenceField.Y(:), davisU(:), davisV(:), ...
        selfU(:), selfV(:), du(:), dv(:), valid(:), ...
        'VariableNames', {'x_mm', 'y_mm', 'davis_u_streamwise_ms', ...
        'davis_v_normal_ms', 'self_u_streamwise_ms', 'self_v_normal_ms', ...
        'difference_u_ms', 'difference_v_ms', 'valid'});

    selfVsDavisSummary = table( ...
        sqrt(mean(du(valid).^2, 'omitnan')), sqrt(mean(dv(valid).^2, 'omitnan')), ...
        mean(abs(du(valid)), 'omitnan'), mean(abs(dv(valid)), 'omitnan'), nnz(valid), ...
        'VariableNames', {'rmse_u_ms', 'rmse_v_ms', 'mean_abs_u_ms', ...
        'mean_abs_v_ms', 'n_valid'});
end

function [windowFields, windowSummary] = computeWindowSizeStudy(cfg)
    pixelToMm = calibratePivImageLikeMain(cfg.selfPiv.calibrationImagePath);
    referenceMaskField = readDavisField(fullfile(cfg.pivProcessedDir, ...
        'Processed15', 'Overlap50MP3AvgStDev', 'B00001.dat'), true);
    rawImage = imread(cfg.selfPiv.rawImagePath);
    halfImageHeight_px = floor(size(rawImage, 1) / 2);
    baseMaskWindow_px = 16;
    xOffset_mm = min(referenceMaskField.X(:)) - 0.5 * baseMaskWindow_px * pixelToMm;
    yOffset_mm = max(referenceMaskField.Y(:)) - ...
        (floor(halfImageHeight_px / baseMaskWindow_px) - 0.5) * baseMaskWindow_px * pixelToMm;

    windowFields = table();
    windowSummary = table();

    for windowSize = cfg.selfPiv.windowSizes_px
        field = computeFixedWindowPivLikeMain(cfg.selfPiv.rawImagePath, windowSize, ...
            cfg.selfPiv.deltaT_s, pixelToMm, cfg.pivStreamwiseSign);
        field = applyDavisMaskToWindowField(field, referenceMaskField, xOffset_mm, yOffset_mm);
        windowFields = [windowFields; field]; %#ok<AGROW>

        valid = field.valid;
        row = table( ...
            windowSize, height(field), ...
            nnz(valid), nnz(valid) / height(field), ...
            mean(field.u_streamwise_ms(valid), 'omitnan'), ...
            std(field.u_streamwise_ms(valid), 'omitnan'), ...
            mean(field.speed_ms(valid), 'omitnan'), ...
            'VariableNames', {'window_size_px', 'n_vectors', ...
            'n_valid_vectors', 'valid_fraction', 'mean_u_streamwise_ms', ...
            'spatial_std_u_streamwise_ms', 'mean_speed_ms'});
        windowSummary = [windowSummary; row]; %#ok<AGROW>
    end
end

function field = applyDavisMaskToWindowField(field, referenceField, xOffset_mm, yOffset_mm)
    imageX = field.x_mm;
    imageY = field.y_mm;
    field.image_x_mm = imageX;
    field.image_y_mm = imageY;
    field.x_mm = imageX + xOffset_mm;
    field.y_mm = imageY + yOffset_mm;

    refX = mean(referenceField.X, 1, 'omitnan');
    refY = mean(referenceField.Y, 2, 'omitnan');
    valid = false(height(field), 1);
    for iPoint = 1:height(field)
        [~, xIdx] = min(abs(refX - field.x_mm(iPoint)));
        [~, yIdx] = min(abs(refY - field.y_mm(iPoint)));
        valid(iPoint) = referenceField.valid(yIdx, xIdx);
    end
    field.valid = valid;
end

function pixelToMm = calibratePivImageLikeMain(calibrationImagePath)
    imread(calibrationImagePath);
    x1 = 330;
    x2 = 1050;
    y1 = 1160;
    y2 = 1150;
    pixelToMm = 80 / sqrt((x2 - x1)^2 + (y2 - y1)^2);
end

function tableOut = computeFixedWindowPivLikeMain(imagePath, windowSize, deltaT_s, pixelToMm, streamwiseSign)
    [img1, img2] = splitPivImageWithoutToolbox(imagePath);
    [height, width] = size(img1);
    nRows = floor(height / windowSize);
    nCols = floor(width / windowSize);
    nVectors = nRows * nCols;

    rowIndex = zeros(nVectors, 1);
    colIndex = zeros(nVectors, 1);
    x_mm = zeros(nVectors, 1);
    y_mm = zeros(nVectors, 1);
    vx = zeros(nVectors, 1);
    vy = zeros(nVectors, 1);

    k = 0;
    for iCol = 1:nCols
        for iRow = 1:nRows
            k = k + 1;
            xPixels = (iCol - 1) * windowSize + 1:iCol * windowSize;
            yPixels = (iRow - 1) * windowSize + 1:iRow * windowSize;

            window1 = double(img1(yPixels, xPixels));
            window2 = double(img2(yPixels, xPixels));
            window1 = window1 - mean(window1(:));
            window2 = window2 - mean(window2(:));

            corrMap = xcorr2(window1, window2);
            [~, vectorizedIndex] = max(corrMap(:));
            [peakY, peakX] = ind2sub(size(corrMap), vectorizedIndex);

            dx_px = peakX - windowSize;
            dy_px = peakY - windowSize;

            rowIndex(k) = iRow;
            colIndex(k) = iCol;
            x_mm(k) = ((iCol - 0.5) * windowSize) * pixelToMm;
            y_mm(k) = ((nRows - iRow + 0.5) * windowSize) * pixelToMm;
            vx(k) = dx_px * pixelToMm / deltaT_s / 1000;
            vy(k) = dy_px * pixelToMm / deltaT_s / 1000;
        end
    end

    u = streamwiseSign * vx;
    tableOut = table( ...
        repmat(windowSize, nVectors, 1), rowIndex, colIndex, x_mm, y_mm, ...
        vx, vy, u, vy, hypot(u, vy), ...
        'VariableNames', {'window_size_px', 'row_index', 'col_index', ...
        'x_mm', 'y_mm', 'vx_image_ms', 'vy_image_ms', ...
        'u_streamwise_ms', 'v_normal_ms', 'speed_ms'});
end

function [img1, img2] = splitPivImageWithoutToolbox(imagePath)
    img = imread(imagePath);
    halfHeight = floor(size(img, 1) / 2);
    img2 = img(1:halfHeight, :);
    img1 = img(halfHeight + 1:2 * halfHeight, :);
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
