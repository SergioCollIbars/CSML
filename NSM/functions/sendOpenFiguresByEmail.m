function sendOpenFiguresByEmail(recipientEmail, subjectLine, consoleText)
% Create a temporary folder
    tempFolder = fullfile(tempdir, 'matlab_figures');
    if ~exist(tempFolder, 'dir')
        mkdir(tempFolder);
    end

    % Get handles to all open figures, sorted by creation number
    figs = findall(0, 'Type', 'figure');
    if isempty(figs)
        warning('No open figures found.');
        return;
    end
    figNumbers = arrayfun(@(f) f.Number, figs);
    [~, sortIdx] = sort(figNumbers);
    figs = figs(sortIdx);

    % Filenames
    pdfFile = fullfile(tempFolder, 'AllFigures.pdf');
    zipFile = fullfile(tempFolder, 'LastTwoFigures.zip');

    % Export all figures to a single PDF
    for k = 1:numel(figs)
        figure(figs(k));
        exportgraphics(figs(k), pdfFile, 'Append', true);
    end

    % Save only the last two figures as .fig files
    lastTwo = figs(max(end-1,1):end);
    figPaths = cell(1, numel(lastTwo));
    for k = 1:numel(lastTwo)
        figFileName = fullfile(tempFolder, sprintf('Figure%d.fig', lastTwo(k).Number));
        savefig(lastTwo(k), figFileName);
        figPaths{k} = figFileName;
    end

    % Zip the .fig files
    zip(zipFile, figPaths);

    % Send email with PDF and ZIP attached
    attachments = {pdfFile, zipFile};
    sendmail(recipientEmail, subjectLine, consoleText, attachments);

    % Clean up
    delete(pdfFile);
    delete(zipFile);
    cellfun(@delete, figPaths);
    rmdir(tempFolder, 's');

    disp('--------------------------------------------------------')
    disp('Email sent!')
end

