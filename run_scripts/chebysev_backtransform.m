addpath('../../../datasets/sugar/data')

% read the data from the files
[spectra, labels, wavelengths] = loadDataset('sugar_neoVNIR1600.csv');
[spectra_val, labels_val, wavelengths_val] = loadDataset('sugar_fieldspec_val.csv');
%spectra = spectra / sqrt( sum(sum(spectra.^2)));

figure(1);
for problem_nr = 0:4;
    
    subplot(2,3,problem_nr+1);
    title(strcat('Problem', int2str(problem_nr)),...
         'FontName','LucidaSans', 'FontWeight','bold');
    hold on;
    % relabel data to correct class
    [newData, newLabels] = reformulateDataset(spectra, labels, problem_nr);
    
    [coefficients, transformationMatrix, backtransformationMatrix] = chebyshev(newData, 20);
    
    backtransform = zeros(1125,20);
    for fvec=1:size(coefficients,1);
        backtransform(fvec,:) = backtransformationMatrix*coefficients(fvec,:)';
    end;
    
    amount_of_classes = length(unique(newLabels));

    %plot mean and variance for every class
    for class_nr = 1:amount_of_classes;
       class_data = backtransform(newLabels == class_nr, :);
       class_mean = mean(class_data);
       class_variance= var(class_data);
       length(class_mean)
       length(class_variance)
       errorbar(1:20, class_mean, class_variance);
       xlabel('wavelength (nm)');
    end;
    
    hold off;

end;
