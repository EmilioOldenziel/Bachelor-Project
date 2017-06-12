addpath('../../../datasets/sugar/data')

% read the data from the files
[spectra, labels, wavelengths] = loadDataset('sugar_neoVNIR1600.csv');
[spectra_val, labels_val, wavelengths_val] = loadDataset('sugar_fieldspec_val.csv');

figure(1);
for problem_nr = 0:4;
    
    subplot(2,3,problem_nr+1);
    title(strcat('Problem', int2str(problem_nr)),...
         'FontName','LucidaSans', 'FontWeight','bold');
    hold on;
    % relabel data to correct class
    [newData, newLabels] = reformulateDataset(spectra, labels, problem_nr);
    
    [coefficients, transformationMatrix, backtransformationMatrix] = chebyshev(newData, 20);

    amount_of_classes = length(unique(newLabels));

    %plot mean and variance for every class
    for class_nr = 1:amount_of_classes;
       class_data = spectra(newLabels == class_nr, :);
       class_mean = mean(class_data);
       class_variance= var(class_data);
       errorbar(wavelengths, class_mean, class_variance);
       xlabel('wavelength (nm)');
    end;
    
    hold off;

end;
