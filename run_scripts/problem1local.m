%import sugerdataset and cheb function
addpath('../../../datasets/sugar/data')
addpath('../../cheb')

%path to lvq with local matrices
addpath('../lvq')

% read the data from the files
[spectra, labels, wavelengths] = loadDataset('sugar_neoVNIR1600.csv');
[spectra_val, labels_val, wavelengths_val] = loadDataset('sugar_neoVNIR1600_val.csv');

% parse them to chebysev polynomials
number_of_coefficients = 20;
[coefficients, transformationMatrix, backtransformationMatrix] = chebyshev(spectra, number_of_coefficients);

%reduce from 9 to 3 classes #problem 1
problem1labels = labels
for i=1:length(labels)
    if problem1labels(i) < 5
        problem1labels(i) = 1;
    elseif problem1labels(i) > 6
        problem1labels(i) = 3;
    else
        problem1labels(i) = 2;
    end
end

%lvq validation with local matrices
[gmlvq_mean,roc_val,lcurves_mean,lcurves_std,param_set] = run_validation(coefficients,problem1labels,50, 10,10);

%save the figures in the results folder
dataset_name = strrep(dataset_filename,'.csv','');
dir_name = strcat(dataset_name, num2str(number_of_coefficients), 'coef', 'problem1local');
dir_path = strcat('../../../results/', dir_name);
mkdir('../../../results/', dir_name);
save_figures(dir_path);