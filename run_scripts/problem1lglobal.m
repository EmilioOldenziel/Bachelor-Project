%import sugerdataset and cheb function
addpath('../../../datasets/sugar/data')
addpath('../../cheb')

%path to origial no nonsense imlementation
addpath('../../no-nonsense-gmlvq-v2-3/')

% read the data from the files]
dataset_filename = 'sugar_neoVNIR1600.csv'
[spectra, labels, wavelengths] = loadDataset(dataset_filename);
[spectra_val, labels_val, wavelengths_val] = loadDataset('sugar_neoVNIR1600_val.csv');

%relable to problem
[newData, newLabels] = reformulateDataset(spectra, labels, 1);

% parse them to chebysev polynomials
number_of_coefficients = 20;
[coefficients, transformationMatrix, backtransformationMatrix] = chebyshev(spectra, number_of_coefficients);

%lvq validate with global matrices
[gmlvq_mean,roc_val,lcurves_mean,lcurves_std,param_set] = run_validation(coefficients,newLabels,50, 10, 10);

%save the figures in the results folder
dataset_name = strrep(dataset_filename,'.csv','');
dir_name = strcat(dataset_name, num2str(number_of_coefficients), 'coef', 'problem1global');
dir_path = strcat('../../../results/', dir_name);
mkdir('../../../results/', dir_name);
save_figures(dir_path);