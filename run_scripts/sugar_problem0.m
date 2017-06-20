%import sugerdataset and cheb function
addpath('../../../datasets/sugar/data')
addpath('../../cheb')

%path to origial no nonsense implementation
addpath('../lvq')

%mode 1 is global mode 4 is local matrices
mode = 4;

% read the data from the files
dataset_filename = 'sugar_neoVNIR1600.csv'
[spectra, labels, wavelengths] = loadDataset(dataset_filename);
[spectra_val, labels_val, wavelengths_val] = loadDataset('sugar_neoVNIR1600_val.csv');

%relable to problem
[newData, newLabels] = reformulateDataset(spectra, labels, 0);

%for number_of_coefficients=20:5:50
    % parse them to chebysev polynomials
    number_of_coefficients = 20;
    [coefficients, transformationMatrix, backtransformationMatrix] = chebyshev(newData, number_of_coefficients);
    strcat(num2str(number_of_coefficients), 'coefficients')

    %lvq validate with global matrices
    [gmlvq_mean,roc_val,lcurves_mean,lcurves_std,param_set,f_measures] = run_validation(coefficients,newLabels,150, 50, 10, [1,2,3,4,5,6,7,8,9],mode);
%end;

