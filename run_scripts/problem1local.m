%import sugerdataset and cheb function
addpath('../../../datasets/sugar/data')
addpath('../../cheb')

%path to origial no nonsense implementation
addpath('../lvq')

% read the data from the files
dataset_filename = 'sugar_neoVNIR1600.csv'
[spectra, labels, wavelengths] = loadDataset(dataset_filename);
[spectra_val, labels_val, wavelengths_val] = loadDataset('sugar_neoVNIR1600_val.csv');

%relable to problem
[newData, newLabels] = reformulateDataset(spectra, labels, 1);
etam = 2;
etap = 1;

%for number_of_coefficients=20:5:50
    %while(etam>(3/1000));
        %while(etap>(3/10000));
        % parse them to chebysev polynomials
        number_of_coefficients = 20;
        [coefficients, transformationMatrix, backtransformationMatrix] = chebyshev(newData, number_of_coefficients);
        strcat(num2str(number_of_coefficients), 'coefficients')
        
        %lvq validate with global matrices
        [gmlvq_mean,roc_val,lcurves_mean,lcurves_std,param_set] = run_validation(coefficients,newLabels,50, 10, 10, [1,2,3],etam, etap);

        %save the figures in the results folder
        dataset_name = strrep(dataset_filename,'.csv','');
        dir_name = strcat(dataset_name, num2str(number_of_coefficients), 'coef', ...
        'problem1local', 'etam', num2str(param_set.etam),'etap', num2str(param_set.etap));
        dir_path = strcat('../../../results/', dir_name);
        mkdir('../../../results/', dir_name);
        save_figures(dir_path);
        %etap = etap/4;
        %end;
    %etam = etam/4;
    %end;
%end;

