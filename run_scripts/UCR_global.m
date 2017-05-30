%import UCRdataset and cheb function
addpath('../../cheb')
addpath('../../../datasets/UCR')
load('UCR_Plane.mat')

%path to origial no nonsense imlementation
addpath('../../no-nonsense-gmlvq-v2-3/')

%parfor number_of_coefficients=20:3:50
    number_of_coefficients = 20;
    % parse them to chebysev polynomials
    [coefficients, transformationMatrix, backtransformationMatrix] = chebyshev(x, number_of_coefficients);
    strcat(num2str(number_of_coefficients), 'coeficients')
    %lvq validate with global matrices
    [gmlvq_mean,roc_val,lcurves_mean,lcurves_std,param_set] = run_validation(coefficients,label,50, 1, 10);
    
    %save the figures in the results folder
    dataset_name = 'UCR_Plane';
    dir_name = strcat(dataset_name, num2str(number_of_coefficients), 'coef', ...
    'UCR_global', 'etam', num2str(param_set.etam),'etap', num2str(param_set.etap));
    dir_path = strcat('../../../results/', dir_name);
    mkdir('../../../results/', dir_name);
    save_figures(dir_path);
    etap = etap/4;
%end