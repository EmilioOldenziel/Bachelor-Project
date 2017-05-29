%import sugerdataset and cheb function
addpath('../../cheb')

%path to origial no nonsense imlementation
addpath('../../no-nonsense-gmlvq-v2-3/')

%parfor number_of_coefficients=20:3:50
    number_of_coefficients = 20;
    % parse them to chebysev polynomials
    [coefficients, transformationMatrix, backtransformationMatrix] = chebyshev(x, number_of_coefficients);
    strcat(num2str(number_of_coefficients), 'coeficients')
    %lvq validate with global matrices
    [gmlvq_mean,roc_val,lcurves_mean,lcurves_std,param_set] = run_validation(coefficients,label,50, 1, 10);
    
%end