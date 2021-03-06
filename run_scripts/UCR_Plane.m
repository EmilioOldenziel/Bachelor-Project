%import UCRdataset and cheb function
addpath('../../cheb')
addpath('../../../datasets/UCR')
load('UCR_Plane.mat')

%path to origial no nonsense imlementation
addpath('../lvq')

%mode 1 is global mode 4 is local matrices
mode = 4;

number_of_coefficients = 20;
% parse them to chebysev polynomials
[coefficients, transformationMatrix, backtransformationMatrix] = chebyshev(x, number_of_coefficients);
strcat(num2str(number_of_coefficients), 'coeficients')
%lvq validate with global matrices
[gmlvq_mean,roc_val,lcurves_mean,lcurves_std,param_set] = run_validation(coefficients,label,50, 10, 15, [1,2,3,4,5,6,7], mode);

%save the figures in the results folder
dataset_name = 'UCR_Plane';
dir_name = strcat(dataset_name, num2str(number_of_coefficients), 'coef', ...
'UCR', 'etam', num2str(param_set.etam),'etap', num2str(param_set.etap));
dir_path = strcat('../../../results/', dir_name);
mkdir('../../../results/', dir_name);
save_figures(dir_path);