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

figure(1);
title('UCR Plane',...
     'FontName','LucidaSans', 'FontWeight','bold');
hold on;

amount_of_classes = length(unique(label));

%plot mean and variance for every class
for class_nr = 1:amount_of_classes;
   class_data = x(label == class_nr, :);
   class_mean = mean(class_data);
   class_variance= var(class_data);
   plot(1:size(x,2), class_mean);
   xlabel('x-axis');
end;
