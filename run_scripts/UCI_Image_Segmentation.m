%import sugerdataset and cheb function
addpath('../../../datasets/UCIimagesegmentation')

%path to origial no nonsense implementation
addpath('../lvq')
dataset_filename = 'segmentation_brick_window.csv';
data = load(dataset_filename);

labels = data(:,1)
data(:,1) = [];

%remove colums 3-5
data(:,3) = [];
data(:,3) = [];
data(:,3) = [];

%mode 1 is global mode 4 is local matrices
mode = 1;

%lvq validate with global matrices
[gmlvq_mean,roc_val,lcurves_mean,lcurves_std,param_set] = run_validation(data,labels,50, 10, 10, [1,2], mode);

%save the figures in the results folder
dataset_name = strrep(dataset_filename,'.csv','');
dir_name = strcat(dataset_name, 'UCI_image_seg_local', 'etam', num2str(param_set.etam),'etap', num2str(param_set.etap));
dir_path = strcat('../../../results/', dir_name);
mkdir('../../../results/', dir_name);
save_figures(dir_path);

