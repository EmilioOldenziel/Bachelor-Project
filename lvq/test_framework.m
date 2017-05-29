
%path to origial no nonsense imlementation
addpath('../lvq')

fvec = (checkerboard(20) > 0.5)+(rand(8*20,8*20)/20);
labels = (fvec(:,1) > 0.5)+1;
fvec=fvec/sqrt(sum(sum(fvec.^2)));
trace(fvec)

%lvq validate with global matrices
[gmlvq_mean,roc_val,lcurves_mean,lcurves_std,param_set] = run_validation(fvec,labels,50, 2, 10);

