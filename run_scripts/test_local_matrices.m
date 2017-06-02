
%intitialize local relevance matrices per prototype and randomly
ndim = 4;
nprots = 3; 

omegas_initialized = zeros(ndim,ndim,nprots);
for iom=1:nprots;
  omega_initialized = eye(ndim);
  % normalization, Trace(Lambda)=1
  omega_initialized=omega_initialized/sqrt(sum(sum(omega_initialized.^2)));
  omegas_initialized(:,:,iom) = omega_initialized;
end

for iom=1:nprots;
    lambda(:,:,iom) = omegas_initialized(:,:,iom)'*omegas_initialized(:,:,iom);
end

lambda

lambda_mean = zeros(ndim,ndim,nprots);
lambda_std = zeros(ndim,ndim,nprots);
for iom=1:nprots
 lambda_mean(:,:,iom) = squeeze(mean(lambda(:,:,iom,:),4));
 lambda_std(:,:,iom)  = sqrt(squeeze(mean(lambda(:,:,iom,:).^2,4))-lambda_mean(:,:,iom).^2);
end

if lambda == lambda_mean;
    1
else
    0
end

