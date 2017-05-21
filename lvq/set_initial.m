
% initialization of prototypes close to class conditional means
% small random displacements to break ties    

function [proti,omi] = set_initial(fvec,lbl,plbl,mode,rndinit)

% input
% fvec  :         feature vectors
% lbl   :         data labels
% plbl  :         prototype labels
% mode  :         0,1 for full matrix       (GMLVQ)
%                 2 for diagonal matrix     (GRLVQ)
%                 3 for prop. to identity   (GLVQ)
%                 4 for local matrices      (LGMLVQ)

  ndim  = size(fvec,2);         % dimension of feature vectors
  nprots= length(plbl);         % total number of prototypes
  
  
  for ic=1:nprots; % compute class-conditional means
      proti(ic,:) = mean(fvec(lbl==plbl(ic),:),1) ;
  end; 
      % displace randomly from class-conditional means
        proti= proti.*(0.99+0.02*rand(size(proti))); 
      % to do: run k-means per class

  %intitialize local relevance matrices per prototype and randomly
  if (mode == 4);
    omegas_initialized = zeros(ndim,ndim,nprots);
    for iom=1:nprots;
      if(rndinit==1)
        omega_initialized = rand(ndim)-0.5;
        omega_initialized = omega_initialized'*omega_initialized;
      else
        omega_initialized = eye(ndim);
      end;
      % normalization, Trace(Lambda)=1
      omega_initialized=omega_initialized/sqrt(sum(sum(omega_initialized.^2)));
      omegas_initialized(:,:,iom) = omega_initialized;
    end
    omi = omegas_initialized;
  else
    % (global) matrix initialization, identity or random
    omi=eye(ndim);          % works for all values of mode if rndinit == 0
    if (mode ~=3 && rndinit==1);  % does not apply for mode==3 (GLVQ)
      omi= rand(ndim)-0.5;       
      omi= omi'*omi;             % square symmetric 
                                  %  matrix of uniform random numbers       
    end; 
    if (mode ==2);
      omi=diag(diag(omi));       % restrict to diagonal matrix 
    end;
    omi=omi/sqrt(sum(sum(omi.^2)));    % normalization, Trace(Lambda)=1
  end
end

