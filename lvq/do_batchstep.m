% perform a single step of batch gradient descent GMLVQ
% with given step size for matrix and prototype updates (input parameter)
% only for one global quadratic omega matrix, potentially diagonal (mode=2)
% optional: null-space correction for full matrix only (mode=1)

function [prot,omat]= do_batchstep (fvec,lbl,proti,plbl,...
                               omegai,etap,etam,mu,mode)

% input: 
% fvec             nvf feature vectors of dim. ndim  fvec(1:nfv,1:ndim);
% lbl              training labels  lbl(1:nfv);
% proti            prototypes before the step
% plbl             prototype labels
% omegai           global matrix or local_matrices before the step
% etap,etam        prototype/matrix learning rate
% mu               mu>0 controls penalty for singular relevance matrix
% mode             0 for full, 1 for matrix with nullspace correction
%                  2 for diagonal matrix only       (GRLVQ)
%                  3 for normalized identity matrix (GLVQ)

% output:
% prot             prototypes after update
% omat             omega matrix after update

ndim = size(fvec,2);    % dimension of feature vectors
nfv = length(lbl);      % number of feature vectors (training samples)
cls = unique(lbl);      % set of class labels
np = size(proti,1);     % number of prototypes

omat = omegai; %omega before step

if(mode ==4)
  lambda = zeros(ndim,ndim,np);
  for iom=1:np
    lambda(:,:,iom) = omat(:,:,iom)'*omat(:,:,iom);
  end
else
    lambda = omat'*omat; % lambda before step
end

% omat=sqrtm(lambda);
prot=proti;                            % prototypes before step

     chp = 0*prot;
     if(mode==4)
        chms_f1 = zeros(ndim,ndim);
        chms_f2 = zeros(ndim,ndim);
     else
       chm = 0*omat;       % initialize change of prot,omega
     end
     for i= 1:nfv;                     % loop through all training examples
       fvi = fvec(i,:); lbi = lbl(i);  % actual example
  
     % calculate squared distances to all prototypes
      dist = nan(np,1);                 % define squared distances
         for j=1:np;                       % distances from all prototypes
             if(mode ==4)
                 dist(j) = (norm(omat(:,:,j)*(fvi-prot(j,:))' ))^2;
            else
                 dist(j) = (norm(omat*(fvi-prot(j,:))' ))^2;
            end
         end;
   
    % find the two winning prototypes
      correct =   find (plbl == lbi);     % all correct prototype indices
      incorrect = find (plbl ~= lbi);     % all wrong   prototype indices
      
      [dJ, JJ] = min (dist (correct));    % correct winner
      [dK, KK] = min (dist (incorrect));  % wrong   winner

      jwin = correct (JJ); kwin = incorrect (KK);   % winner indices 
      wJ = prot(jwin,:);   wK = prot(kwin,:);       % winning prototypes

      % GMLVQ prototype update for one example fvi
     
      DJ = (fvi-wJ)';  DK = (fvi-wK)';        % displacement vectors 
      norm_factor = (dJ+dK)^2;                % denominator of prefactor
      if(mode ==4)
        dwJ = -(dK/norm_factor)*lambda(:,:,jwin)*DJ;      % change of correct winner
        dwK =  (dJ/norm_factor)*lambda(:,:,kwin)*DK;      % change of incorrect winner

        f1 = ( dK/norm_factor)*(omat(:,:,jwin)*DJ)*DJ';   % term 1 of matrix change
	    f2 = (-dJ/norm_factor)*(omat(:,:,kwin)*DK)*DK';   % term 2 of matrix change
      else
        dwJ = -(dK/norm_factor)*lambda*DJ;    % change of correct winner
        dwK =  (dJ/norm_factor)*lambda*DK;    % change of incorrect winner

        f1 = ( dK/norm_factor)*(omat*DJ)*DJ';   % term 1 of matrix change
	    f2 = (-dJ/norm_factor)*(omat*DK)*DK';   % term 2 of matrix change
      end
      % matrix update, single (global) matrix omat for one example
	  
      % negative gradient update added up over examples
      chp(jwin,:) = chp(jwin,:) - dwJ';  % correct   winner summed update
      chp(kwin,:) = chp(kwin,:) - dwK';  % incorrect winner summed update
      if(mode==4)
          chms_f1 = chms_f1 - f1;
          chms_f2 = chms_f2 - f2; 
      else
	      chm = chm - (f1 + f2);             % matrix summed update
      end
    end; % end of one loop through (sum over) all examples

    if(mode~=4)
      % singularity control: add  derivative of penalty term times mu 
      if(mu>0)
          chm = chm + mu*pinv(omat');
      end;
      
      % compute normalized gradient updates (length 1)
      % separate nomralization for prototypes and the matrix
      % computation of actual changes, diagonal matrix imposed here if nec.
      n2chw = 0;               % zero initial value of sum
      for ni=1:np;             % loop through (sum over) set of prototypes
          n2chw = n2chw + dot(chp(ni,:),chp(ni,:)); 
      end;                     % total length of concatenated prototype vec.
    
      if (mode==2)             % if diagonal matrix used only 
          chm=diag(diag(chm)); % reduce to diagonal changes
      end;
      n2chm = sum(sum(chm.^2));% total 'length' of matrix update
        
      % final, normalized gradient updates after 1 loop through examples
      prot = prot  + etap * chp/sqrt(n2chw);
      omat = omat  + etam * chm/sqrt(n2chm);
      
      if (mode==2)                   % if diagonal matrix only 
          omat = diag(diag(omat));  % reduce to diagonal matrix
      end;                           % probably obsolete as chm diagonal
        
      %  nullspace correction using Moore Penrose pseudo-inverse
      if (mode==1)
          xvec=[fvec;prot];               % concat. protos and fvecs
          omat= ((omat*xvec')*pinv(xvec')); % corrected omega matrix 
      end;
        
      if (mode==3);            
          omat=eye(ndim);  % reset to identity regardless of gradients
      end;
      
      % normalization of omega, corresponds to Trace(lambda) = 1
      omat = omat / sqrt( sum(sum(omat.^2)) );
    else
      % mode 4
      if (mu>0)
        chms_f1 = chms_f1 * mu*pinv(omat(:,:,jwin))';
        chms_f2 = chms_f2 * mu*pinv(omat(:,:,kwin))';
      end

      n2chw = 0;               % zero initial value of sum
      for ni=1:np;             % loop through (sum over) set of prototypes
          n2chw = n2chw + dot(chp(ni,:),chp(ni,:)); 
      end;
      n2chm_f1 = sum(sum(chms_f1.^2));% total 'length' of matrix update
      n2chm_f2 = sum(sum(chms_f2.^2));% total 'length' of matrix update

      prot = prot  + etap * chp/sqrt(n2chw);

      xvec=[fvec;prot]; % concat. protos and fvecs

      %update, correct and normalize winner matrix
      omat(:,:,jwin) = omat(:,:,jwin) + etam * chms_f1/sqrt(n2chm_f1);
      omat(:,:,jwin)= ((omat(:,:,jwin)*xvec')*pinv(xvec')); % corrected omega matrix
      omat(:,:,jwin) = omat(:,:,jwin) / sqrt( sum(sum(omat(:,:,jwin).^2)));


      %update, correct and normalize loser matrix
      omat(:,:,kwin) = omat(:,:,kwin) + etam * chms_f2/sqrt(n2chm_f2);
      omat(:,:,kwin)= ((omat(:,:,kwin)*xvec')*pinv(xvec')); % corrected omega matrix
      omat(:,:,kwin) = omat(:,:,kwin) / sqrt( sum(sum(omat(:,:,kwin).^2)));

      %update, correct and normalize all loser matrices
      % for iom=1:length(incorrect)
      %   omat(:,:,incorrect(iom)) = omat(:,:,incorrect(iom)) + etam * f2/sqrt(n2chm_f2);
      %   omat(:,:,incorrect(iom))= ((omat(:,:,incorrect(iom))*xvec')*pinv(xvec')); % corrected omega matrix 
      %   omat(:,:,incorrect(iom)) = omat(:,:,incorrect(iom)) / sqrt( sum(sum(omat(:,:,incorrect(iom)).^2)));
      % end


    end
    
end  % one full, normalized gradient step performed, return omat and prot




