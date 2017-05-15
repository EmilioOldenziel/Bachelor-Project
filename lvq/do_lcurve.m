
% gmlvq with global matrix only, square Omega, potentially diagonal
% batch gradient descent with step size control
% following a modified Papari procedure (unpublished) 
% perform training based using fvec, lbl
% evaluate performance measures after each step
% with respect to training set and validation set

function [w,omega,cftra,tetra,cwtra,auctra,cfval,teval,cwval,aucval]= ...
   do_lcurve(fvec,lbl,fvecval,lblval,plbl,totalsteps, mode)
% input arguments    
% fvec, lbl              training data, feature vectors and labels
% fvecval, lblval        validation data, feature vectors and labels
% plbl: labels assigned to prototypes, also specifies number per class
% e.g. plbl=[1,1,2,2,3] for 5 prototypes with labels 1,2,3  
% mode 
  % 0 for matrix without regularization
  % 1 for matrix with null-space correction
  % 2 for diagonal updates only  
% rndinit
  % 0 for initialization of relevances as identity matrix 
  % 1 for randomized initialization of relevance matrix
% totalsteps: number of batch gradient steps to be performed

% general algorithm settings and parameters of the Papari procedure
[~,~,mode,rndinit, m_update_step_size, prototype_update_step_size, mu, decfac, incfac, n_original] =...
                                        set_parameters(fvec, mode);
% m_update_step_size:             step size of matrix updates
% prototype_update_step_size:     step size of prototype updates
% mu   :                          control parameter of penalty term
% decfac:                         factor for decrease of step sizes for oscillatory behavior
% incfac:                         factor for increase of step sizes for smooth behavior
% n_original:                     number of copies in Papari procedure

nfv = size(fvec,1);          % number of feature vectors in training set
nfvval = size(fvecval,1);    % number of feature vectors for validation
ndim = size(fvec,2);         % dimension of feature vectors
number_of_classes = length(unique(lbl));  % number of classes 
number_of_prototypes = length(plbl);       % total number of prototypes
          
% initialize prototypes and omega
  [prototypes_initial,omega_intitial] =  set_initial(fvec,lbl,plbl,mode,rndinit);
  prototypes=prototypes_initial;  omega=omega_intitial;   % initial values

  for i=1:number_of_prototypes
    omega_intitial{i}
  end

%copies of prototypes are stored in prototypes_original
prototypes_original = zeros(n_original,size(prototypes,1),size(prototypes,2));

if(mode==4)
  omega_original = cell(n_original,size(omega,2));
else
%  copies of omegas are stored in omega_original
  omega_original   = zeros(n_original,size(omega,1) , size(omega,2));
end
  
   % learning curves, perfomance w.r.t. training and validation set
 cftra =  NaN(totalsteps,1); cfval = cftra; % cost fucnction and equivalent
 tetra =  NaN(totalsteps,1); teval = tetra; % total errors
 cwtra =  NaN(totalsteps,number_of_classes); cwval = cwtra; % class-wise errors
 auctra = NaN(totalsteps,1); aucval = auctra;   % auc(roc) 
  
 
  % perform the first n_original steps of gradient descent and compute
  % performance
  for inistep=1: n_original;
      [prototypes,omega]= do_batchstep(fvec,lbl,prototypes,plbl,omega,prototype_update_step_size,m_update_step_size,mu,mode); 
      prototypes_original(inistep,:,:)= prototypes;
      omega_original  (inistep,:,:)= omega;
      
      omega=omega/sqrt(sum(sum(omega.*omega))); 
      % compute costs without penalty term here
      [costf,~,marg,score]    = compute_costs(fvec,lbl,prototypes,plbl,omega,0);
      [costval,~,margval,scoreval]= ...
                         compute_costs(fvecval,lblval,prototypes,plbl,omega,0); 
      cftra(inistep)= costf; cfval(inistep)= costval;
      
      tetra(inistep)= sum(marg>0)/nfv; 
      teval(inistep) = sum(margval>0)/nfvval; 
      for icls=1:number_of_classes;
          cwtra(inistep,icls) = sum(marg(lbl==icls)>0)/sum(lbl==icls); 
          cwval(inistep,icls) = sum(margval(lblval==icls)>0)/sum(lblval==icls);
      end;    
      
      % roc with respect to class 1 versus all others only
     [~,~,auroc,thresholds] = compute_roc(lbl>1,score);
     auctra(inistep)=auroc; 
     [~,~,auroc,thresholds] = compute_roc(lblval>1,scoreval); 
     aucval(inistep)=auroc;
      
      
  
  end;
 
   [~,~,~,score]    = compute_costs(fvec,lbl,prototypes,plbl,omega,mu); 
  
% initial steps of Papari procedure complete, now remaining steps:
 
for jstep=(n_original+1):totalsteps;  
 % calculate mean positions over latest steps
 protmean = squeeze(mean(prototypes_original,1)); 
 omega_mean = squeeze(mean(omega_original,1));
 omega_mean=omega_mean/sqrt(sum(sum(omega_mean.^2))); 
 % note: normalization does not change cost function value
 %       but is done for consistency
 
 
 
 % compute cost functions for mean prototypes, mean matrix and both 
[costmp,~,~,score ] = compute_costs(fvec,lbl,protmean,plbl,omega, 0);
[costmm,~,~,score ] = compute_costs(fvec,lbl,prototypes,    plbl,omega_mean,mu); 
% [costm, ~,~,score ] = compute_costs(fvec,lbl,protmean,plbl,omega_mean,mu); 

% remember old positions for Papari procedure
ombefore=omega; 
protbefore=prototypes;
 
 % perform next step and compute costs etc.
[prototypes,omega]= do_batchstep (fvec,lbl,prototypes,plbl,omega,prototype_update_step_size,m_update_step_size,mu,mode);  
[costf,~,~,score] = compute_costs(fvec,lbl,prototypes,plbl,omega,mu); 

% by default, step sizes are increased in every step
 m_update_step_size=m_update_step_size*incfac; % (small) increase of step sizes
 prototype_update_step_size=prototype_update_step_size*incfac; % at each learning step to enforce oscillatory behavior 

% costfunction values to compare with for Papari procedure
% evaluated w.r.t. changing only matrix or prototype
[costfp,~,~,score] = compute_costs(fvec,lbl,prototypes,plbl,ombefore,0);
[costfm,~,~,score] = compute_costs(fvec,lbl,protbefore,plbl,omega,mu); 
   
% heuristic extension of Papari procedure
% treats matrix and prototype step sizes separately
 if (costmp <= costfp ); % decrease prototype step size and jump
                         % to mean prototypes
     prototype_update_step_size = prototype_update_step_size/decfac;
     prototypes = protmean;
 end; 
 if (costmm <= costfm ); % decrease matrix step size and jump
                         % to mean matrix
     m_update_step_size = m_update_step_size/decfac;
     omega = omega_mean;   
 end
                          
 % update the copies of the latest steps, shift stack 
 for iicop = 1:n_original-1;
   prototypes_original(iicop,:,:)=prototypes_original(iicop+1,:,:);
   omega_original(iicop,:,:)  =omega_original  (iicop+1,:,:);
 end;
 prototypes_original(n_original,:,:)=prototypes;  omega_original(n_original,:,:)=omega;
 
 % determine training and test set performances 
 % and calculate cost function without penalty terms
[costf,~,marg,score] = compute_costs(fvec,lbl,prototypes,plbl,omega,0); 
[costval,~,margval,scoreval]= ...
                     compute_costs(fvecval,lblval,prototypes,plbl,omega,0);

       
 cftra(jstep)= costf; cfval(jstep)= costval;
 tetra(jstep)= sum(marg>0)/nfv; teval(jstep) = sum(margval>0)/nfvval; 
 
 for icls=1:number_of_classes;
     cwtra(jstep,icls) = sum(marg(lbl==icls)>0)/sum(lbl==icls); 
     cwval(jstep,icls) = sum(margval(lblval==icls)>0)/sum(lblval==icls);
 end; 
 
 % roc with respect to class 1 versus all others only
 [~,~,auroc,thresholds] = compute_roc(lbl>1,score);
 auctra(jstep)=auroc; 
 [~,~,auroc,thresholds] = compute_roc(lblval>1,scoreval); 
 aucval(jstep)=auroc;
 
end; % end of totalsteps gradient steps (jrun)

 % output final lvq configuration after totalsteps steps: 
 w  =prototypes;   omega =omega;

 % end of batch gradient descent