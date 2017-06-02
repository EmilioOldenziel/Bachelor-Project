function [costf,crout,marg,score]=compute_costs(fvec,lbl,prot,plbl,omat,mu,mode)

% calculates gmlvq cost function, labels and margins for a set of 
% labelled feature vectors, given a particular lvq system

% arguments: 
% fvec             nvf feature vectors of dim. ndim  fvec(1:nfv,1:ndim);
% lbl              data labels  lbl(1:nfv);
% prot,plbl        prototypes, prototype labels
% omat             matrix omega

% output:
% costf  glvq-costs per example 
% marg   margins of classifying the examples
% crout  crisp classifier outputs 
% score  distance based 'scores' 

 nfv= size(fvec,1); ndim = size(fvec,2);  % # and dim. of feature vectors
 np = length(plbl); 
 
  costf =0; marg=zeros(1,nfv); score=marg; crout=zeros(1,nfv);
  if(mode==4)
    for iom=1:np
      omat(:,:,iom) = omat(:,:,iom)/sqrt(sum(sum(omat(:,:,iom).*omat(:,:,iom)))); % normalized omat
    end
  else
    omat = omat/sqrt(sum(sum(omat.*omat)));
  end;
    for iii=1:nfv;                        % loop through examples
       fvc = fvec(iii,:); lbc=lbl(iii);
       distl = nan(np,1);
          for jk=1:np;
              if(mode==4)
                distl(jk) = norm(omat(:,:,jk)*(fvc-prot(jk,:))')^2;
              else
                distl(jk) = norm(omat*(fvc-prot(jk,:))')^2;
              end;
          end;
       % find the two winning prototypes for example iii
       correct   = find (plbl == lbc); 
       incorrect = find (plbl ~= lbc);
       [dJJ, JJJ] = min (distl (correct));
       [dKK, KKK] = min (distl (incorrect));
       JJ=correct(JJJ); KK=incorrect(KKK);   % winner indices 

       costf = costf + (dJJ-dKK)/(dJJ+dKK)/nfv;
    
       marg(iii) = (dJJ-dKK)/(dJJ+dKK);  % gmlvq margin of example iii
       
       
       % un-normalized difference of distances
       if (lbc==1) 
          % score(iii)= 1./(1+exp((dKK-dJJ)/2));  % "the larger the better"
            score(iii)= dKK-dJJ; 
          % distdiff(iii)=dKK-dJJ;
          % score (iii) = 0.5* (1+marg(iii)); 
       else 
          % score(iii)= 1./(1+exp((dJJ-dKK)/2));  % "the larger the worse" 
            score(iii)= dJJ-dKK; 
          % distdiff(iii)=dJJ-dKK;
          % score (iii) = 0.5* (1-marg(iii)); 
       end
           
       crout(iii)= plbl(JJ)*(dJJ<=dKK) + plbl(KK)*(dJJ>dKK); 
                   % the class label according to nearest prototype                  
     end;
     
     % add penalty term 
     if (mu>0);
        if(mode==4)
          costf=costf-mu/2*log(det(omat(:,:,JJ)*omat(:,:,JJ)'))/nfv;
        else
          costf=costf-mu/2*log(det(omat*omat'))/nfv; 
        end;
     end; 
     
end

