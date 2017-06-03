function [f1_micro, f1_macro]=compute_f_measure(fvec,lbl,prot,plbl,omat,mu,mode)

    % calculates LVQ f-measure, recoil and precision

    % arguments: 
    % fvec             nvf feature vectors of dim. ndim  fvec(1:nfv,1:ndim);
    % lbl              data labels  lbl(1:nfv);
    % prot,plbl        prototypes, prototype labels
    % omat             matrix omega

    % output:
    % score:  overall f1 score
    % recoil: amount of what is classified
    % precision: amount of what is classified correctly

    nfv= size(fvec,1); ndim = size(fvec,2);  % # and dim. of feature vectors
    n_classes = length(unique(lbl));
    np = length(plbl); 

    crout=zeros(1,nfv);

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

        crout(iii)= plbl(JJ)*(dJJ<=dKK) + plbl(KK)*(dJJ>dKK); 
        % the class label according to nearest prototype         
    end;

    true_positives = zeros(1,nfv);
    false_positives = zeros(1,nfv);
    false_negatives = zeros(1,nfv);

    for iii=1:n_classes
        correct_fvecs = (lbl'==iii);
        true_positives(iii) = sum((crout == iii)==correct_fvecs);
        false_positives(iii) = sum((crout == iii)~=correct_fvecs);
        false_negatives(iii) = sum((crout ~= iii)==correct_fvecs);
    end

    precision_micro = sum(true_positives) / sum(true_positives+false_positives);
    recall_micro = sum(true_positives) / sum(true_positives+false_negatives);

    precision_macro = sum(true_positives / (true_positives+false_positives)) / n_classes;
    recall_macro = sum(true_positives / (true_positives+false_negatives)) / n_classes;

    f1_micro = 2*(precision_micro*recall_micro) / (precision_micro+recall_micro);
    f1_macro = 2*(precision_macro*recall_macro) / (precision_macro+recall_macro);
end

