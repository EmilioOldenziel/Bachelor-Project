
% visualize the trained LVQ system in terms of prototypes
% and relevance matrix

function display_gmlvq(prot,lambda,plbl,ndim) 

% input
% prot  : array of prototype vectors
% lambda: relevance matrix (ndim x ndim)
% plbl  : prototype labels
% ndim  : dimension of prototype vectors 

 nprots=length(plbl);          % number of prototypes
 ncls  =length(unique(plbl));  % number of classes
 
 % define colors for up to seven classes
 symbstrings = ['b';'g';'r';'c';'m';'y';'k']; 

% geometry of figure showing prototypes
if nprots<=6; rows=3;  else rows=4;   end; 

col=4;  % number of columns
posdiag=col;              % position of diagonal element bar plot
poslam =rows*col;         % position of off-diagonal element matrix plot

posiw=[];                 % positions of prototype bar plots
for i=1:rows;
    posiw=[posiw,(1:col-1)+(i-1)*(col)];
end;
   
   for iw=1:nprots        
     % display prototypes as bar plots
     subplot(rows,col,((iw-1)*4)+1); 
       hold on;
       bar(prot(iw,:),'FaceColor',symbstrings(plbl(iw))); 
       title(['prot. ',num2str(iw), ', class ',num2str(plbl(iw))], ... 
          'FontName','LucidaSans', 'FontWeight','bold'); 
           axis([0.3 ndim+0.7 -max(max(abs(prot))) +max(max(abs(prot)))]); 
       grid on; 
       box on; 
   
  
   
    % display diagonal matrix elements as bar plot
    subplot(rows,col,((iw-1)*4)+2); 
      hold on;
      bar(sort(eig(lambda(:,:,iw)),'descend')); 
      title(strcat('eigenvalues of rel.mat. ','prot. ',num2str(iw)), ... 
          'FontName','LucidaSans', 'FontWeight','bold'); 
      xlabel('feature number'); 
      grid on;  axis 'auto y'; box on;
      axis([0.3 ndim+0.7 0 (0.01+max(diag(lambda(:,:,iw))))]); 
      hold off;
    
    subplot(rows,col,((iw-1)*4)+3); 
      hold on;
      bar(diag(lambda(:,:,iw))); 
      title(strcat('rel. matrix, diag.','prot. ',num2str(iw)), ... 
          'FontName','LucidaSans', 'FontWeight','bold'); 
      xlabel('feature number'); 
      grid on;  axis 'auto y'; box on;
      axis([0.3 ndim+0.7 0 (0.01+max(diag(lambda(:,:,iw))))]); 
      hold off;
      
      
    % display off-diagonal matrix elements as matrix
    subplot(rows,col,((iw-1)*4)+4); 
        lambdaoff = lambda(:,:,iw).*(1-eye(ndim));     % zero diagonal 
        imagesc(lambdaoff); colormap(summer); 
        axis square; 
        xlabel('off-diag. el.', ... 
          'FontName','LucidaSans', 'FontWeight','bold'); 
        if(col==2); 
          colorbar('Location','EastOutside');
        else
          colorbar('Location','SouthOutside');
        end;
        hold off;

   end; 
       
     

    
    
    
    
  