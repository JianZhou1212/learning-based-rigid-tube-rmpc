%*****************************************************************
%% primUB: compute a upper bound for the exact primal
%%         optimal value. 
%%
%%*****************************************************************

  function UB = primUB(blk,At,C,b,X,y,Z,mu) 

  if (nargin < 8); mu = 1.1; end
  AX = AXfun(blk,At,[],X);
  Rp = b-AX; 
  UB0 = blktrace(blk,C,X); 
  pert = mu*norm(Rp)*norm(y); 
  for p = 1:size(blk,1)
     pblk = blk(p,:);
     if strcmp(pblk{1},'s')
        eigXp = eig(full(X{p}));        
        maxeigZ = max(eig(full(Z{p}))); 
     elseif strcmp(pblk{1},'l')
        eigXp = full(X{p});        
        maxeigZ = max(Z{p}); 
     end      
     pert = pert + (mu*maxeigZ)*sum(max(0,-eigXp));     
  end
  UB = UB0 + pert; 
  fprintf('\n primal obj = %-10.9e',UB0); 
  fprintf('\n valid UB   = %-10.9e\n',UB);   
%%*****************************************************************