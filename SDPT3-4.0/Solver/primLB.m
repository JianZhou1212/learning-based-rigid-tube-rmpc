%%*****************************************************************
%% primLB: compute a lower bound for the exact primal
%%         optimal value. 
%%
%%*****************************************************************

  function LB = primLB(blk,At,C,b,X,y,Z,mu) 

  if (nargin < 7); mu = 1.1; end
  Aty = Atyfun(blk,At,[],[],y);
  Znew = ops(C,'-',Aty); 

  LB0 = b'*y; 
  pert = 0; 
  for p = 1:size(blk,1)
     pblk = blk(p,:);
     if strcmp(pblk{1},'s')        
        eigtmp = eig(Znew{p}); 
        idx = find(eigtmp < 0); 
        maxeigX = max(eig(full(X{p}))); 
     elseif strcmp(pblk{1},'l')
        eigtmp = Znew{p};
        idx = find(eigtmp < 0); 
        maxeigX = max(X{p}); 
     end      
     numneg = length(idx); 
     if (numneg) 
        mineig = min(eigtmp(idx)); 
        fprintf('\n numneg = %3.0d,  mineigZnew = %- 3.2e',numneg,mineig);
        pert = pert + (mu*maxeigX)*sum(eigtmp(idx)); 
     end
  end
  LB = LB0 + pert; 
  fprintf('\n dual obj = %-10.9e',LB0);
  fprintf('\n valid LB = %-10.9e\n',LB);  
%%*****************************************************************
