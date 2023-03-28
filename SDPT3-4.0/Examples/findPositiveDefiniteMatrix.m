%%***********************************************************************
%% Test whether a positive definite matrix exists in the
%% range space of At.
%% (topt,yopt) = argmax{ t | -At*y + t*I <= 0, At*y <= alpha*I } 
%% where alpha is a given positive constant. 
%% If the optimal value topt is positive, then there exist a 
%% positive definite matrix At*yopt in the range space of At. 
%%***********************************************************************
%% Copyright (c) 2018 by
%% Kim-Chuan Toh
%%***********************************************************************

   function [bblk,AAt,C,b] = findPositiveDefiniteMatrix(blk,At)

   numblk = size(blk,1);   
   bblk = blk; 
   bblk(numblk+[1:numblk],:) = blk;
   AAt = cell(2*numblk,1); C = cell(2*numblk,1);
   normA = 0; 
   for p=1:numblk
       normA = normA+norm(At{p},'fro'); 
   end
   const = normA; 
   for p=1:numblk
       pblk = blk(p,:); 
       n = sum(pblk{2});       
       if strcmp(pblk{1},'s')
          n2 = sum(pblk{2}.*(pblk{2}+1))/2; 
          II = svec(pblk,speye(n,n),1);
          ZZ = sparse(n2,1);
          C{p,1} = sparse(n,n);
          C{numblk+p,1} = const*speye(n,n);
       elseif strcmp(pblk{1},'l');
          II = ones(n,1);
          ZZ = sparse(n,1);
          C{p,1} = ZZ;
          C{numblk+p,1} = const*II;   
       end
       AAt{p,1} = [-At{p},II]; 
       AAt{numblk+p,1} = [At{p},ZZ];     
   end
   m = size(At{1},2);
   b = [zeros(m,1); 1]; 
%%*****************************************************
