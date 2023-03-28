%%*****************************************************************
%% schurmat_lblk: compute A*D*A'
%%*****************************************************************
%% SDPT3: version 4.0
%% Copyright (c) 1997 by
%% Kim-Chuan Toh, Michael J. Todd, Reha H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*****************************************************************

   function [schur,UU,EE] = schurmat_lblk(blk,At,par,schur,UU,EE,p,dd)

   global idxdenAl
   
   iter = par.iter; 
   n = sum(blk{p,2});
   m = size(At{p},2); 

   if (iter==1) 
      idxdenAl{p} = checkdense(At{p}'); 
   end
   idxden = idxdenAl{p};  
   len = length(idxden);    
   ddsch = dd{p}; 
   removed_densecol = 0; 
   if (len > 0) && (len < 0.5*m && m > 500) %%2020-Apr-20
      Ad = At{p}(idxden,:)' *spdiags(sqrt(ddsch(idxden)),0,len,len); 
      UU = [UU, Ad];
      if isempty(EE)
         count = 0; 
      else
         count = max(max(EE(:,1)),max(EE(:,2))); 
      end
      tmp = count + [1:len]'; 
      EE = [EE; [tmp, tmp, -ones(len,1)] ]; 
      ddsch(idxden) = zeros(len,1); 
      removed_densecol = 1; 
   end
   if (len > 500) && ~removed_densecol
      Ad = full(At{p}(idxden,:)'*spdiags(sqrt(ddsch(idxden)),0,len,len));
      ddsch(idxden) = 0; %%must put here
      schurtmp = Ad*Ad' + At{p}'*spdiags(ddsch,0,n,n)*At{p}; 
   else
      schurtmp = At{p}' *spdiags(ddsch,0,n,n) *At{p}; 
   end
   schur = schur + schurtmp;
%%*******************************************************************
