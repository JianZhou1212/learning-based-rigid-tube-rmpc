%%***********************************************************************
%% Given At*yinput that is approximately <= C
%% Find an y that is dual feasible: At*y <= C, and y is close to yinput
%% min { norm(y-yinput) | At*y <= C } which is equivalent to 
%% max s 
%% s.t. (0,At)*[s;y] + Z = C, Z >= 0
%%      diag([-1,ones(m,1)])*(s;y) + [t;v] = [0;yinput], [t;v] >= 0 (socp)
%%***********************************************************************
%% Copyright (c) 2018 by
%% Kim-Chuan Toh
%%***********************************************************************

   function [bblk,AAt,CC,bb] = findDualFeasiblePoint(blk,At,C,b,yinput)

   numblk = size(blk,1);   
   bblk = blk; 
   AAt = At; 
   CC = C; 
   for p=1:numblk
       pblk = blk(p,:); 
       n = sum(pblk{2});       
       if strcmp(pblk{1},'s')
          n2 = sum(pblk{2}.*(pblk{2}+1))/2; 
          ZZ = sparse(n2,1);
          AAt{p,1} = [ZZ,At{p}];
       elseif strcmp(pblk{1},'l');
          ZZ = sparse(n,1);
          AAt{p,1} = [ZZ,At{p}];
       end 
   end
   m = size(At{1},2);   
   options=1;
   if (options==1) %%minimize norm(y-yinput,inf) %%better than option=2
      AAt{numblk+1,1} = [-sparse(ones(m,1)), speye(m,m)];
      CC{numblk+1,1} = yinput;
      bblk{numblk+1,1} = 'l'; bblk{numblk+1,2} = m; 
      AAt{numblk+2,1} = [-sparse(ones(m,1)), -speye(m,m)];
      CC{numblk+2,1} = -yinput;
      bblk{numblk+2,1} = 'l'; bblk{numblk+2,2} = m;      
      bb = [-1; zeros(m,1)];        
   elseif (options==2) %%minimize norm(y-yinput,2)
      AAt{numblk+1,1} = spdiags([-1; ones(m,1)],0,m+1,m+1);  
      CC{numblk+1,1} = [0; yinput];
      bblk{numblk+1,1} = 'q'; bblk{numblk+1,2} = m+1; 
      bb = [-1; zeros(m,1)];       
   end
%%*****************************************************
