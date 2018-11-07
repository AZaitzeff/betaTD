function [intervals,count,sizes] = decomposez(binseq)
n=size(binseq,2);
M=800;
intervals = zeros(10,M);
sizes=zeros(10,1);
count = 0; % Number of patches.

pos = 1;
flag = 0; % 1 if inside a patch of 1's.
patchsize = 0; % Size of the patch of 1's.
patch = zeros(1,n); % Used to record the locations in the patch.
while pos <= n

  if binseq(pos) == 1
    if flag == 0 % This is a new patch of 1's.
      flag = 1;
      count = count + 1;
      if count>10
          count=10;
          'warning count'
          break
      end
      patchsize = 1;
      patch(patchsize) = pos; % Record the position.
    else % We're in an existing patch of 1's.
      patchsize = patchsize + 1;
      patch(patchsize) = pos; % Recod the position.
    end
  end

  if ( binseq(pos) == 0 ) || ( pos == n )
    if flag == 1 % If we have just exited a patch of 1's...
      flag = 0;
      if patchsize>M
          'Warning patch to big'
      end
      pathlength=min(patchsize,M);
      intervals(count,1:pathlength)=patch(1:pathlength);
      sizes(count)=pathlength;
      patch(1:patchsize) = 0;
    end
  end

  pos = pos + 1;

end