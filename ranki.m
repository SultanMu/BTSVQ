function [out] = ranki (in)
N = length (in);
in2 = [1:1:N];
in = [in2' in'];
%out  = zeros (1,N);
j = 1;
for i = 1:N
   if in(i,2) == 1
      out([j]) = in(i,1);
      j = j+1;
   else
      xx= 0;
   end
end

      