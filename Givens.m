function [c,s,r] = Givens(x,y)
%
%  computes c,s,r such that
%
%  [c  s] [x] = [r]
%  [-s c] [y]   [0]
%
%  with c*c + s*s = 1;
%
   if (y == 0),
      c = 1; s = 0; r = x;
   else 
      if (abs(x) >= abs(y))
         t = y/x; r = sqrt(1 + t*t);
         c = 1/r;
         s = t*c;
         r = x*r;
      else 
         t = x/y; r = sqrt(1 + t*t);
         s = 1/r;
         c = t*s;
         r = y*r;
      end
   end 
      
