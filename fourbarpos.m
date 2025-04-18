%%**********************************************************************
% 
%  fourbarpos.m
%  A function that solves a four bar position analysis problem. Only solves
%  for real roots. 
%
%   Input:  l_vec = vector of link lengths
%           th_1 = ground angle (degrees)
%           th_2 = input angle (degrees)
%           delta = crossed or open (open = -1, crossed = +1)
%
%   Output: th_vec = vector of angles including the input angles, th_1 and
%   th_2 (in degrees). Always between 0 and 360
% 
%  Meredith Symmank
% 
%%**********************************************************************

function [th_vec] = fourbarpos(l, th_1, th_2, delta)
   
   % Find my K constants from my link lengths
   K1 = l(1)/l(2);
   K2 = l(1)/l(4);
   K3 = (l(2)^2 - l(3)^2 + l(4)^2 + l(1)^2)/(2 * l(2) * l(4));
   K4 = l(1)/l(3);
   K5 = (l(4)^2 - l(1)^2 - l(2)^2 - l(3)^2)/(2 * l(2) * l(3));

   % Quadratic constants with the input angle theta 2
   A = -K1 - K2*cosd(th_2) + K3 + cosd(th_2);
   B = -2*sind(th_2);
   C = K1 - K2*cosd(th_2) + K3 - cosd(th_2);
   D = cosd(th_2) - K1 +K4*cosd(th_2) + K5;
   E = -2*sind(th_2);
   F = K1 + (K4 - 1)*cosd(th_2) + K5;

   % Check for imaginary roots. We won't handle those here.
   disc_4 = B^2-4*A*C;
   disc_3 = E^2-4*D*F;
   if disc_4 < 0 || disc_3 < 0
       error("This function does not handle imaginary roots")
   end
   
   if delta == -1 %open 
       th_4 = 2*atan2d((-B + delta*sqrt(B^2-4*A*C)),(2*A));
       th_3 = 2*atan2d((-E + delta*sqrt(E^2-4*D*F)),(2*D));

   elseif delta == 1 %crossed
       th_4 = 2*atan2d((-B - delta*sqrt(B^2-4*A*C)),(2*A));
       th_3 = 2*atan2d((-E - delta*sqrt(E^2-4*D*F)),(2*D));

   else 
        fprintf("Please enter the correct type of circuit: open: -1; crossed: +1");
   end 

   % Get it between 0 and 360
   th_4 = mod(th_4,360);

   % Get it between 0 and 360
   th_3 = mod(th_3,360);

   % Output vector
   th_vec = [th_1, th_2, th_3, th_4];
end