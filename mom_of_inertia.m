function [mom_of_inertia] = mom_of_inertia(a,b,h, thickness, dens)
mom_of_inertia = thickness*dens* (b^3*h-b^2*h*a+b*h*a^2+b*h^3)/36;
end

