function E_LJ =lennard_jones(x,sigma,epsilon)
  % E_LJ=lennard_jones(x,sigma,epsilon) 
  % returns the energy based upon a lennard jones potential, E_LJ
  % based upon the bond length, x, and two parameters,,sigma and epsilon
  E_LJ = 4*epsilon*((sigma./x).^12-(sigma./x).^6);
end

