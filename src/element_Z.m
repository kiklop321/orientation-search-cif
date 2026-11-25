function Z = element_Z(symbol)
%ELEMENT_Z  Atomic number of a chemical element ([] if unknown).
%
%   Z = element_Z('Ti')   ->  22

persistent map
if isempty(map)
    s = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar', ...
         'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', ...
         'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe', ...
         'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', ...
         'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', ...
         'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf'};
    map = containers.Map(s, num2cell(1:numel(s)));
end
if isKey(map, symbol)
    Z = map(symbol);
else
    Z = [];
end
end
