function A_RemovePiston = Remove_Piston(A)
A_RemovePiston = A-mean(A(~isnan(A)));
end