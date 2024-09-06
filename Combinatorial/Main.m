clear
clc

for tour = 1:50
    [Result_BA1C(tour).Cost , Result_BA1C(tour).NFE, Result_BA1C(tour).OptCost] = Function_BA1_NN_comparativeLocal('Eil51', 51);
end

save Hamming_Eil_51_city_size_with_Tri_comp_Local Result_BA1C
