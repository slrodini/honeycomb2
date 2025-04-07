import numpy as np

data = np.loadtxt(
    "/mnt/d/Dropbox/ThinkPad/Documents/Code/honeycomb/dat/results/full_kernels_1_to_10000/result_N120_M25_E_up_nf5.dat"
)
result = "#include <honeycomb2/honeycomb2.hpp>\n"
result += "std::vector<Honeycomb::RnC::Triplet> honeycomb_points = {\n"

for i in range(data.shape[0]):
    result += "{"
    result += f"{data[i][2]:.16e}, {data[i][3]:.16e}, {data[i][4]:.16e}"
    result += "}, \n"

result += "};\n"
print(result)
