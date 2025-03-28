import numpy as np

g2_exact = np.loadtxt("../build/g2.dat")
g2_inter = np.loadtxt("../build/g2_interp.dat")

if g2_exact.shape[0] != g2_inter.shape[0]:
    print(
        "Different values of xBj used in the two computations, automatic detection of similar xBj not yet implemented"
    )
    exit(-1)

if np.amax(np.abs(g2_exact[:, 0] - g2_inter[:, 0])) > 1.0e-15:
    print(
        "Different values of xBj used in the two computations, automatic detection of similar xBj not yet implemented"
    )
    exit(-1)


m = (
    np.abs(g2_exact[0, 1] - g2_inter[0, 1])
    if np.abs(g2_exact[0, 1]) < 1.0e-10
    else np.abs(1.0 - g2_inter[0, 1] / g2_exact[0, 1])
)
i_m = 0

m_abs = np.abs(g2_exact[0, 1] - g2_inter[0, 1])
i_m_abs = 0

for i in range(1, g2_exact.shape[0]):
    m1 = (
        np.abs(g2_exact[i, 1] - g2_inter[i, 1])
        if np.abs(g2_exact[i, 1]) < 1.0e-10
        else np.abs(1.0 - g2_inter[i, 1] / g2_exact[i, 1])
    )
    m1_abs = np.abs(g2_exact[i, 1] - g2_inter[i, 1])
    if m1 > m:
        m = m1
        i_m = i
    if m1_abs > m_abs:
        m_abs = m1_abs
        i_m_abs = i

print(
    "Max relative difference: ",
    m,
    " for [x, g2]: [",
    g2_exact[i_m, 0],
    ", ",
    g2_exact[i_m, 1],
    "]",
)
print(
    "Max relative difference: ",
    m,
    " for [x, g2]: [",
    g2_exact[i_m, 0],
    ", ",
    g2_inter[i_m, 1],
    "]",
)

print(
    "Max absolute difference: ",
    m_abs,
    " for [x, g2]: [",
    g2_exact[i_m_abs, 0],
    ", ",
    g2_exact[i_m_abs, 1],
    "]",
)

print(
    "Max absolute difference: ",
    m_abs,
    " for [x, g2]: [",
    g2_exact[i_m_abs, 0],
    ", ",
    g2_inter[i_m_abs, 1],
    "]",
)
