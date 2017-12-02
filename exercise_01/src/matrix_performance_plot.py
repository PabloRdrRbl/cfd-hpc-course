import numpy as np
import matplotlib.pyplot as plt

# matplotlib conf
plt.style.use('seaborn-paper')

# Data

n = [50, 200, 500, 1000, 2000]  # System size                                  

baseline_time = np.array([0.2213001E-04, 0.1197290E-03, 0.1239267E-01, 0.9722320E-01, 0.8599598E+00])
baseline_gflops = np.array([0.1129688E+02, 0.1670438E+02, 0.2017322E+02, 0.2057122E+02, 0.1860552E+02])

simple_time = np.array([0.1015782E-03, 0.9546518E-03, 0.1390355E+00, 0.2110359E+01, 0.2767759E+02])
simple_gflops = np.array([0.2461157E+01, 0.2095005E+01, 0.1798102E+01, 0.9477059E+00, 0.5780849E+00])

blocked_time = np.array([0.1008892E-03, 0.8807611E-03, 0.9812286E-01, 0.7834278E+00, 0.6509825E+01])
blocked_gflops = np.array([0.2477966E+01, 0.2270763E+01, 0.2547826E+01, 0.2552884E+01, 0.2457824E+01])

# Figure: Size vs. Time
plt.figure(1)

plt.semilogy(n, baseline_time, '--o', label='Baseline')
plt.semilogy(n, simple_time, '--o', label='Simple')
plt.semilogy(n, blocked_time, '--o', label='Blocked ($nn=50$)')

plt.xlabel('System size $n$')
plt.ylabel('Elapsed time $[s]$')

plt.title('System Size vs. Elapsed Time for different matmul implementations')

plt.legend()
plt.grid()

plt.savefig('elapsed_time.png')

# Figure: Size vs. GFLOPS
plt.figure(2)

plt.semilogy(n, baseline_gflops, '--o', label='Baseline')
plt.semilogy(n, simple_gflops, '--o', label='Simple')
plt.semilogy(n, blocked_gflops, '--o', label='Blocked ($nn=50$)')

plt.xlabel('System size $n$')
plt.ylabel('GFLOPS')

plt.title('System Size vs. Performance for different matmul implementations')

plt.legend()
plt.grid()

plt.savefig('performance.png')
