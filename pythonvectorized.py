import numpy as np
import time

N_darts = 1024*1024
N_trials = 1000
pis = []

t_1 = time.time()

for trial in range(0, N_trials):
    x = np.random.rand(N_darts)
    y = np.random.rand(N_darts)
    r = np.sqrt(x*x + y*y)
    count = r[ np.where(r <= 1) ]
    pi = 4.0*len(count)/N_darts
    pis.append(pi)
t_2 = time.time()
print("Time to compute", N_trials, "estimates of pi:", t_2 - t_1, "s")
print("The average value of pi was found to be: ", np.average(pis))
print("The standard deviation of pi was: ", np.std(pis))

f = open('pythonvectorized.dat', 'w')
for pi in pis:
    f.write('%f' % pi)
f.close()
    
    
    
 
