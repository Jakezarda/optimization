import math
import random
import time

N_darts = 1048576
N_trials = 1000
pis = []

t_1 = time.time()
for trial in range(0, N_trials):
    count = 0
    for i in range(0, N_darts):
        x = random.random()
        y = random.random()
        r = math.sqrt(x*x + y*y)
        if (r < 1 ):
            count += 1
    pi = 4.0*count/N_darts
    pis.append(pi)
t_2 = time.time()
print("Time to compute", N_trials, "estimates of pi:", t_2 - t_1, "s")

f = open('pythondumb.dat', 'w')
for pi in pis:
    f.write('%f' % pi)
f.close()
