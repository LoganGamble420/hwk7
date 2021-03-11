import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time

def f(x,c):
    return 1-np.exp(-1*c*x)

def relaxation(start_guess=1, func_to_relax=None, tolerance=1e-6):
    """Function that computes the root via the fixed point
    (relaxation) method. inputs are a starting guess,
    a function to use, any function arguments, and a tolerance
    to exit the function when successive approximations are
    less than this value"""
    x1=start_guess
    x2=f(x1,c)
    while(np.abs(x2-x1) > tolerance):
        x1 = x2
        x2 = f(x2,c)
    return x2

if __name__ == "__main__":
    print("""a. Write a program to solve this equation for x using the relaxation method for the case c = 2. Calculate your solution to an accuracy of six decimal places.  (3 pts)""")
    c=2
    print("The solution for part a is:")
    print(relaxation(func_to_relax=f))
    print("""b) Modify your program to calculate the solution for values of R0 from 0 to 5 in steps of 0.01 and make a plot of P as a function of R0. Save this plot as a png file as an output to your program, with properly labeled axes and titles and such. You should see a clear transition from a regime in which P = 0 to a regime of nonzero P. This is another example of a phase transition. In physics this transition is known as the percolation transition; in epidemiology it is the epidemic threshold.  (7 pts)""")

    #loop over r0 values from 0 to 5
    r0_values = np.arange(0, 5.0, 0.01)
    solutions = []
    for r0 in r0_values:
        c=r0    
        answer = relaxation(func_to_relax=f)
        solutions.append(answer)

    #save the output data
    output_textfile = 'problem6_1_data.txt'
    np.savetxt(output_textfile,
               np.array(np.vstack((r0_values, solutions))).T,
               delimiter = ', ', header='R_0, Probability of epidemic',
               fmt = ('%.2f', '%.3e'))

data=pd.read_csv('problem6_1_data.txt',delimiter=',')
r0=data['# R_0']
prob=data[' Probability of epidemic']
plt.plot(r0,prob)
plt.xlabel('r0')
plt.ylabel('Probability of epidemic')
plt.show()
