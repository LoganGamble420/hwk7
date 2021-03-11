import numpy as np

 
def f(r,G=6.647e-11,M=1.989e30,m=5.972e24,R=1.496e11,w=1.991e-7): #when lflag=True we are in L2, when false we are in L1
    if r>R:
        return ((G*M)/r**2)+((G*m)/(R-r)**2)-(w**2)*r
    if r<R:
        return ((G*M)/r**2)-((G*m)/(R-r)**2)-(w**2)*r

def derivative(r,G=6.647e-11,M=1.989e30,m=5.972e24,R=1.496e11,w=1.991e-7): #when lflag=true we are in L2, when false l1
    if r>R:
        return ((-2*G*M)/r**3)+((2*G*m)/(R-r)**3)-(w**2)
    if r<R:
        return ((-2*G*M)/r**3)-((2*G*m)/(R-r)**3)-(w**2)

guess=1.0
answer=0

def newtons_method(x,func=None,df=None,max_iterations=None):
    '''Function to compute newtons method
        it takes in x as the first guess'''
    guess = x
    for i in range(max_iterations):
       fx = func(guess)
       fprime = df(guess)
       guess = guess - fx/fprime
    return guess

def secant_method(x, x2, func=None, deriv=None, max_iterations=100):
    '''Function  '''
    firstguess, secondguess = x, x2
    for i in range(max_iterations):
        fx1 = f(firstguess)
        fx2= f(secondguess)
        firstguess = secondguess - f(secondguess) *((secondguess-firstguess)/(f(secondguess)-f(firstguess)))

    return firstguess


rL1 = newtons_method(1e11,func=f,df=derivative,max_iterations=100) 
rL2= newtons_method(1e12,func=f,df=derivative,max_iterations=100) 
rSL1 = secant_method(1e11, 1.1e11,func=f,deriv=derivative,max_iterations=100)
rSL2 = secant_method(1.5e11, 2e11,func=f,deriv=derivative,max_iterations=100)


R = 1.496e11

answer= R - rL1
answer2= R + answer
answerS= R - rSL1
answerS2= R + answerS

print('Newtons Method: ')
print('\n L1 in km : ' , answer/1000)
print('\n L2 in km : ' , answer2/100000)
print('\nSecants Method: ')
print('\n L1 in km : ' , answerS/1000)
print('\n L2 in km : ' , answerS2/100000)






