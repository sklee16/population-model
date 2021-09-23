"""
Created on Tue Feb 26 10:02:00 2019
MT2507 Computing Project

Part 1: A Single Species Population Model: Biochemical Switch
Part 2: A Modified 2-Species Population Model

@author: sl288
"""
#Import Packages
from numpy import *
from pylab import *

"""Part 1 Question 1"""

#Even though the exponential can be written with **, for a faster program
#exponent is written as x * x instead of x ** 2, for example

#Biochemical switch model with a chemical substance s decaying at a linear rate of -rx 
def f(x):
    s = 0.01
    r = 0.4
    return s - (r * x) + ((x * x)/(1 + (x * x)))

#The derivative of the f(x)
def der_f(x):
    r = 0.4
    return -r + ((2*x)/ ((x*x) + 1))- ((2 * x * x * x)/(((x * x) + 1) * ((x * x ) + 1)))


def estimate_guess():
    x = arange(-0.5, 2.5, 0.001)
    plot (x, f(x))
    xlabel('x')
    ylabel('f(x)')
    title('f(x) values of model for a Biochemical Switch on a Limited Range')
    plt.grid(True)
    show()
    print("""
   The function has total of three positive roots. Specifically, the 
   roots lies between 0 and 2.5. Based on the graph above, 2 roots 
   lie between 0.0 and 0.5, and another root between 2.0 and 2.5. 
   Since the first root is very close to 0.0, first initial guess has 
   been chosen as 0.05. Secondly, another root exists close to 0.5,
   thus the initial guess has been decided as 0.45. Finally, the third
   root is approximately located at 2.10. Hence, The roots of the f(x) 
   are approximately at 0.05, 0.45, and 2.10.""")
    
#Newton-Rhapson Method
def newton_rhapson(initial_guess):
    
    print ('Newton-Rhapson Method')
    print ('Initial Point x0: %.2f'%(initial_guess))
    
    
    #accuracy of four decimal places
    error = 0.0001
    iteration = 0
    guess = initial_guess

    while 1:
        func_value = f(guess)
        der_value = der_f(guess)
        x_k1 = guess - (func_value/der_value)
        
        #iterate to find steady state is within the error indicated
        if abs(x_k1 - guess) < error:
            sol = x_k1
            print('n = %i, x%i = %.4f, x%i - x%i = %.4f'%(iteration, iteration, guess, 
                                            iteration, iteration - 1, x_k1 - guess))
            print('')
            print('Newton Method Result: %.4f'%(sol))
            break
        
        if iteration == 0:
            print('n = %i, x%i = %.4f, x%i = %.4f'%(iteration, iteration, guess, 
                                            iteration, x_k1))
        else:
            print('n = %i, x%i = %.4f, x%i - x%i = %.4f'%(iteration, iteration, guess, 
                                            iteration, iteration - 1, x_k1 - guess))
        iteration += 1
        guess = x_k1


"""Part 1 Question 2"""
def runge_kutta(t, x, h, n_steps):
    x_values = []
    t_values = []
    
    x_values.append(x)
    t_values.append(t)


    #Runge_Kutta Fourth Order Method Derivation
    for iteration in range(n_steps + 1):
        tn_1 = t + h
        
        #compute the slope at the beginning of the time step
        k_1 = h * f(x)
        
        #estimate of the slope at the midpoint of the previous slope
        k_2 = h * f(x + (0.5 * k_1))
        k_3 = h * f(x + (0.5 * k_2))
        k_4 = h * f(x +  k_3)
        
        #estimate the slope at the endpoint by averaging all the slopes
        xn_1 = x + ((1 / 6) * (k_1 + (2 * k_2) + (2 * k_3) + k_4))
        
        
        x_values.append(xn_1)
        t_values.append(tn_1)
        t = tn_1
        x = xn_1

    plot(t_values, x_values, 'k-')
    xlabel('t')
    ylabel('x(t)')
    title('Fourth Order Runge-Kutta Method with Initial Conditions of ' +
          'x(0) = 0.45 and x(0) = 0.5')
    

#Part 1  
def part1_graph():
    print('')
    print ('Part 1: A Single Spieces Population Model')
    print('')
    estimate_guess()
    print ('')
    
    #Part 1 Problem 1
    print('Part 1 Problem 1')
    print('')
    newton_rhapson(0.05)
    print('________________________________________________________')
    print ('')
    newton_rhapson(0.45) 
    print('________________________________________________________')
    print ('')
    newton_rhapson(2.10)
    print('________________________________________________________')
    print ('')
    print ('')
    
    #Part 1 Problem 2
    print('Part 1 Problem 2')
    print('')
    print('With step size of 0.01 for 500 iterations, you can see that the '+
          'solution with initial condition x(0) = 0.45 is decreasing, ' +
          'whereas x(0) =0.5. However, to clearly see the behavior of ' +
          'of the solutions, h = 0.05 and 250 total number of steps was chosen.')
    #Part a: x(0) = 0.45
    runge_kutta(0, 0.45, 0.05, 250)
    #Part b: x(0) = 0.5
    runge_kutta(0, 0.5, 0.05, 250)
    show()


"""Part 2"""
def comp_species_dx(x, y):
    "Computes the dy/dt model of y(t) species"
    #Given parameters
    a1 = 1.0
    b1 = 1.0
    c1 = 0.1
    d1 = 0.1
    
    f = x * (a1 - (b1 * x) - (c1 * exp(d1 * y)))
    
    return f

def comp_species_dy(x, y):
    "Computes the dy/dt model of y(t) species"
    #Given parameters
    a2 = 1.0
    b2 = 1.0
    c2 = 0.2
    d2 = 0.3
    f = y * (a2 - (b2 * y) - (c2 * exp(d2 * x)))
    
    return f
    
"""Part 2 Question 1"""
#Newton-Rhapson Method
        
def steady_state(x):
    """Computes the function values of set of equations 
    for the fourth steady state"""
    
    f = zeros(len(x))
    a1 = 1.0
    b1 = 1.0
    a2 = 1.0
    b2 = 1.0
    c1 = 0.1
    d1 = 0.1
    c2 = 0.2
    d2 = 0.3
    
    #Systems of equations that satisfies the fourth steady state
    f[0] = a1 - (b1 * x[0]) - (c1 * exp(d1 * x[0]))
    f[1] = a2 - (b2 * x[1]) - (c2 * exp(d2 * x[1]))
    return f

def jacobian(x):
    """Computes the Jacobian matrix for Newton Raphson 2D 
        scheme at each x values"""
    b1 = 1.0
    b2 = 1.0
    c1 = 0.1
    d1 = 0.1
    c2 = 0.2
    d2 = 0.3
    
    jacobian_mat = zeros((2, 2))
    
    #Note: f1 = a1 - (b1 * x[0]) - (c1 * exp(d1 * x[0]))
    #      f2 = a2 - (b2 * x[1]) - (c2 * exp(d2 * x[1]))
    #Following are derivatives with respect to a variable that is indicated
    #at the end of the function name (e.g. f1_x is a derivative of f1 with 
    #respect to x)
    f1_x = -b1
    f1_y = -c1 * d1 * exp(d1 * x[1])
    f2_x = -c2 * d2 * exp(d2 * x[0])
    f2_y = -b2
    
    jacobian_mat[0][0] = f1_x
    jacobian_mat[0][1] = f1_y
    jacobian_mat[1][0] = f2_x
    jacobian_mat[1][1] = f2_y
    return jacobian_mat


def newton_raphson_2d(guess):
    """Solves the system of equations for f(x) = 0 by using the 
        Newton-Raphson 2D method.
        Note: steady_state is a vector."""
        
    #accuracy of four decimal places
    eps = 0.0001
    iteration = 0
    
    while (1):
        func_value = steady_state(guess)
        x_k1 = guess - (np.dot(jacobian(guess), func_value))
        
        #Found x and y are within the error when it satisfies the
        #following if statement
        if abs(x_k1[0] - guess[0]) < eps and abs(x_k1[1] - guess[1]) < eps:
            return x_k1

        iteration += 1
        
        #No solution found by the Newton-Raphson Scheme
        if iteration >= 100:
            print('Too many iterations')
            
        guess = x_k1
        
def steady_states(print_value):
    """Calculates all four steady states for the modivied 2-specifies
    population model"""
    
    #Given parameters
    a1 = 1.0
    a2 = 1.0
    b1 = 1.0
    b2 = 1.0
    c1 = 0.1
    c2 = 0.2
    
    state_1 = [0, 0]
    state_2 = [0, (a2 - c2)/ b2]
    state_3 = [(a1 - c1)/ b1, 0]
    x = [1.0, 1.0]
    state_4 = newton_raphson_2d(x)

    if (print_value): #Print steady states values only when required
        print('The steady states are:')
        print(str(state_1) + ', ' + str(state_2) + ', '+  str(state_3) + ', ' + str(state_4))
        print('')
        
    plot(state_1[0],state_1[1], '*b')
    plot(state_2[0],state_2[1], '*g')
    plot(state_3[0],state_3[1], '*r')
    plot(state_4[0],state_4[1], '*c')
    
"""Part 2 Question 2"""        
def runge_kutta_multi(t, x, y, h, n_steps, color):
    """Solves the systems of ODES using the fourth order Runge-Kutta method"""
    
    x_values = []
    y_values = []
    x_values.append(x)
    y_values.append(y)



    for iteration in range(n_steps):
        #compute the slope at the beginning of the time step
        k_1 = h * comp_species_dx(x, y)
        l_1 = h * comp_species_dy(x, y)
        
        #estimate of the slope at the midpoint of the previous slope
        k_2 = h * comp_species_dx(x + (0.5 * k_1), y + (0.5 * l_1))
        l_2 = h * comp_species_dy(x + (0.5 * k_1), y + (0.5 * l_1))
        k_3 = h * comp_species_dx(x + (0.5 * k_2), y + (0.5 * l_2))
        l_3 = h * comp_species_dy(x + (0.5 * k_2), y + (0.5 * l_2))
        k_4 = h * comp_species_dx(x + k_3, y + l_3)
        l_4 = h * comp_species_dy(x + k_3, y + l_3)

        #estimate the slope at the endpoint by averaging all the slops
        xn_1 = x + ((1 / 6) * (k_1 + (2 * k_2) + (2 * k_3) + k_4))
        yn_1 = y + ((1 / 6) * (l_1 + (2 * l_2) + (2 * l_3) + l_4))
        
        x_values.append(xn_1)
        y_values.append(yn_1)
        x = xn_1
        y = yn_1

        plot(x_values, y_values, color)
        xlabel('x(t)')
        ylabel('y(t)')
        title('Modified 2-Species Population Method Solved ' 
              + 'Using Fourth Order Runge-Kutta Method')

def part2_graph():
    print('Part 2 Problem 1')    
    print('')
    steady_states(True)
    print('Part 2 Problem 2')    
    
    print('')
    print('About 3-4 trajectories are plotted for each of the phase space ' +
          'trajector to illustrate the solutions of the system of ODE. ' + 
          'For each of the solutions computed with Runge-Kutta method, ' +
          'step size of 0.01 was chosen to illustrate a precise solution, ' +
          'with a total number of steps ranging from 250 to 400 ' +
          'to clearly see the behvaior of solutions, especially with ' +
          ' saddle phase solutions. ')
    #trajectories for steady state at (0, 0)
    runge_kutta_multi(0, 0.2, 0.03, 0.01, 250, 'g-')
    runge_kutta_multi(0, -0.04, 0.1, 0.01, 250, 'g-')
    runge_kutta_multi(0, -0.03, -0.02, 0.01, 250, 'g-')
    runge_kutta_multi(0, 0.01, -0.05, 0.01, 250, 'g-')

    #trajectories for steady state at (0, a2 - c2 / b2)
    runge_kutta_multi(0, 0.01, 1.4, 0.01, 400, 'b-')
    runge_kutta_multi(0, -0.02, 1.8, 0.01, 300, 'b-')
    runge_kutta_multi(0, 0.01, 0.2, 0.01, 400, 'b-')

    #trajectories for steady state at (a1-c1 / b1, 0)
    runge_kutta_multi(0, 1.6, 1.4, 0.01, 250, 'r-')
    runge_kutta_multi(0, 1.7, 0.8, 0.01, 250, 'r-')    
    runge_kutta_multi(0, 0.6, 1.5, 0.01, 250, 'r-')

    #trajectories for steady state that satifies the set of equations solved 
    #by newton-raphson method
    runge_kutta_multi(0, 0.6, -0.02, 0.01, 400, 'c-')
    runge_kutta_multi(0, 1.1, -0.05, 0.01, 300, 'c-')
    runge_kutta_multi(0, 1.1, 0.15, 0.01, 400, 'c-')

    show()

#Display the results for part 1 and part 2
part1_graph()
part2_graph()