# -*- coding: utf-8 -*-
import numpy
import scipy
import matplotlib.pyplot as pyplot
import scipy.integrate
import random


''' a_0, l_0 hjælper med at lave enhederne om til noget mere fornuftigt og putte det i naturlige enheder sådan at 1/4pie_0 = 1 og q = 1 og derfor gøre det mere overskueligt'''
N = 7
spatial_points = 300
T = 2*numpy.pi #4/(2.4189*10**(-17)) #*1.6*10**(-11)/(6.58*10**(-34))
time_step = float(T)/float(N)
m = 1 #0.5*10**6
q = 1 #0.0854 1.6*10**(-19)
a_0 = 0.529177*10**(-10)
e_0 = 8.854*10**(-12)
r_0 = 1.25*10**(-15)
a_0n = 2.68*10**(-4)
l_0 = 1 #1.932705*10**(-7)
start_point = -3
end_point = 3
met_step = float(20000000)
number_times = float(1)
removed_iterations = 5
screening_parameter = 3.36
movement = 0.5

def f(grid,new,old):
    if ((abs(grid[old])+abs(grid[new]))/2.) > 5:
        V = 1000000
    else:
        V = -(2*q**2/(((abs(grid[old]))+(abs(grid[new])))/2.))*(1-numpy.exp(-screening_parameter*(((abs(grid[old]))+abs(grid[new])))/2.))
    return (0.5*m*(((abs(grid[new])-abs(grid[old])))/(time_step))**2 + V)#*((grid[new]+grid[old])/2.)
    #return V

def inter(grid1,grid2,new,old):
    r1 = (abs(grid1[new]) + abs(grid1[old]))/2.
    r2 = (abs(grid2[new]) + abs(grid2[old]))/2.
    g = abs(r2+r1)
    if g == 0:
        V = 10**(25)
    else:
        V = 1/(g)*(1-numpy.exp(-screening_parameter*g))
    return V

def E(N,spatial_points,T,time_step,m):
    '''Creating the sum over all possible path, so hence find the energy of all existing paths'''
    
    x1 = numpy.linspace(start_point,end_point,spatial_points) #There are seven dimensions and hence need seven possible set of lokations
    '''
    x2 = numpy.linspace(start_point,end_point,spatial_points)
    x3 = numpy.linspace(start_point,end_point,spatial_points)
    x4 = numpy.linspace(start_point,end_point,spatial_points)
    x5 = numpy.linspace(start_point,end_point,spatial_points)
    x6 = numpy.linspace(start_point,end_point,spatial_points)
    x7 = numpy.linspace(start_point,end_point,spatial_points)
    '''
    #grid = numpy.meshgrid(x1,x2,x3,x4,x5,x6,x7) #The grid of all possible paths in space
    #Epsilon = numpy.zeros((spatial_points,spatial_points,spatial_points,spatial_points,spatial_points,spatial_points,spatial_points))
    #Epsilon = numpy.zeros((spatial_points,met_step))
    point = numpy.random.randint(0,spatial_points,N)
    point1 = []
    point2 = []
    point11 = []
    point22 = []
    #point3 = []
    yes = []
    no = []
    chance = []
    nochance = []
    p = 0
    accept = []
    for i in range(N):
        #point1.append(x1[point[i]])
        #point2.append(x1[point[i]])
        #point3.append(x1[point[i]])
        point1.append(r_0)
        point2.append(r_0)
        point11.append(r_0)
        point22.append(r_0)
        #point3.append(0)
    #point2 = point1
    Epsilon = []
    #n = 1 #The dimension to start from, the 0 is excluded, due to the way of integration
    for j in range(int(met_step)):
       Epsilon_p = []
       for h in range(N):
           #if h < removed_iterations:
           #    displacement = (2*random.random()-1)/100.
           #else:
           #    displacement = (2*random.random()-1)/400.
           displacement1 = (2*random.random()-1)*movement
           displacement2 = (2*random.random()-1)*movement
           #print displacement
           #point3 = point1
           #point2 = point1
           #print point1[h],point2[h]
           point2[h] = point2[h] + displacement1
           point22[h] = point22[h] + displacement2
           #print point1[h],point2[h]
           #print 0
           #print point1
           #print point2
           #print point
           #print point2
           '''i1 = float(point[0])
           i2 = float(point[1])
           i3 = point[2]
           i4 = point[3]
           i5 = point[4]
           i6 = point[5]
           i7 = point[6]
           j1 = point2[0]
           j2 = point2[1]
           j3 = point2[2]
           j4 = point2[3]
           j5 = point2[4]
           j6 = point2[5]
           j7 = point2[6]'''
           energy1 = []
           energy2 = []
           #for i in range(0):
           '''if i == 0:
                   energy1.append(f(point3,0,5))# + 1/(4*numpy.pi*x1[point[i+1]]+0.000000000000001)
                   energy2.append(f(point2,0,5))# + 1/(4*numpy.pi*x1[point2[i+1]+0.000000000000001)
           '''
               #else:
               #energy1.append(f(point1,h+i,h+i-1))# + 1/(4*numpy.pi*x1[point[i+1]]+0.000000000000001)
               #energy2.append(f(point2,h+i,h+i-1))# + 1/(4*numpy.pi*x1[point2[i+1]+0.000000000000001)
           if h == 0:
               E1 = f(point1,h,6) + f(point1,h+1,h) + f(point11,h,6) + f(point11,h+1,h) + inter(point1,point11,h,6) + inter(point1,point11,h+1,h)
               E2 = f(point2,h,6) + f(point2,h+1,h) + f(point22,h,6) + f(point22,h+1,h) + inter(point2,point22,h,6) + inter(point2,point22,h+1,h)
           else:
               if h == 6:
                   E1 = f(point1,h,h-1) + f(point11,h,h-1) + inter(point1,point11,h,h-1)
                   E2 = f(point2,h,h-1) + f(point22,h,h-1) + inter(point2,point22,h,h-1)
               else:
                   E1 = f(point1,h,h-1) + f(point1,h+1,h) + f(point11,h,h-1) + f(point11,h+1,h) + inter(point1,point11,h,h-1) + inter(point1,point11,h+1,h)
                   E2 = f(point2,h,h-1) + f(point2,h+1,h) + f(point22,h,h-1) + f(point22,h+1,h) + inter(point2,point22,h,h-1) + inter(point2,point22,h+1,h)
           #print h,h-1
           #print point1
           #print point2
           point = point1 + point11
           point_n = point2 + point22
           #E1,E2 = sum(energy1),sum(energy2)
           deltaE = E2 - E1
           #diff = E1/E2
           #print deltaE
           if deltaE < 0:
               yes.append(1)
               nochance.append(1)
               if j < removed_iterations:
                   point1 = list(point2)
                   point11 = list(point22)
               else:
                   #Epsilon.append(point2[h])
                   Epsilon_p.append(point_n[h])
                   point1 = list(point2)
                   point11 = list(point22)
           else:
               if numpy.exp(-deltaE) > random.random():
                   yes.append(1)
                   chance.append(1)
                   accept.append(numpy.exp(-deltaE))
                   if j < removed_iterations:
                       point1 = list(point2)
                       point11 = list(point22)
                   else:
                       #Epsilon.append(point2[h])
                       Epsilon_p.append(point_n[h])
                       point1 = list(point2)
                       point11 = list(point22)
               else:
                   no.append(1)
                   #print no
                   if j < removed_iterations:
                       point2 = list(point1)
                       point22 = list(point11)
                   else:
                       #Epsilon.append(point1[h])
                       Epsilon_p.append(point_n[h])
                       point2 = list(point1)
                       point22 = list(point11) 
           #print Epsilon[i1,i2,i3,i4,i5,i6,i7]
       if j > removed_iterations:
            if p == 300: #int((spatial_points/6.)**2):
                #print int((spatial_points/6.)**2)
                Epsilon = Epsilon + Epsilon_p
                p = 0
                #print Epsilon_p
                #print Epsilon_p
            else:
                p = p + 1
            #print p
           #print deltaE
    print "accept",sum(yes)/float((sum(no)+sum(yes)))*100,"%"
    print "chance",sum(chance)/(float((sum(no)+sum(yes))))*100,"%"
    print "nochance",sum(nochance)/float(sum(yes)+sum(no))*100,"%"
    #print deltaE
    
    return Epsilon


def G():
    '''Creating unnormalized integral'''
    #G = numpy.zeros((number_times,spatial_points))
    G = []
    for k in range(int(number_times)):
        print "progress",k/number_times*100, "%"
        epsilon = E(N,spatial_points,T,time_step,m)
        G = G + epsilon
        '''epsilon_masked = numpy.ma.masked_equal(epsilon,0)
        #integrand = numpy.exp(-time_step*epsilon_masked)
        for i in range(spatial_points):
            array1 = epsilon_masked[i,:]
            array_no0 = array1.compressed()
            z = sum(array_no0)
            size = numpy.shape(array_no0)
            G[k,i] = size[0] #*(met_step/float(size[0]))'''
    #G_summed = G.sum(0)
    #print G
    return G

#E = E(N,spatial_points,T,time_step,m)
#accept = E(N,spatial_points,T,time_step,m)[1]
#G,accept = E
#G = G(N,spatial_points,T,time_step,m)[0]
#accept = F(N,spatial_points,T,time_step,m)[1]
#print G
def z(x,G):
    '''Normalize my integral'''
    #x = numpy.linspace(start_point,end_point,spatial_points)
    z = []
    for i in range(int(numpy.shape(G)[0])-1):
        '''Trapezian rule'''
        A = (G[i] + G[i+1])*0.5*(x[i+1]-x[i])
        z.append(A)
        #print A
    s = sum(z)
    return s
#z = z(N,spatial_points,T,time_step,m)

#Feynmann = G/z

#print Feynmann

def totprob(N,spatial_points,T,time_step,m):
    '''Finding total probability of particle being somewhere (is a safetycheck)'''
    x = numpy.linspace(start_point,end_point,spatial_points)
    totprob = 0
    for i in range(spatial_points-1):
        A = (Feynmann[i] + Feynmann[i+1])*0.5*(x[i+1]-x[i])
        totprob = totprob + A
    return totprob

position = numpy.linspace(start_point,end_point,spatial_points) #The x-axis for the numerical solution
position2 = numpy.linspace(-5,5,500) #The x-axis for the analytical solution
'''
oscilator = []
Creating the analytical solution
for i in range(10000):
    oscilator.append((numpy.exp(-0.5*(position2[i])**2)/numpy.pi**0.25)**2)
'''

radial1 = []
nonradial = []
for i in range(500):
    radial1.append((4/1.30752)*numpy.pi*numpy.exp(-1.6875*abs(position2[i])*2)*position2[i]**2)
    nonradial.append(numpy.exp(-1.6875*abs(position2[i])*2)/0.592593)

#oscilatorplot = oscilator/sum(oscilator) 

#print scipy.integrate.quad(lambda x: (numpy.exp(-0.5*x**2)/numpy.pi**0.25)**2,-numpy.inf,numpy.inf)[0]

#totalprob = scipy.integrate.quad(lambda x: (numpy.exp(-0.5*x**2)/numpy.pi**0.25)**2,-numpy.inf,numpy.inf)[0],totprob(N,spatial_points,V,T,time_step,m)
#print 'Probtheory=,Probmine=',totalprob

#print numpy.histogram(G,bins=10,normed=True)[1]
#print numpy.histogram(G,bins=10,normed=True)[0]
#print numpy.histogram(G,bins=10,normed=True)

def radial():
    G_1 = G()
    hi = numpy.histogram(G_1,bins=spatial_points,normed=True)
    x_a_0 = hi[1]
    y_a_0 = hi[0]
    x_final = []
    y_final = []
    for i in range(spatial_points):
        x_a = (x_a_0[i]+x_a_0[i+1])/2.
        y_final.append(y_a_0[i]*4*numpy.pi*x_a**2)
        x_final.append(x_a)
    return hi, x_final,y_final,G_1

#x_final,y_final = radial()

def bohr():
    hi1, x_final,y_final,G_2 = radial()
    x_height = [0] + [0] + x_final + [0] + [0]
    y_height = [0] + [0] + y_final + [0] + [0]
    x_bohr = 0
    y_bohr = 0
    x_h1 = []
    y_h1 = []
    for i in range(spatial_points):
        y_h = (y_height[i] + y_height[i+1] + y_height[i+2] + y_height[i+3] + y_height[i+4])/5.
        x_h = x_height[i+2]
        x_h1.append(x_h)
        y_h1.append(y_h)
        #hi1[i] = hi1[i]*0.5 + hi1[i+1]*0.5
        if y_h > y_bohr:
            y_bohr = y_h
            x_bohr = x_h
    return abs(x_bohr),x_h1,y_h1,G_2,y_final

x_bohr,x_final,y_final,G_3,y_f = bohr()

y_fsum = z(x_final,y_final)
ysum = z(x_final,y_f)
y_final /= y_fsum
y_f /= y_fsum

cc = z(position2,radial1)

print 'radial',cc

'''

bohr_list = []
screening = []
for i in range(30):
    bohr_radii = bohr()[0]
    print "bohr, screening=", bohr_radii,screening_parameter
    bohr_list.append(bohr_radii)
    screening.append(screening_parameter)
    screening_parameter += 0.2

print 'Bohr=',bohr_list
print 'Screening Parameter=',screening
'''
print "Radius =",x_bohr,"a_0"

#print x_final
#print y_final
#y_final = y_a_0*4*numpy.pi*x_a_0**2

'''
pyplot.figure()
pyplot.plot(screening,bohr_list)
pyplot.xlabel('Screening Parameter')
pyplot.ylabel('Radius of atom (a_0)')
pyplot.legend()
pyplot.show()
'''
'''
pyplot.figure()
#pyplot.plot(position,Feynmann,color='blue',label='Feynmann integral')
#pyplot.plot(position2,oscilator,color='red',linestyle='-',label='Analytical solution')
pyplot.plot(x_final,y_f)

pyplot.xlabel('Position (a_0)')
pyplot.ylabel('Probability')
pyplot.legend()
#pyplot.xlim([-7,7])
pyplot.show()
'''

pyplot.figure()
#pyplot.plot(position,Feynmann,color='blue',label='Feynmann integral')
#pyplot.plot(position2,oscilator,color='red',linestyle='-',label='Analytical solution')
pyplot.hist(G_3,bins=spatial_points,normed=True,label='Nummerical')
pyplot.plot(position2,nonradial,color='black',linestyle='--',label='Analytical')
pyplot.xlabel('Position (a_0)')
pyplot.ylabel('Probability')
pyplot.legend()
#pyplot.xlim([-7,7])
pyplot.show()

pyplot.figure()
#pyplot.plot(position,Feynmann,color='blue',label='Feynmann integral')
#pyplot.plot(position2,oscilator,color='red',linestyle='-',label='Analytical solution')
pyplot.hist(G_3,bins=spatial_points,normed=True,label='Nummerical')
pyplot.plot(position2,nonradial,color='black',linestyle='--',label='Analytical')
pyplot.xlabel('Position (a_0)')
pyplot.ylabel('Probability')
pyplot.legend()
pyplot.xlim([0,5])
pyplot.show()


'''
pyplot.figure()
pyplot.hist(accept,bins=10,normed=True)
pyplot.show()
'''

pyplot.figure()
pyplot.plot(x_final,y_final,color='red',label='moving average')
pyplot.plot(x_final,y_f,marker='x',linestyle='none',label='Nummerical')
pyplot.plot(position2,radial1,color='black',linestyle='--',label='Analytical')
pyplot.xlabel('Position (a_0)')
pyplot.ylabel('Probability')
pyplot.legend()
pyplot.show()


pyplot.figure()
pyplot.plot(x_final,y_final,color='red',label='Moving average')
pyplot.plot(x_final,y_f,marker='x',linestyle='none',label='Nummerical')
pyplot.plot(position2,radial1,color='black',linestyle='--',label='Analytical')
pyplot.xlabel('Position (a_0)')
pyplot.ylabel('Probability')
pyplot.legend()
pyplot.xlim([0,5])
pyplot.show()


'''
pyplot.figure()
pyplot.plot(position,Feynmann,color='blue',label='Feynmann integral')
pyplot.plot(position2,oscilator,color='red',linestyle='dashed',label='Analytical solution')
pyplot.xlabel('Position')
pyplot.ylabel('Probability')
pyplot.xlim([0,2])
pyplot.legend()
pyplot.show()
'''

