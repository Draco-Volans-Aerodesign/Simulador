import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
from plane import *
from numpy import (roots, polyfit, linspace, sqrt, pi)
import math
from scipy.spatial import transform as spt
import pygame
from pygame.locals import *
import matplotlib.style as mplstyle
mplstyle.use('fast')


import OpenGL
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *


P = Plane(0.02)
T = []
a = []
b = []
c = []
a0 = []
x = []
y = []
h = []
ft_carga = []
u = []
v = []
w = []
p = []
q = []
r = []

fig, axs = plt.subplots(7, 1, sharex = True)
# fig = plt.figure()
# axs = p3.Axes3D(fig)
# axs.set_xlim(-2000,2000)
# axs.set_ylim(-2000,2000)
# axs.set_zlim(0,100)
# axs.autoscale(False)

empuxo = 31
Profundor = -0.02217
Aileron = 0
Leme = 0
Pezao = 0
j = 0
comando1 = 0
comando2 = 0
comando3 = 0

Pf = []
Ar = []
Le = []
Dr = []

Vento_Externo = []
# VEE = np.random.rand(1)[0]
# if 0 <= VEE < 0.27:
#     Vento_Externo.append(np.random.rand(1)[0]*1)
# elif 0.27 <= VEE < 0.89:
#     Vento_Externo.append(np.random.rand(1)[0]*6 + 1)
# elif 0.89 <= VEE <= 1:
#     Vento_Externo.append(np.random.rand(1)[0]*15 + 7  ) 


P.upDate(empuxo, 0, Profundor, 0, 0)
k = 0


def animate(i):
    global Vento_Externo
    # if P.tempo % 5 == 0 or i == 0:
    #     VEE = np.random.rand(1)[0]
    #     if 0 <= VEE < 0.27:
    #         Vento_Externo.append((np.random.rand(1)[0]*1)/1.944* np.cos(np.random.rand(1)[0]*pi*2))
    #     elif 0.27 <= VEE < 0.89:
    #         Vento_Externo.append((np.random.rand(1)[0]*6 + 1)/1.944* np.cos(np.random.rand(1)[0]*pi*2))
    #     elif 0.89 <= VEE <= 1:
    #         Vento_Externo.append((np.random.rand(1)[0]*15 + 7)/1.944 * np.cos(np.random.rand(1)[0]*pi*2))     
    # else:
    #     Vento_Externo.append(Vento_Externo[-1])
    Vento_Externo.append(0)
    # Vestol = 11.62135069247685*1.05
    Throttle = 1
    # Throttle = 7.744609370292276/trac_disp(P.getVelocidade(0, 0, 0))
    global Profundor
    global Aileron
    global Pezao
    global comando1
    global comando2
    global comando3
    # P.salva_dados()
    T.append(P.tempo)
    a.append(P.get_euler()[0])
    b.append(P.get_euler()[1])
    c.append(P.get_euler()[2])
    a0.append(P.get_angs()[0])
    x.append(P.get_pos()[0])
    y.append(P.get_pos()[1])
    h.append(-P.get_pos()[2])
    # ft_carga.append(P.fator_carga())
    u.append(P.get_velocidades_teste([3.7*np.sin(b[-1]),3.7*np.cos(b[-1]),0])[0])
    v.append(P.get_velocidades()[1])
    w.append(-P.get_velocidades_teste([3.7*np.sin(b[-1]),3.7*np.cos(b[-1]),0])[2])
    p.append(P.get_vel_ang()[0])
    q.append(P.get_vel_ang()[1])
    r.append(P.get_vel_ang()[2])
    P.set_interval_integrate(0.3)
    
    
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            plt.close()
            
            
        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_q:
                comando1 -= 2*np.pi/180
            if event.key == pygame.K_w:
                comando1 -= 0.5*np.pi/180
            if event.key == pygame.K_e:
                comando1 += 0.5*np.pi/180
            if event.key == pygame.K_r:
                comando1 += 2*np.pi/180
                
            if event.key == pygame.K_a:
                comando2 -= 2*np.pi/180
            if event.key == pygame.K_s:
                comando2 -= 0.5*np.pi/180
            if event.key == pygame.K_d:
                comando2 += 0.5*np.pi/180
            if event.key == pygame.K_f:
                comando2 += 2*np.pi/180  
                
            if event.key == pygame.K_z:
                comando3 -= 2*np.pi/180
            if event.key == pygame.K_x:
                comando3 -= 0.5*np.pi/180
            if event.key == pygame.K_x:
                comando3 -= 0.05*np.pi/180
            if event.key == pygame.K_c:
                comando3 += 0.05*np.pi/180
            if event.key == pygame.K_c:
                comando3 += 0.5*np.pi/180
            if event.key == pygame.K_v:
                comando3 += 2*np.pi/180  
                

             
    Pf.append(Profundor - comando1)   
    Ar.append(Aileron - comando2)
    Le.append(Leme - comando3)     
        
    
    ### os eixos do opengl estão trocados
    
    glRotatef(p[-1]*0.2, 1, 0, 0)
    glRotatef(q[-1]*0.2, 0, 0, 1)
    glRotatef(r[-1]*0.2, 0, 1, 0)
    glTranslatef(u[-1]*0.2, w[-1]*0.2,v[-1]*0.2)

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    Cube()
    pygame.display.flip() 
    
    if i%1 == 0:  
        # axs.clear()
        # axs.plot(x,y,h)
        
        axs[0].clear()
        axs[1].clear()
        axs[2].clear()
        axs[3].clear()
        axs[4].clear()
        axs[5].clear()
        axs[6].clear()
        axs[0].set_title('Profundor')
        axs[1].set_title('Alpha')
        axs[2].set_title('X')
        axs[3].set_title('Altura')
        axs[4].set_title('Y')
        axs[5].set_title('Vh')
        axs[6].set_title('Velocidade angular')
        axs[0].plot(T, np.array(Le)*180/pi, label = 'Profundor')
        axs[1].plot(T, np.array(a)*180/pi)
        axs[1].plot(T, np.array(a0)*180/pi)
        axs[2].plot(T, x)
        axs[3].plot(T, h)
        axs[4].plot(T, y)  
        axs[5].plot(T, w)
        axs[6].plot(T, np.array(r)*180/pi)
    

    empuxo = trac_disp(P.getVelocidade(3.7*np.sin(b[-1]),3.7*np.cos(b[-1]),0))*Throttle
    if empuxo <= 0:
        empuxo = 0.


        
    P.upDate(empuxo, Ar[-1], Pf[-1], Le[-1], 0, VE = (3.7*np.sin(b[-1]),3.7*np.cos(b[-1]),0))

    
def trac_disp(x):   
    trac = np.polyval([-0.000159, 0.009837, -0.230541, 1.191796, 30.026713], x)
    return trac


vertices = (
    (10, -10, -10),
    (10, 10, -10),
    (-30, 10, -10),
    (-30, -10, -10),
    (10, -10, 10),
    (10, 10, 10),
    (-30, -10, 10),
    (-30, 10, 10)
    )

edges = (
    (0,1),
    (0,3),
    (0,4),
    (2,1),
    (2,3),
    (2,7),
    (6,3),
    (6,4),
    (6,7),
    (5,1),
    (5,4),
    (5,7)
    )


def Cube():   
    glBegin(GL_LINES)
    for edge in edges:
        for vertex in edge:
            glVertex3fv(vertices[vertex])
    glEnd()
    

def inicializar():  
    pygame.init()
    display = (800,600)
    screen = pygame.display.set_mode(display, pygame.DOUBLEBUF|pygame.OPENGL)
    gluPerspective(60, (display[0]/display[1]), 0.1, 5000.0)
    glTranslatef(-1500,0.0, -1500)
    
inicializar()



ani = animation.FuncAnimation(fig, animate, interval = 1)

