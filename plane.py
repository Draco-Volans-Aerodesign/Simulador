from numpy import cos, sin, tan, pi
from scipy.integrate import odeint
import numpy as np 

#velocidade do som (m/s)
vsom = 343

#conversão angular
ca = 180/pi

def gravidade(H):
        mi = 3.986e14
        RT = 6378000
        return mi/((RT + H)**2)

def densidadeAr(H):
        g0 = 9.80665
        A0 = -6.5e-3
        T0 = 288.15
        p0 = 1.225
        R  = 287.043
        return p0*((1 + A0*H/T0)**(-(1 + g0/(A0*R))))

def equations(z, t, g, Fx, Fy, Fz, L, M, N, c):
        u  = z[0]
        v  = z[1]
        w  = z[2]
        fi = z[3]
        th = z[4]
        ps = z[5]
        p  = z[6]
        q  = z[7]
        r  = z[8]
        x  = z[9]
        y  = z[10]
        z  = z[11]

        #VELOCIDADES
        du  =  r*v - q*w - g*sin(th)         + Fx
        dv  = -r*u + p*w + g*sin(fi)*cos(th) + Fy
        dw  =  q*u - p*v + g*cos(fi)*cos(th) + Fz



        #ANGULOS DE EULER
        dfi =  p + tan(th)*(q*sin(fi) + r*cos(fi))
        dth =  q*cos(fi) - r*sin(fi)
        dps = (q*sin(fi) + r*cos(fi))/cos(th)

        #VELOCIDADES ANGULARES
        dp  = (c[1]*r + c[2]*p)*q + c[3]*L + c[4]*N
        dq  =  c[5]*p*r - c[6]*(p**2 - r**2) + c[7]*M
        dr  = (c[8]*p - c[2]*r)*q + c[4]*L + c[9]*N

        #POSICÃO
        dx  =  u*cos(th)*cos(ps) + v*(sin(fi)*sin(th)*cos(ps) - cos(fi)*sin(ps)) + w*( cos(fi)*sin(th)*cos(ps) + sin(fi)*sin(ps))
        dy  =  u*cos(th)*sin(ps) + v*(cos(fi)*cos(ps) - sin(fi)*sin(th)*sin(ps)) + w*(-sin(fi)*cos(ps) + cos(fi)*sin(th)*sin(ps))
        dz  = -u*sin(th)         + v*(sin(fi)*cos(th))                           + w*( cos(fi)*cos(th))
        # if z <= 0.:
        #     z = 0.
        #     if dz <=0:
        #         dz = 0.
        #     if dw <= 0.:
        #         dw = 0.
        #     if dq <= 0.:
        #         dq = 0. 
        return [du, dv, dw, dfi, dth, dps, dp, dq, dr, dx, dy, dz]
        

class Plane:
        global vsom, ca;
        def __init__(self, deltaT):
                ##### MOMENTOS DE INERCIA (Kg m2)
                self.__I = {'xx' : 0.621, 'yy' : 0.114, 'zz' : 0.717, 'xy': 0.0, 'xz' : -0.002, 'yz' : 0.0}
                GAMA = self.__I['xx']*self.__I['zz'] - self.__I['xz']**2
                
                self.__c1 = ((self.__I['yy'] - self.__I['zz'])*self.__I['zz'] - self.__I['xz']**2)/GAMA
                self.__c2 = ( self.__I['xx'] - self.__I['yy'] + self.__I['zz'])*self.__I['xz']/GAMA
                self.__c3 =   self.__I['zz']/GAMA
                self.__c4 =   self.__I['xz']/GAMA
                self.__c5 = ( self.__I['zz'] - self.__I['xx'])/self.__I['yy']
                self.__c6 =   self.__I['xy']/self.__I['yy']
                self.__c7 =   1/self.__I['yy']
                self.__c8 = ((self.__I['xx'] - self.__I['yy'])*self.__I['xx'] + self.__I['xz']**2)/GAMA
                self.__c9 =   self.__I['xx']/GAMA

                self.__c  = (0, self.__c1, self.__c2, self.__c3, self.__c4, self.__c5, self.__c6, self.__c7, self.__c8, self.__c9)
                
                ##### INTERVALO DE TEMPO
                self.__dT = deltaT/1000
                
                ##### ANGULOS DE EULER, RADIANOS
                self.__fi    = 0
                self.__theta =0.07288146413558501
                self.__psi   = 0

                ##### ANGULOS AERODINÂMICOS, RADIANOS (apenas definindo-os)
                self.__alpha = 0.07288146413558519
                self.__beta  = 0

                self.__Dalpha = 0 #variação de alpha
                self.__Dbeta  = 0 #variação de beta

                ##### VELOCIDADES ANGULARES, RADIANOS/SEGUNDO
                self.__p = 0 #rolamento 
                self.__q = 0 #arfagem
                self.__r = 0 #guinada

                ##### VELOCIDADES EM RELAÇÃO AO EIXO DO AVIÃO, METROS/SEGUNDO
                self.__u = 0.1
                self.__v = 0
                self.__w = -0.

                ##### ACELERAÇÕES EM RELAÇÃO AO EIXO DO AVIÃO
                self.__au = 0
                self.__av = 0
                self.__aw = 0

                ##### POSIÇÃO DO CG (cordenadas globais)
                self.__x = 0
                self.__y = 0
                self.__z = -0.

                ##### CORDA GEOMÉTRICA MÉDIA
                self.__c_ = 0.609

                ##### ENVERGADURA DA ASA
                self.__b = 2.428

                ##### ANGULO DE INCLINAÇÃO DO MOTO EM RELAÇÃO AO PLANO HORIZONTAL DO AVIÃO
                self.__alpha_t = 0

                ##### ÁREA TOTAL DA ASA (m2)
                self.__S = 1.45

                ##### ÁREA DA ASA VARIDA PELO 'WASH'
                self.__SD = 0.0

                ##### MASSA DO AVIÃO (kg)
                self.__massa = 13.6

                ##### COEFICIENTES DE FORÇA (apenas definindo-os)
                self.__CL = 0 #SUSTENTAÇÂO
                self.__CD = 0 #ARRASTO
                self.__CY = 0 #SIDEFORCE

                ##### COEFICIENTES DE MOMENTO (apenas definindo-os)
                self.__Cl = 0 #ROLAMENTO
                self.__Cm = 0 #ARFAGEM
                self.__Cn = 0 #GUINADA
                
                ##### FORÇAS (apenas definindo)
                self.__Fw = [0, 0, 0]
                self.__F  = [0, 0, 0]
                self.__M  = 0.
                
                ##### ESTADO DO AVIÃO (VOO, VOO_ESTOL, DECOLAGEM, POUSO)
                self.ESTADO = 'VOO'

                ##### SALVA DADOS
               #  self.arq = {}
               #  self.arq['vel']        = open('DADOS/vel_u_v_w.txt'   , 'w')
               #  self.arq['pos']        = open('DADOS/pos_x_y_z.txt'   , 'w')
               # #self.arq['acel']       = open('DADOS/acel_x_y_z.txt'  , 'w')
               #  self.arq['vel_ang']    = open('DADOS/velA_p_q_r.txt'  , 'w')
               #  self.arq['angs_aero']  = open('DADOS/angs_ar_a_b.txt' , 'w')
               #  self.arq['forces_aero']= open('DADOS/forces_D_L_Y.txt', 'w')
               #  self.arq['angs_euler'] = open('DADOS/angs_eu_th_fi_psi.txt', 'w')

                ##### TEMPO CORRIDO PELA SIMULAÇÃO
                self.tempo = 0

        def upDate(self, empuxo, deltaA, deltaP, deltaL, deltaR, VE = (0, 0, 0)):

                self.__deltaP = deltaP
                self.__VE = VE
                ##### GRAVIDADE EM RELAÇÃO AO EIXO DO AVIÃO
                g  = 9.81
                # Ga = g*np.array([-sin(self.__theta), sin(self.__fi)*cos(self.__theta), cos(self.__fi)*cos(self.__theta)])

                ##### VENTO EXTERNO
                uE, vE, wE = self.ventoExterno(VE)
                
                ##### VELOCIDADE EM RELAÇÃO AO EIXO ABSOLUTO
                # Va = self.get_velocidades_teste([uE, vE, wE])
                
                Va = np.array([self.__u + uE, self.__v + vE, self.__w + wE], dtype = float)
                # Vte = self.get_velocidades()
                
                # print(Va)
                ##### ANGULOS AERODINÂMICOS
                # gamma_a  =  np.arctan2(-Va[2], Va[0])
                # gamma_b  =  np.arcsin(Va[1]/np.sqrt(sum(Va*Va)))

                # self.__alpha_aux = self.__theta - gamma_a
                # self.__beta_aux = self.__psi - gamma_b
                # self.__Dalpha = (self.__alpha_aux - self.__alpha)/self.__dT
                # self.__Dbeta  = (self.__beta_aux  -  self.__beta)/self.__dT

                # self.__alpha = np.copy(self.__alpha_aux)
                # self.__beta  = self.__beta_aux
                
                # print(id(self.__alpha), id(self.__alpha_aux), id(self.__theta), id(gamma_a))
                # print(np.arctan2(Vp[2], Vp[0]))
                
                alpha_ =  np.arctan2(Va[2], Va[0])
                beta_  =  np.arcsin(Va[1]/np.sqrt(sum(Va*Va)))

                self.__Dalpha = (alpha_ - self.__alpha)/self.__dT
                self.__Dbeta  = (beta_  - self.__beta )/self.__dT

                self.__alpha = alpha_
                self.__beta  = beta_ 

                ##### VELOCIDADE TOTAL
                VT = self.getVelocidade(uE, vE, wE)
                # print(VT,'t')
                ##### PAGAR OS COEFICIENTES DE FORÇA E MOMENTO
                self.__setCL(self.__alpha , (deltaP), VT)
                self.__setCD(self.__alpha , (deltaP), VT)
                self.__setCY(self.__beta  , deltaA, deltaL, VT)
                self.__setCl(self.__beta  , deltaA, deltaL, VT)
                self.__setCm(self.__alpha , (deltaP), VT)
                self.__setCn(self.__beta  , deltaA, deltaL, VT)
                
                
                # if self.__alpha > 9.25*np.pi/180 and self.__z >= -self.__b/2:
                #     self.__stall(alpha_, deltaP, VT)
                        
                # elif self.__alpha > 19*np.pi/180 and self.__z < -self.__b/2:
                #     self.__stall(alpha_, deltaP, VT)
                        
                
                ##### EMPUXO, OU TRAÇÃO EM RELAÇÃO AO EIXO NO AVIÃO
                T = empuxo*np.array([cos(self.__alpha_t), 0,  -sin(self.__alpha_t)])

                ##### FORÇAS AERODINÂMICAS EM RELAÇÃO AO EIXO DO VENTO
                q_ = self.getPressaoDinamica(uE, vE, wE) #PRESSÃO DINÂMICA
                D  =  q_*self.__S*self.__CD*np.array([-1,  0,  0]) #DRAG
                Y  =  q_*self.__S*self.__CY*np.array([ 0,  1,  0]) #SIDEFORCE
                L  =  q_*self.__S*self.__CL*np.array([ 0,  0, -1]) #LIFT

                Fw = D+L+Y
                self.__Fw = Fw
                
                ### EU ALTEREI O SINAL DE FW[0]
                ### O SINAL DE Fw[2] está certo na primeira equação? -? porque não mais?
                Fa = np.array([Fw[0]*cos(self.__alpha)*cos(self.__beta) - Fw[1]*cos(self.__alpha)*sin(self.__beta) - Fw[2]*sin(self.__alpha),\
                               Fw[0]*sin(self.__beta)             + Fw[1]*cos(self.__beta)                               ,\
                               Fw[0]*sin(self.__alpha)*cos(self.__beta) + Fw[1]*sin(self.__alpha)*sin(self.__beta) + Fw[2]*cos(self.__alpha)])
                
                     
                ##### Atrito com o solo
                if self.__z != 0:
                    REC = [0,0,0]
                else:
                    atritoDinamico = 0.044
                    REC = atritoDinamico*((self.__massa*9.81)-L[2])*np.array([cos(self.__beta),sin(self.__beta),0])
                    
                ##### MOMENTOS AERODINÂMICOS EM RE
                L_ =  q_*self.__S*self.__b *self.__Cl #ROLLING  MOMENT
                M  =  q_*self.__S*self.__c_*self.__Cm #PITCHING MOMENT
                N  =  q_*self.__S*self.__b *self.__Cn #YAWING   MOMENT
                
                #### Correcao do momento devido ao motor
                M -=  0.0072*empuxo
                
                #### Caso esteja em solo, o momento se modifica:
                if self.__z <= 0.07:
                    M -=  self.__massa*g*(0.0272) 
                

                ##### SOMA DAS FORÇAS
                F = ((Fa+T-REC)/self.__massa)
                self.__F = F*self.__massa
                
                self.__M = M
                
                #CONDICAO INICIAL
                z0 = [self.__u, self.__v, self.__w, self.__fi, self.__theta, self.__psi, self.__p, self.__q, self.__r, self.__x, self.__y, self.__z]
                tspan = [self.tempo, self.tempo + self.__dT]
                dados = (g, F[0], F[1], F[2], L_, M, N, self.__c)

                resul = odeint(equations, z0, tspan, args = dados)

                #PROXIMA ITERAÇÃO
                self.__u     = resul[1][0]
                self.__v     = resul[1][1]
                self.__w     = resul[1][2]
                self.__fi    = resul[1][3]
                self.__theta = resul[1][4]
                self.__psi   = resul[1][5]
                self.__p     = resul[1][6]
                self.__q     = resul[1][7]
                self.__r     = resul[1][8]
                self.__x     = resul[1][9]
                self.__y     = resul[1][10]
                self.__z     = resul[1][11]
                
                ## Tentando impedir o avião de entrar no chão. Se a velocidade for abaixo da decolagem ou o avião estiver no chão, ele é impedido de adquirir uma posição negativa, uma velocidade
                ## negativa vertical e de rotacionar para dentro do solo. Lembrando que o sentido do solo é +z.

                if self.__z >= -0.001 or (self.__z > -0.1 and self.__x < 47):
                    if self.__z > -0.001 or (self.__z > -0.1 and self.__x < 47):
                          self.__z = 0.
                    if self.__w > 0. :
                        self.__w = 0.
                    if self.__theta <= 0.0726732186437:
                        if self.__q < 0. :
                            self.__q = 0.
                        self.__theta = 0.0726732186437

                ### calculo do stall. Considera 0.5 como a altura onde as propriedades deixam de ser as de solo.
                # if self.__theta > 0.14 and self.__z >= -self.__b/2:
                #     self.__q = self.__q*(0.14/self.__theta)
                        
                # elif self.__theta > 0.323 and self.__z < -self.__b/2:
                #     self.__q = self.__q*(0.323/(self.__theta*1.1))
                
                
                
                
                
                ##### PASSO DE TEMPO
                self.tempo += self.__dT
                
        
        ##### Coeficientes, FORÇAS  ###################################################################################################
        ##### SUSTENTAÇÂO                                                                                                             #
        def __setCL(self, alpha, deltaP, VT):
                ## troca entre as propriedades de solo e as aereas. Considera 0.5
            
                if self.__z < -self.__b/2:
                    CLalpha   =  3.6927 #rad-1
                    CLzero    =  0.06
                else:
                    CLalpha   =  6.72*(self.__z + self.__b/2 )/(self.__b/2) + 3.6927*(self.__z)/(-self.__b/2) #rad-1
                    CLzero    =  0.124*(self.__z + self.__b/2 )/(self.__b/2) + 0.06*(self.__z)/(-self.__b/2)
                CLdelta   =  1.428 #rad-1
                CLq       =  4.58 
                CLalphaP  =  0
                self.__CL = CLzero + CLalpha*alpha + CLdelta*deltaP \
                                   + (self.__c_/VT)*(CLq*self.__q + CLalphaP*self.__Dalpha)
                                                                                                                                      #
        ##### ARRASTO                                                                                                                 #
        def __setCD(self, alpha, deltaP, VT):  
                if self.__z < -self.__b/2:
                    self.__CD =  0.0122 - 0.00279*self.__CL + 0.0703*self.__CL*self.__CL
                else:              
                    ## CD baseado na analise de elementos finitos                                                                          
                    self.__CD =  ((0.0167 - 0.00122*self.__CL + 0.0267*self.__CL*self.__CL)*(self.__z + self.__b/2 )/(self.__b/2) 
                    + (0.0122 - 0.00279*self.__CL + 0.0703*self.__CL*self.__CL)*(self.__z)/(-self.__b/2))
                                                                                                                                                      #
        ##### SIDEFORCE                                                                                                               #
        def __setCY(self, beta, deltaA, deltaL, VT):                                                                                    #
                CYbeta    = -0.191295                                                                                                       #
                CYdeltaA  = -0.046238
                CYdeltaR  =  0.000000
                CYp       = -0.042224
                CYr       =  0.084857
                self.__CY = CYbeta*beta + CYdeltaA*deltaA + CYdeltaR*deltaL\
                                        + (self.__b/VT)*(CYp*self.__p + CYr*self.__r)
        ###############################################################################################################################
        
        ##### Coeficientes, MOMENTOS ##################################################################################################
        ##### ROLAMENTO                                                                                                               #
        def __setCl(self, beta, deltaA, deltaL, VT):                                                                                    #
                Clbeta    = -0.127056                                                                                                      #
                Clp       = -0.404480
                Clr       =  0.103636
                CldeltaA  = -0.273702
                CldeltaR  =  0.012953
                self.__Cl = Clbeta*beta + (self.__b/VT)*(Clp*self.__p + Clr*self.__r) \
                                        + CldeltaA*deltaA + CldeltaR*deltaL                    #
                                                                                                                                      #
        ##### ARFAGEM                                                                                                                 #
        def __setCm(self, alpha, deltaP, VT):
                if self.__z < -self.__b/2:
                    Cmalpha   = -0.18322
                    Cmzero    =  0.002776
                else:
                    Cmalpha   = -0.5449*(self.__z + self.__b/2 )/(self.__b/2) + -0.18322*(self.__z)/(-self.__b/2)
                    Cmzero    = -0.003964*(self.__z + self.__b/2 )/(self.__b/2) + 0.002776*(self.__z)/(-self.__b/2)
                Cmq       = -1.020635 #-0.015401 
                CmalphaP  =  0.0
                Cmdelta   = -0.477217
                self.__Cm = Cmalpha*alpha + Cmzero + (self.__c_/VT)*(Cmq*self.__q + CmalphaP*self.__Dalpha)\
                                          + Cmdelta*deltaP
                
        ##### GUINADA                                                                                                                 #
        def __setCn(self, beta, deltaA, deltaL, VT):                                                                                    #
                Cnbeta   =  0.022216                                                                                                        #
                Cnp      = -0.013246
                Cnr      = -0.012476
                CndeltaA = -0.001604
                CndeltaR = -0.006818
                self.__Cn = Cnbeta*beta + (self.__b/VT)*(Cnp*self.__p + Cnr*self.__r) \
                                        + CndeltaA*deltaA + CndeltaR*deltaL                     #
        ###############################################################################################################################
        
        def __stall(self, alpha, deltaP, VT):
            ### Prototipo. A ideia é espelhar o comportamento aerodinamico em relação ao alpha máximo.
            
            if self.__z < -self.__b/2:
                CLalpha   =  3.6927 #rad-1
                CLzero    =  0.06
                alpha =  2*9.25*np.pi/180 - alpha
            else:
                CLalpha   =  6.72*(self.__z + self.__b/2)/(self.__b/2) + 3.6927*(self.__z)/(-self.__b/2) #rad-1
                CLzero    =  0.124*(self.__z + self.__b/2 )/(self.__b/2) + 0.06*(self.__z)/(-self.__b/2)
                alpha =  2*19*np.pi/180 - alpha
                if alpha > np.pi/2:
                    alpha = np.pi/2
            CLdelta   =  1.42655 #rad-1
            CLq       =  4.56392 
            CLalphaP  =  0
            
            
            self.__CL = CLzero + CLalpha*alpha + CLdelta*deltaP \
                               + (0.5*self.__c_/VT)*(CLq*self.__q + CLalphaP*self.__Dalpha)          

            if self.__z < -self.__b/2:
                Cmalpha   = -0.512797
                Cmzero    = -0.003704
            else:
                Cmalpha   = -0.5449*(self.__z + self.__b/2 )/(self.__b/2) + -0.18322*(self.__z)/(-self.__b/2)
                Cmzero    = -0.003964*(self.__z + self.__b/2 )/(self.__b/2) + 0.002776*(self.__z)/(-self.__b/2)
            Cmq       = -1.017122 #-0.015401 
            CmalphaP  =  0.0
            Cmdelta   = -0.477217
            self.__Cm = Cmalpha*alpha + Cmzero + (self.__c_/VT)*(Cmq*self.__q + CmalphaP*self.__Dalpha)\
                                      + Cmdelta*deltaP
                                      
        def ventoExterno(self, VE):
                th = self.__theta
                fi = self.__fi
                ps = self.__psi
                uE  =  VE[0]*cos(th)*cos(ps) + VE[1]*(sin(fi)*sin(th)*cos(ps) - cos(fi)*sin(ps)) + VE[2]*( cos(fi)*sin(th)*cos(ps) + sin(fi)*sin(ps))
                vE  =  VE[0]*cos(th)*sin(ps) + VE[1]*(cos(fi)*cos(ps) - sin(fi)*sin(th)*sin(ps)) + VE[2]*(-sin(fi)*cos(ps) + cos(fi)*sin(th)*sin(ps))
                wE  = -VE[0]*sin(th)         + VE[1]*(sin(fi)*cos(th))                           + VE[2]*( cos(fi)*cos(th))
                return uE, vE, wE
        
        def getVelocidade(self, uE, vE, wE):
                return np.sqrt((self.__u + uE)**2 + (self.__v + vE)**2 + (self.__w + wE)**2)

        def getPressaoDinamica(self, uE, vE, wE):
                pp = 1.113
                return 0.5*pp*(self.getVelocidade(uE, vE, wE)**2)

        def getCoeficienteEmpuxo(self, empuxo, uE, vE, wE):
                q_ = self.getPressaoDinamica(uE, vE, wE)
                return empuxo/(q_*self.__SD)

        def getNumeroMach(self, uE, vE, wE):
                return self.getVelocidade(uE, vE, wE)/vsom

        def salva_dados(self):
                self.arq['vel'].write('%f  %f  %f\n'%(self.__u, self.__v, self.__w))
                self.arq['pos'].write('%f  %f  %f\n'%(self.__x, self.__y, self.__z))
                #self.arq['acel'].write('%f  %f  %f\n'%(self.__au, self.__av, self.__aw))
                self.arq['vel_ang'].write('%f  %f  %f\n'%(self.__p, self.__q, self.__r))
                self.arq['angs_aero'].write('%f  %f\n'%(self.__alpha, self.__beta))
                self.arq['forces_aero'].write('%f %f %f\n'%(self.__Fw[0], self.__Fw[2], self.__Fw[1]))
                self.arq['angs_euler'].write('%f %f %f\n'%(self.__theta, self.__psi, self.__fi))
                self.arq['corrente'].write('%f \n'%(self.get_servo()))

        def fecha_arquivos(self):
                for nome in self.arq:
                        self.arq[nome].close()

        def get_velocidades(self):
                return self.__u, self.__v, self.__w

        def get_velocidades_teste(self, VE):
                th = self.__theta
                fi = self.__fi
                ps = self.__psi
                
                u  =  (self.__u  + VE[0])*cos(th)*cos(ps) + (self.__v+ VE[1])*(sin(fi)*sin(th)*cos(ps) - cos(fi)*sin(ps)) + (self.__w + VE[2])*( cos(fi)*sin(th)*cos(ps) + sin(fi)*sin(ps))
                v  =  (self.__u  + VE[0])*cos(th)*sin(ps) + (self.__v+ VE[1])*(cos(fi)*cos(ps) - sin(fi)*sin(th)*sin(ps)) + (self.__w + VE[2])*(-sin(fi)*cos(ps) + cos(fi)*sin(th)*sin(ps))
                w  =  -(self.__u + VE[0])*sin(th)         + (self.__v+ VE[1])*(sin(fi)*cos(th))                           + (self.__w + VE[2])*( cos(fi)*cos(th))
                return np.array([u,v,w], dtype = float)

        def get_pos(self):
                return self.__x, self.__y, self.__z

        def get_acel(self):
                return self.__F
            
        def get_momento(self):
            return self.__M

        def get_cm(self):
            return self.__Cm
        
        def get_vel_ang(self):
                return self.__p, self.__q, self.__r

        def get_angs(self):
                return self.__alpha, self.__beta

        def get_forces_aero(self):
                return self.__Fw

        def get_euler(self):
                return self.__theta, self.__psi, self.__fi

        def fator_carga(self):
                g  = gravidade(0)
                return -self.__Fw[2]/(self.__massa*g)

        def set_interval_integrate(self, dt):
                self.__dT = dt/10
                
        def get_servo(self):
            VE = self.__VE
            forca = 1.428*self.__deltaP*self.__S*self.__CL*self.getPressaoDinamica(VE[0],VE[1],VE[2])
            corrente = 0.4 + (0.4286*9.81)*forca*0.005
            
            return corrente

                
