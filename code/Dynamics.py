import numpy as np
from functions import f1,f2,fz,History
import zipfile
import os




class Dynamics:
    

    def __init__(self,model,Numeric_method,Network):
        
        if  model["model_name"]=="epileptor" and model["model_dimension"]==6: 
            self.I1=model["I1"]
            self.I2=model["I2"]
            self.r=model["r"]
            self.gamma=model["gamma"]
            self.tau2=model["tau2"]
            self.m=model["m"]
            self.x00=model["x0"] 
        elif model["model_name"]=="epileptor" and model["model_dimension"]==2:
            self.I1=model["I1"]
            self.r=model["r"]

        self.model_name=model["model_name"]
        self.model_dim=model["model_dimension"]
        self.Numeric_Method=Numeric_method

        self.h=Network["h"]
        self.t0=Network["t0"]
        self.W=Network["W"]
        self.Tau=Network["Tau"]
        self.N=Network["N"]
        self.M=Network["M"]
        self.EZs=Network["EZs"]
        self.d_down=20
        self.save_time_series=False
        self.path_to_time_series="./"
        self.file_name="temp.txt"
        self.M_down=int(self.M/self.d_down)+1
        self.mHist=0


        self.x0=np.ones(self.N)*(-2.2)
        self.t=self.t0        
        self.t1=self.t0
        self.t2=0
        self.k0=np.zeros([self.t0+self.mHist+1,self.N])   
        self.X=np.zeros([self.model_dim,self.t0+self.mHist+1,self.N])
        self.Y=np.zeros([self.M_down,self.N])
        self.T_axis=np.zeros([self.M_down])
        self.hist=np.zeros(self.N)
        self.stimulus=np.zeros([self.model_dim,self.t0+self.mHist+1,self.N])
        self.XX=np.zeros([self.model_dim,self.N])  
        self.dX=np.zeros([self.model_dim,self.N])
        self.noise=np.zeros([self.model_dim,self.t0+self.mHist+1,self.N])
        self.stds=np.zeros([self.model_dim])
        self.seed=0
        self.refractory=False
        self.seizure_status=np.zeros(self.N)
        self.std_x1=np.zeros(self.N)
        self.refractory_time=np.zeros(self.N)
        self.termination_time=np.ones(self.N)*self.M*self.h
        self.seizure_01=np.zeros(self.N).astype(int)
        self.seizure_onsets=np.zeros(self.N)
        self.q=0
        

        if self.Numeric_Method=="Runge_Kutta_2":  

            self.k1=np.zeros([self.model_dim,self.t0+self.mHist+1,self.N]) 
            self.k2=np.zeros([self.model_dim,self.t0+self.mHist+1,self.N])

        elif self.Numeric_Method=="Euler":  

            self.k1=np.zeros([self.model_dim,self.t0+self.mHist+1,self.N])
        

    def One_step(self):
        
        if self.Numeric_Method=="Runge_Kutta_2" and self.model_dim==6:
            self.RK2()
            self.X[:,self.t2,:]=self.X[:,self.t1,:]+(self.h/4.0)*\
                ( self.k1[:,self.t1,:] + 3*self.k2[:,self.t1,:] ) + np.sqrt(self.h)*self.noise[:,self.t1,:]
            if(np.mod(self.t,self.d_down)==0):
                t_down=int(self.t/self.d_down)
                self.Y[t_down,:]=self.X[3,self.t2,:]-self.X[0,self.t2,:]
                self.T_axis[t_down]=self.t*self.h

        elif self.Numeric_Method=="Euler" and self.model_dim==6: 
            self.Euler()   
            self.X[:,self.t2,:]=self.X[:,self.t1,:]+(self.h)*\
                ( self.k1[:,self.t1,:]  ) + np.sqrt(self.h)*self.noise[:,self.t1,:]
            if(np.mod(self.t,self.d_down)==0):
                t_down=int(self.t/self.d_down)
                self.Y[t_down,:]=self.X[3,self.t2,:]-self.X[0,self.t2,:]
                self.T_axis[t_down]=self.t*self.h



     

    def run(self,t_EZ_onset,EZ_x0):
        if self.save_time_series==True:
            f=open(self.path_to_time_series+self.file_name+'.txt','w')
            g=open(self.path_to_time_series+"T_"+self.file_name+".txt",'w')

        np.random.seed(self.seed)
        for i in range(self.model_dim):
            self.noise[i]=np.random.normal(0.0, self.stds[i],  [self.t0+self.mHist+1,self.N])
        while (self.t<self.M):
            
            if self.t == int(0.1*t_EZ_onset/self.h) + self.t0 :
                self.x0=self.x00
            if self.t == int(t_EZ_onset/self.h) + self.t0 :    
                self.x0[self.EZs]=EZ_x0
                

            if self.refractory==True and self.t>int(t_EZ_onset/self.h) + self.t0:
                self.x0_to_refractory(t_EZ_onset)
            
            self.One_step()
            for i in range(self.model_dim):
                self.noise[i,self.t2,:]=np.random.normal(0.0, self.stds[i],  [self.N])
            
            if self.save_time_series==True and self.t>int(t_EZ_onset/self.h)+self.t0:
                np.savetxt(f, self.X[3,self.t2,:]-self.X[0,self.t2,:], fmt='%1.6f' ,newline=" ")
                f.write("\n")
                np.savetxt(g, np.array([self.t*self.h]), fmt='%1.6f' , newline=" ")
                g.write("\n")
                
            self.t1=self.t1+1
            self.t2=self.t2+1
            if(self.t1>self.t0+self.mHist):
                self.t1=0
            if(self.t2>self.t0+self.mHist):
                self.t2=0
            if(self.t==self.M-1 and self.q==0):
                a=np.where(self.seizure_01==1)
                self.q=1
                if(np.size(a)>0):
                    self.Y=np.append( self.Y , np.zeros([int(0.5*self.M_down)+1,self.N]) ,axis=0)
                    self.T_axis=np.append( self.T_axis , np.zeros([int(0.5*self.M_down)+1]) )
                    self.M=int(1.5*self.M)
            self.t=self.t+1
        if self.save_time_series==True:
            f.close()
            g.close()
           

              

        
    def Euler(self):
    
        self.XX=self.X[:,self.t1,:]
        self.hist=History(self.X[0,:,:],self.k0,\
            self.noise[0,:,:],0,self.Tau,self.W,self.t0+self.mHist,self.h,self.t1)
        self.F()
        self.k1[:,self.t1,:]=self.dX


    def RK2(self):
        
        
        self.XX=self.X[:,self.t1,:]
        self.hist=History(self.X[0,:,:] , self.k0[:,:] , self.noise[0,:,:] , 0.0 , self.Tau , self.W , self.t0+self.mHist , self.h , self.t1 )
        
        self.F()
        self.k1[:,self.t1,:]=self.dX

        self.XX=self.X[:,self.t1,:] + self.h * self.k1[:,self.t1,:] + np.sqrt(self.h) * self.noise[:,self.t1,:]
        self.hist=History(self.X[0,:,:] , self.k1[0,:,:] , self.noise[0,:,:] , 2.0/3.0 , self.Tau , self.W , self.t0+self.mHist , self.h , self.t1)
        
        self.F()
        self.k2[:,self.t1,:]=self.dX
        

       
        
    def F(self):
        if self.model_name=="epileptor" and self.model_dim==6:
            
            # X[0]=x1 , X[1]=y1 , X[2]=z , X[3]=x2 , X[4]=y2, X[5] 
            self.dX[0,:] = self.XX[1] - f1(self.XX[0],self.XX[3],self.XX[2]) - self.XX[2] + self.I1 #+self.stimulus[0,self.t,:]
            self.dX[1,:] = 1.0 - 5.0 * self.XX[0]**2 - self.XX[1] #+ self.stimulus[1,self.t,:]
            self.dX[2,:]= (self.r) * ( 4.0* (self.XX[0]-self.x0) - self.XX[2] - self.hist ) #+ self.stimulus[2,self.t,:] ) 
            self.dX[3,:]= self.XX[3] - self.XX[4] - self.XX[3]**3 + self.I2 + 0.002 * self.XX[5] - 0.3 * ( self.XX[2] -3.5 )# + self.stimulus[3,self.t,:]
            self.dX[4,:]= (1.0/self.tau2) * ( f2( self.XX[3] ) - self.XX[4]  ) #+ self.stimulus[4,self.t,:] )
            self.dX[5,:]= self.XX[0] -  self.gamma*self.XX[5]
           
        elif self.model_name=="epileptor" and self.model_dim==2:
            self.dX[0,:] = -self.XX[0]**3 - 2*self.XX[0]**2 + 1 - self.XX[1] + self.I1 #+self.stimulus[0,self.t,:]
            self.dX[1,:]= (self.r) * ( 4.0* (self.XX[0]-self.x0) - self.XX[1] - self.hist ) #+ self.stimulus[1,self.t,:] )

            
            
    
    def Set_Initial_Condition(self,Y0):
        for i in range(self.model_dim):
            self.X[i,0:self.t0+1,:]=Y0[i]
        self.Y[0:int(self.t0/self.d_down)+1,:]=Y0[3]-Y0[1]

      

    def Set_noise(self,noise_stds,seed):
        self.stds=noise_stds
        self.seed=seed
        
 

    def Set_stimulus(self,t_stim_on,d_stim,stimulus_strength):
        t_stim_off=t_stim_on+d_stim
        M_stim_on=int(t_stim_on/self.h)
        M_stim_off=int(t_stim_off/self.h)
        self.stimulus[0,M_stim_on:M_stim_off,self.EZs]=stimulus_strength[0]
        print(self.stimulus.shape)
        

    def x0_to_refractory(self,t_EZ_onset):
        e=1
        th_of=-80.0
        th_on=-120.0
        if(self.t>self.t0+int(t_EZ_onset/self.h)):
            a=np.where(self.seizure_01==0)
            
            for i in a[0]:
                if(self.X[5,self.t1,i]>th_on):
                    self.seizure_01[i]=1
                    self.seizure_onsets[i]=self.t*self.h

            a=np.where(self.seizure_01==1)
            for i in a[0]:
                if(self.X[5,self.t1,i]<th_of+e and self.X[5,self.t1,i]>th_of-e):
                    t3=self.t2-50
                    if(t3<0):
                        t3=self.t0+t3
                    if(self.X[5,t3,i]-self.X[5,self.t1,i]>0.01*e):
                        self.seizure_01[i]=2
                        self.termination_time[i]=self.t
                        self.x0[i]=-2.2
                        self.refractory_time[i]=(self.t)*self.h
                        self.W[i,:]=0
                        self.W[:,i]=0
   





 # def x0_to_refractory1(self,t_EZ_onset):
    #     # e=2.0
    #     # for i in range(self.N):
    #     #     if (self.refractory_time[i]<0.0001 and self.t>self.t0+int(t_EZ_onset/self.h)):
    #     #         if(self.X[5,self.t-1,i]<-e and self.X[5,self.t-1,i]>-2*e):
    #     #             if(np.mean(self.X[0,self.t-100:self.t-1,i])-self.gamma*self.X[5,self.t-1,i]<0):
    #     #                 self.x0[i]=-3.0
    #     #                 # self.X[:,self.t,:]=self.X[:,self.t0,:]
    #     #                 self.refractory_time[i]=(self.t)*self.h
    #     #                 # self.W[i,:]=0
    #     #                 # self.W[:,i]=0


    #     time_bin=300
    #     m=int(time_bin/self.h)
    #     std_th=0.08
    #     t_reset_to=int(t_EZ_onset/self.h)-100
        
    #     for i in range(self.N):
    #         if (self.refractory_time[i]<0.0001):
    #             self.std_x1[i]=np.std(self.X[0,self.t-m:self.t,i])
    #             if self.std_x1[i]>std_th and self.seizure_status[i]==0:
    #                 self.seizure_status[i]=1
    #             elif self.std_x1[i]<=std_th and self.seizure_status[i]==1:
    #                 self.seizure_status[i]=2
    #                 self.x0[i]=-3.0
    #                 self.X[:,self.t,:]=self.X[:,t_reset_to,:]
    #                 self.refractory_time[i]=(self.t)*self.h
    #     #         # print(self.t*self.h)


    # def x0_to_refractory2(self,t_EZ_onset):
    #     e=1
    #     if(self.t>self.t0+int(t_EZ_onset/self.h)):
    #         a=np.where(self.seizure_01<0.5)
            
    #         for i in a[0]:
                
    #             if(self.X[5,self.t-1,i]<e and self.X[5,self.t-1,i]>-e):
    #                 if(np.mean(self.X[0,self.t-100:self.t-1,i])-self.gamma*self.X[5,self.t-1,i]<0):
    #                     self.seizure_01[i]=1
    #                     self.termination_time[i]=self.t

    #         a=np.where(self.seizure_01>0.5)
    #         for i in a[0]:
    #             if(self.refractory_time[i]>0.001 and self.h*(self.t-self.termination_time[i])>300):
    #                 self.x0[i]=-3.5
    #                 self.refractory_time[i]=(self.t)*self.h
    






    # def RK4(self): 
    

    #     self.XX=self.X[:,self.t-1,:]
    #     self.hist=History(self.X[0,self.t-self.t0:self.t,:],self.k0,\
    #         self.noise[0,self.t-self.t0:self.t,:],0,self.Tau,self.W,self.t0,self.h)
    #     self.F()
    #     self.k1[:,self.t-1,:]=self.dX

    #     self.XX=self.X[:,self.t-1,:]+0.5*self.h*self.k1[:,self.t-1,:] +0.5* np.sqrt(self.h)*self.noise[:,self.t-1,:]
    #     self.hist=History(self.X[0,self.t-self.t0:self.t,:],self.k1[0,self.t-self.t0:self.t,:],\
    #         self.noise[0,self.t-self.t0:self.t,:],0.5,self.Tau,self.W,self.t0,self.h)
    #     self.F()
    #     self.k2[:,self.t-1,:]=self.dX

    #     self.XX=self.X[:,self.t-1,:]+0.5*self.h*self.k2[:,self.t-1,:] +0.5* np.sqrt(self.h)*self.noise[:,self.t-1,:]
    #     self.hist=History(self.X[0,self.t-self.t0:self.t,:],self.k2[0,self.t-self.t0:self.t,:],\
    #         self.noise[0,self.t-self.t0:self.t,:],0.5,self.Tau,self.W,self.t0,self.h)
    #     self.F()
    #     self.k3[:,self.t-1,:]=self.dX

    #     self.XX=self.X[:,self.t-1,:]+self.h*self.k3[:,self.t-1,:] + np.sqrt(self.h)*self.noise[:,self.t-1,:]
    #     self.hist=History(self.X[0,self.t-self.t0:self.t,:],self.k3[0,self.t-self.t0:self.t,:],\
    #         self.noise[0,self.t-self.t0:self.t,:],1.0,self.Tau,self.W,self.t0,self.h)
    #     self.F()
    #     self.k4[:,self.t-1,:]=self.dX
