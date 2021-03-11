#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:39:04 2020

@author: andrewmartin
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
from datetime import datetime
import scipy.interpolate as spinterp

#import sys #import exit


#
# Use current value of Vc to update internal osmolality Ci
#
def update_Ci(Vc, C0, V0, Vb):
    Ci = C0*(V0-Vb)/(Vc-Vb)
    return Ci
    
#
# Update time-derivative of total cell volume Vc
# No solute case
#
def dVcdt_without_solute(y, t, Lp,Ci,Ce,A,R,T,C0,V0,Vb):
    #unpack y and parameters
   # [Lp, Ce, A, R, T, C0, V0, Vb] = params
    Ci = update_Ci(y[0], C0, V0, Vb)   # the catch here is that it is not a differential equaiton
    #print("Ci :", Ci, C0, V0, Vb, y[0])
    output = Lp*A*(Ci-Ce)*R*T
    #print("dVcdt :", output, Lp, A, Ci, Ce, R, T)
    return [output] 

#
# Update internal solute osmolality
#
def update_CiS(Vs, C0, V0, Vb):
    CiS = C0*(V0-Vb)/(Vs-Vb)
    return CiS

#
# update time-derivative of solute volume Vs
#
def dVsdt(y,Ps,A,CeS,C0,V0,Vb,VbarS):
    #CiS = update_CiS(y[1], C0, V0, Vb)
    CiS = y[1]/(VbarS*(y[0]-y[1]-Vb))
    #print("CiS :", CiS, y[1], C0, V0, Vb)
    dNsdt = Ps*A*(CeS-CiS)
    output = dNsdt*VbarS
    return output

#
# Solve for the change in total cell water volume
#
def dVcdt_without_solute_S(y, t, Lp,Ci,Ce,A,R,T,C0,V0,Vb,VbarS):
    Ci = update_Ci(y[0]-y[1], C0, V0, Vb)  
    CiS = y[1]/(VbarS*(y[0]-y[1]-Vb))
    output = Lp*A*(Ci+CiS-Ce)*R*T
    return [output] 

#
# solve for change in total cell volume
#
def dVcdt_solute(y, t, Lp,Ps,Ci,Ce,A,R,T,C0,V0,Vb,CeS,VbarS):
        f = dVcdt_without_solute_S(y, t, Lp,Ci,Ce,A,R,T,C0,V0,Vb,VbarS)
        g = dVsdt(y,Ps,A,CeS,C0,V0,Vb,VbarS)
        #print( "f, g", f, g)
        output = f[0] + g
        return [output,g]

class shrinkSwellToolBox():
    
    def __init__(self,ss=0):
        
        # simulation type
        self.shrinkswell = ss
        self.check_shrinkswell()
        
        #time parameters
        self.total_time = 1.0
        self.nt = 10
        self.dt = 0.0
        
        #volume parameters (what units?)
        self.Vb = 1.0
        self.V0 = 1.0
        self.VbarS = 1.0
        
        # motilities
        self.C0 = 1.0
        self.Ce = 1.0   #extracellular motility to water
        self.CeS = 1.0  #extracellular motility to CPA
                
        # cell membrane area (constant)
        self.A = 1.0
        
        #temperature
        self.T = 1.0
                
        # radius (micron)
        self.radius = 10
        
        # physical constants
        self.R = 1.0
        
        # variables that are solved dynamically
        self.Ci = 1.0   #intracellular motiliity (water)
        self.CiS = 1.0  #intracellular motility (CPA)
        self.Vc = 1.0
        self.Vs = 1.0
        
        # Fitting parameters
        self.Ps = 0.0
        self.PsStart = 0.0
        self.PsEnd = 1.0
        self.PsN = int(100)
        
        # cell membrane permiability to water
        self.Lp = 0.0
        self.LpStart = 0.0
        self.LpEnd = 1.0
        self.LpN = int(100)
        
        # Vc data
        self.VcdataFile = "none.txt"   # a text file (.txt or .csv) containging the data
        self.Vcdata = np.ones(self.LpN)
        self.errArr = np.ones(self.LpN)+0.5
        
        self.VcdataCPA = np.ones(self.PsN)
        self.errArrCPA = np.ones(self.LpN)+0.5

        
        # output parameters
        self.outpath = "None"
        self.tag = "None"
        self.logname = "log.txt"
        
        self.R = 8.207 *1e-2  # L atm K^-1 mol^-1

        self.shrinkswell = 1
        self.simORfit = 1
        self.usedata = False
        self.datapath = "None"
        self.data_filename = "None"
        self.plot_ext = "png"
        
    
    # write out all the input paramters to file
    
    # load data
    def load_data( self,fname, delimiter=",", skiprows=0):
        
        self.Vcdata = np.loadtxt( fname, delimiter=delimiter, skiprows=skiprows)

    def timepoints(self):
        self.t = np.arange(self.nt)*self.total_time/self.nt
   
    def solve_no_solute_case(self):
        
        #test
        #Ci = update_Ci(self.Vc, self.C0, self.V0, self.Vb)
        #print("Ci:", Ci)
        #Vctest = dVcdt_without_solute([self.Vc], 0.0, self.Lp,self.Ci,self.Ce,self.A,
        #                             self.R,self.T,self.C0,self.V0,self.Vb)
        #print("Vctest:", Vctest, self.Vc)
        
        y0 = self.Vc
        params = (self.Lp, self.Ci, self.Ce, self.A, self.R, self.T, self.C0, self.V0, self.Vb)
        
        self.VcArr = spi.odeint( dVcdt_without_solute, y0, self.t, args=params )
        self.VcArr = np.squeeze(self.VcArr)


    def Vc_error( self, Vcdata, VcCurrent):
        
        norm1 = np.sum(np.abs(VcCurrent)**2)
        norm2 = np.sum(np.abs(Vcdata)**2)
        error = np.sum( np.abs(VcCurrent - Vcdata)**2 )\
                    / np.sqrt( norm1*norm2)
        return error            
        

    def fitLp(self):
        
        Vcstore = self.Vc
        LpStep = (self.LpEnd - self.LpStart)/self.LpN
        errArr = np.ones(self.LpN)*10000.0
        for i in np.arange(self.LpN):
            self.Lp = self.LpStart + i*LpStep
            self.Vc = Vcstore
            self.solve_no_solute_case()
            VcArr = np.copy(self.VcArr)
            errArr[i] = self.Vc_error( self.Vcdata, VcArr )
            
#            plt.figure()
#            plt.plot( self.Vcdata)
#            plt.plot(VcArr)
#            plt.draw()
#            plt.show()
#            exit()
            
            if errArr[i] == np.min(errArr):
                VcArrStore = np.squeeze(VcArr)
                self.LpFitted = self.Lp
        self.VcArr = VcArrStore        
        self.errArr = errArr
        
        
    def solve_with_solute(self):
        
        #test
        #Ci = update_Ci(self.Vc, self.C0, self.V0, self.Vb)
        #print(Ci)
        #Vctest = dVcdt_wateronly([self.Vc], 0.0, self.Lp,self.Ci,self.Ce,self.A,
        #                             self.R,self.T,self.C0,self.V0,self.Vb)
        
        y0 = [self.Vc, self.Vs]
        params = (self.Lp, self.Ps,self.Ci, self.Ce, self.A, self.R, \
                      self.T, self.C0, self.V0, self.Vb,self.CeS,\
                          self.VbarS)
        y_out = spi.odeint( dVcdt_solute, y0, self.t, args=params )
        self.VcArrCPA = y_out[:,0]
        self.VsArr = y_out[:,1]


    def fitPs(self):
        
        PsStep = (self.PsEnd - self.PsStart)/self.PsN
        errArr = np.ones(self.PsN)*10000000.0
        for i in np.arange(self.PsN):
            self.Ps = self.PsStart + i*PsStep
            self.solve_with_solute()
            VcArrCPA = np.copy(self.VcArrCPA)
            errArr[i] = self.Vc_error( self.Vcdata, VcArrCPA )
            
            if errArr[i] == np.min(errArr):
                VcArrStore = np.squeeze(VcArrCPA)
                self.PsFitted = self.Ps
        self.VcArrCPA = VcArrStore        
        self.errArrCPA = errArr

    def check_outputpath(self):
        if self.outpath == "None":
                    print("OUTPATH has not been set.\
                          Please set a path to write the output files.")
                    raise Exception("Outpath not set")
        
        if self.tag == "None":
                    print("TAG has not been set.\
                          Please set tag to a short string to uniquely name the output files.")
                    raise Exception("Tag not set")

    
        #print("Exiting... outpath and/or tag not set.")
        
    #
    #  Writes all the input parameters to a log file
    #
    def write_all_params_to_file(self, name="None", script="fitShrinkSwell.py"):
        
        if name=="None":
            if (self.outpath == "None") or (self.tag=="None"):
                self.check_outputpath()
            f = open( self.outpath+self.tag+"_"+self.logname, 'w')
        else:
            f = open( name, 'w' )
        f.write( "# log of input parameters (shrinkSwellTools.py)\n")

        if script != "None":
            f.write( "# generated by "+script+"\n")

        a = self.__dict__
        for d, e in a.items():
            if d in ['Vcdata','VcdataCPA','errArr','errArrCPA']:
                print("found one", d)
                continue
            f.write( d+" = "+str(e)+"\n" )

        f.close()


    #
    # Plot a curve and output
    #
    def plot_and_save_fig(self, data, x=None, ylabel=None, xlabel=None,\
                          title=None,legend=None,suffix="plot", ext="png"):
        
        plt.figure()
        p = []
        for d in data:
            ptmp, = plt.plot(x,d) 
            p.append(ptmp)
        
        if ylabel!=None:
            plt.ylabel(ylabel)
            
        if xlabel!=None:
            plt.xlabel(xlabel)
            
        if title!=None:
            plt.title(title)
            
        if legend != None:
            plt.legend( p , legend)
        
        plt.savefig(self.outpath+self.tag+"_"+suffix+ext)
    
    #
    # Write data out to csv file
    #
    def write_fitted_curve_to_csv(self,data,t,fitp=None,script=None,\
                                  suffix="fit", ext="csv"):
        
        now=datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        
        self.check_outputpath()
        
       
        if len(data)!=len(t):
            print("write_fitted_curve_to_csv in shinkSwellTools.py")
            if script!=None: print("Called by :", script, "suffix:", suffix)
            print("Data points and time do not have equal number of values", len(data), len(t) )
            print("csv file of fitted data has not been written.")
        else:
            
            f = open(self.outpath+self.tag+"_"+suffix+"."+ext,'w')
        
            f.write( "# fitted data output (shrinkSwellTools.py)\n")
        
            if script!=None:
                f.write( "# generated by "+script+"\n")
                        
            f.write( "# "+dt_string+"\n")
        
            if fitp!=None:
                f.write( "# "+fitp["pname"]+" = "+str(fitp["value"])+"\n")
                f.write( "# error = "+str(fitp["error"])+"\n")
    
            for i in range(0,len(data)):
                f.write(str(t[i])+", "+str(data[i])+"\n")
    
            f.close()
            
    def check_shrinkswell(self):
        
        print(self.shrinkswell)
        if (self.shrinkswell!=1)and(self.shrinkswell!=2):
            print("Please set shrinkswell to 1 or 2")
            input("Click in the console, press enter and the code will quit")
            #sys.tracebacklimit = 0
            raise Exception('exit (ignore all the traceback text above)')
            #sys.tracebacklimit = 1000
            
    def fitp_dict( self,pname, value, units, error ):
        fitp= { "pname":pname,\
               "value":value,\
               "units":units,\
               "error":error}
        return fitp
    
    
    def run(self): #, shrinkswell,simORfit,usedata,datapath,data_filename,plot_ext):
        
        shrinkswell = self.shrinkswell
        simORfit = self.simORfit
        usedata = self.usedata
        datapath = self.datapath
        data_filename = self.data_filename
        plot_ext = self.plot_ext
        
        
        #
        # data is read from file if required and interpolated onto a regular spaced time points
        #
        if usedata == True:
            if datapath!=None:
                Vcdata = np.loadtxt(datapath+data_filename, delimiter=",") 
            else:
                Vcdata = np.loadtxt(data_filename, delimiter=",") 
            
            # subtract the minumum time
            Vcdata[:,0] += -Vcdata[0,0]
            
            if self.total_time > max(Vcdata[:,0]):
                self.total_time = max(Vcdata[:,0])    
                
        # set the time points for the simulation
        self.timepoints()
        #print(Vcdata.shape)
        #print(max(self.t),max(Vcdata[:,0]))
        
        if usedata == True:
            Vcinterp = spinterp.interp1d(Vcdata[:,0], Vcdata[:,1], kind='linear')
            self.Vcdata = Vcinterp( self.t )
            self.Vcdata *= self.V0 / self.Vcdata[0]
                
            if self.total_time > max(Vcdata[:,0]) :
                self.total_time = max(Vcdata[:,0])
            
        
        
        
        #
        # SIMULATE
        #
        if simORfit == 1:
            
            if shrinkswell==1:
                self.solve_no_solute_case()
                dataSim = self.VcArr/self.V0
        
                if not(usedata):
                    self.plot_and_save_fig([dataSim], x=self.t,\
                                       ylabel="Volume V/V_0", xlabel="Time (s)",\
                                       title="Predicted Shrinking Volume curve",\
                                       legend=None,suffix="shrink_Lp_sim_plot", ext=plot_ext)
                else:
                    d = [dataSim, self.Vcdata/self.Vcdata[0]]
                    legend = ["Simulated", "Data"]
                    self.plot_and_save_fig( d, x=self.t,\
                                       ylabel="Volume V/V_0", xlabel="Time (s)",\
                                       title="Predicted Shrinking Volume curve",\
                                       legend=legend,suffix="shrink_Lp_sim_plot", ext=plot_ext)
            
                self.write_fitted_curve_to_csv(dataSim,self.t, script="fitShrinkSwell.py",\
                                               suffix="Lp_Sim", ext="csv")
            
            elif shrinkswell==2:
                self.solve_with_solute()
                dataSim = self.VcArrCPA/self.V0
        
                if not(usedata):
                    self.plot_and_save_fig([dataSim], x=self.t,\
                                       ylabel="Volume V/V_0", xlabel="Time (s)",\
                                       title="Predicted Shrink-Swell Volume curve",\
                                       legend=None,suffix="shrinkSwell_Ps_sim_plot", ext=plot_ext)
                else:
                    d = [dataSim, self.Vcdata/self.Vcdata[0]]
                    legend = ["Simulated", "Data"]
                    self.plot_and_save_fig( d, x=self.t,\
                                       ylabel="Volume V/V_0", xlabel="Time (s)",\
                                       title="Predicted Shrink-Swell Volume curve",\
                                       legend=legend,suffix="shrinkSwell_Ps_sim_plot", ext=plot_ext)
                    
                self.write_fitted_curve_to_csv(dataSim,self.t, script="fitShrinkSwell.py",\
                                               suffix="Ps_Sim", ext="csv")
        
            else:
                raise Exception('Please set shinkswell to 1 or 2. (ignore all the traceback text above)')
        
        #
        # FIT
        #
        elif simORfit == 2:
            
            if not(usedata):
                    print("Please set usedata=True for fitting. No data has been read in.")
                    raise Exception('exit (ignore all the traceback text above)')
            
            if shrinkswell == 1:
                
                #
                # Fit the data
                #
                self.fitLp()
                fittedData = self.VcArr/self.V0
                print("Lp fitted value:", self.LpFitted*60)    
                
                #
                # plot and write the fitted Volume values to file (Lp Fit)
                #
                d = [fittedData, self.Vcdata/self.Vcdata[0]]
                legend = ["Fitted", "Data"]
                self.plot_and_save_fig( d, x=self.t,\
                                       ylabel="Volume V/V_0", xlabel="Time (s)",\
                                       title="Fitted Shrinking Volume curve",\
                                       legend=legend,suffix="shrink_Lp_fitted_plot", ext=plot_ext)
                
                fitp = self.fitp_dict( "Fitted Lp", self.LpFitted*60, "(micron/min/atm)", np.min(self.errArr) )
                self.write_fitted_curve_to_csv(fittedData,self.t,fitp=fitp,script="fitShrinkSwell.py",\
                                          suffix="Lp_fit", ext="csv")
                
                #
                # plot and write the error information to file (Lp fit)
                #
                lpvals = 60 * np.arange(self.LpN)*(self.LpEnd-self.LpStart)/self.LpN
                self.plot_and_save_fig( [self.errArr], x=lpvals,\
                                       ylabel="Lp fitting error", xlabel="Lp value (micron/min/atm)",\
                                       title="Errors in fitted Lp curve",\
                                       suffix="shrink_Lp_fitted_errors", ext=plot_ext)
                
                self.write_fitted_curve_to_csv(self.errArr,lpvals,fitp=fitp,script="fitShrinkSwell.py",\
                                          suffix="Lp_errors", ext="csv")
         
            elif shrinkswell == 2:
                   
                self.fitPs()
                fittedData = self.VcArrCPA/self.V0
                print("Fitted Ps value:", self.PsFitted*60.)
                print("Minimum error value:", np.min(self.errArrCPA))
        
                #
                # plot and write the fitted Volume values to file (Lp Fit)
                #
                d = [fittedData, self.Vcdata/self.Vcdata[0]]
                legend = ["Simulated", "Data"]
                self.plot_and_save_fig( d, x=self.t,\
                                       ylabel="Volume V/V_0", xlabel="Time (s)",\
                                       title="Predicted Shrink-Swell Volume curve",\
                                       legend=legend,suffix="shrinkSwell_Ps_fitted_plot", ext=plot_ext)
                
                fitp = self.fitp_dict( "Fitted Ps", self.PsFitted*60, "(micron/min)", np.min(self.errArrCPA) )
                self.write_fitted_curve_to_csv(fittedData,self.t,fitp=fitp,script="fitShrinkSwell.py",\
                                          suffix="Ps_fit", ext="csv")
                
                #
                # plot and write the error information to file (Ps fit)
                #
                psvals = 60 * np.arange(self.PsN)*(self.PsEnd-self.PsStart)/self.PsN
                self.plot_and_save_fig( [self.errArrCPA], x=psvals,\
                                       ylabel="Lp fitting error", xlabel="Ps value (micron/min)",\
                                       title="Errors in fitted Ps curve",\
                                       legend=legend,suffix="shrink_Ps_fitted_errors", ext=plot_ext)
                
                self.write_fitted_curve_to_csv(self.errArrCPA,psvals,fitp=fitp,script="fitShrinkSwell.py",\
                                          suffix="Ps_errors", ext="csv")
                
            else:
                raise Exception('Please set shinkswell to 1 or 2. (ignore all the traceback text above)')
        
