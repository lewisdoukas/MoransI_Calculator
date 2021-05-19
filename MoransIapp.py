#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#+----------+-----------------------------------------------------------------+
#|   TITLE  | Moran's Index Calculator                                        |
#+----------+-----------------------------------------------------------------+
#|  AUTHOR  | Ilias Doukas                                                    |
#+----------+-----------------------------------------------------------------+
#|  CONTACT | hliasduke@gmail.com                                             |
#+----------+-----------------------------------------------------------------+
#|  DETAILS |This GUI application calculates the spatial autocorrelation using|
#|          |Moran's I, depending on various ways of calculating spatial lags.|
#+----------+-----------------------------------------------------------------+

import geopandas as gpd
import matplotlib, tooltip
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats 
import tkinter as tk
from pandastable import Table
from tkinter import ttk
from PIL import Image, ImageTk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import (NavigationToolbar2Tk)
matplotlib.use('TkAgg')



class MainWindow:
    def __init__(self, master):
        self.master = master
        self.master.geometry("500x500+400+50")
        self.master.title(" Moran's I ")
        self.frame = tk.Frame(self.master, bg= '#565051')
        self.frame.place(relheight=1, relwidth=1)
        
        self.lfeat = tk.StringVar(value='')
        self.labelfeat = tk.Label(self.frame, textvariable= self.lfeat, bg= '#565051')
        self.labelfeat.place(relx=0, rely=0.12, relwidth= 1, relheight= 0.05)
        self.feat = ''
        self.method = ''
        self.neighbormethod = ''
        self.neighN = 0
        self.buffer = tk.DoubleVar(value= 0)
        self.Dist = 0
        self.power = 0
        self.MoransI = 0
        self.perm = 0
               
        self.imageFolder = ImageTk.PhotoImage(Image.open('./imgs/folder.png'))
        self.imageAttrTable = ImageTk.PhotoImage(Image.open('./imgs/table.png'))
        self.imageOutliers = ImageTk.PhotoImage(Image.open('./imgs/outliers.png'))
        self.imageMethod = ImageTk.PhotoImage(Image.open('./imgs/method.png'))
        self.imageDistance = ImageTk.PhotoImage(Image.open('./imgs/distance.png'))
        self.imageCalc = ImageTk.PhotoImage(Image.open('./imgs/calc.png'))
        self.imageResults = ImageTk.PhotoImage(Image.open('./imgs/atr.png'))
        self.imageSave = ImageTk.PhotoImage(Image.open('./imgs/save.png'))
        self.imageRefresh = ImageTk.PhotoImage(Image.open('./imgs/refresh.png'))
        self.imageHelp = ImageTk.PhotoImage(Image.open('./imgs/help.png'))    
        self.imageOk = ImageTk.PhotoImage(Image.open('./imgs/ok.png'))
        self.imageBack = ImageTk.PhotoImage(Image.open('./imgs/back.png'))
        self.imageOutliers2 = ImageTk.PhotoImage(Image.open('./imgs/outliers2.png'))
      
        self.button = tk.Button(self.frame, image= self.imageFolder, command=lambda: self.readshp(), bg= '#565051', activebackground= '#565051', relief='flat', highlightthickness=0, bd=0)
        self.button.place(relx=0.015, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button, "Open file", 40, 45)
        
        self.button2 = tk.Button(self.frame, image= self.imageAttrTable, command=lambda: self.callwinFeat(),  state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button2.place(relx=0.135, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button2, "Select attribute", 40, 45)
        
        self.button3 = tk.Button(self.frame, image= self.imageOutliers, command=lambda: self.callWinOutliers(), state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button3.place(relx=0.255, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button3, "Check for outliers", 40, 45)
        
        self.button4 = tk.Button(self.frame, image= self.imageMethod,  command=lambda: self.callwinMethod(), state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button4.place(relx=0.375, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button4, "Select method", 40, 45)
        
        self.button5 = tk.Button(self.frame, image= self.imageDistance,  command=lambda: self.callwinNeighbor(), state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button5.place(relx=0.495, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button5, "Select number of neighbors - distance", 40, 45)
       
        self.button6 = tk.Button(self.frame, image= self.imageCalc, command=lambda: self.calcMoransI(), state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button6.place(relx=0.615, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button6, "Calculate Moran's I", 40, 45)
        
        self.button7 = tk.Button(self.frame, image= self.imageResults,  command=lambda: self.callwinResults(), state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button7.place(relx=0.735, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button7, "Show results", 40, 45)
        
        self.button8 = tk.Button(self.frame, image= self.imageRefresh, bg= '#565051', state= 'disabled', command=lambda: self.ResetApp(), activebackground= '#565051', relief='flat', highlightthickness=0, bd=0)
        self.button8.place(relx=0.855, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button8, "Reset application", 40, 45)   

#-----------------------------------------------------------------------------#
#------------------------- R E A D   S H A P E F I L E -----------------------#
#-----------------------------------------------------------------------------#
    def readshp(self):
        self.filename = tk.filedialog.askopenfilename(initialdir= ".", filetypes=[("ESRI shape files", "*.shp")])
        
        if self.filename:
            try:
                del self.canvas  
            except:
                pass
            self.polygons = gpd.GeoDataFrame.from_file(self.filename)
            
            self.names = self.polygons.columns
            self.atr = list(self.names[1:-1])
            
#-------------------------- P L O T   S H A P E F I L E ----------------------#     
            self.f, self.ax = plt.subplots(figsize=(5,4), dpi=100)
            self.f.tight_layout(pad=0)
            self.ax.set_aspect('equal')
            
            self.ax = self.polygons.plot(ax= self.ax, color= '#FFE4B5', edgecolor= 'k', clip_on=False)
            self.f.set_facecolor("#98AFC7")
            self.ax.axis('off')
            plt.close()
            
            self.ax.fmt_xdata = lambda x: "{:.3f}".format(x)  
            self.ax.fmt_ydata = lambda x: "{:.3f}".format(x)
    
            self.canvas = FigureCanvasTkAgg(self.f, master=self.frame)
            self.canvas.get_tk_widget().place(relx= 0.15, rely= 0.19, relheight= 0.65, relwidth= 0.7)
            self.canvas.draw()
            
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.frame)
            for button in self.toolbar.winfo_children():
                button.config(background='#565051')
            self.toolbar.config(background='#565051')
            self.toolbar._message_label.place(relx=0.04, rely=0.02)
            self.toolbar.place(relx= 0.26, rely= 0.85, relheight=0.19)
            self.toolbar.update()

            self.button2.config(state= 'normal')
            self.button8.config(state= 'normal')
            return(self.atr)
        
#-----------------------------------------------------------------------------#        
#--- P L O T   C E N T R O I D S  &  C A L C U L A T E   D I S T A N C E S ---#
#-----------------------------------------------------------------------------# 
    def calcDistances(self):
        self.centroids = self.polygons.centroid
        #self.centroids.to_file("hdMoransIcentroids.shp")
        
        try:
            del self.canvas  
        except:
            pass
 
        self.f, self.ax = plt.subplots(figsize=(5,4), dpi=100)
        self.ax.set_aspect('equal')
        self.f.tight_layout(pad=0)
            
        self.ax = self.polygons.plot(ax= self.ax, color= '#FFE4B5', edgecolor= 'k', clip_on=False)
        self.centroids.plot(ax= self.ax, color= 'green', markersize=10)
        self.f.set_facecolor("#98AFC7")
        self.ax.axis('off')
        plt.close()
        
        self.ax.fmt_xdata = lambda x: "{:.3f}".format(x)  
        self.ax.fmt_ydata = lambda x: "{:.3f}".format(x)

        self.canvas = FigureCanvasTkAgg(self.f, master=self.frame)
        self.canvas.get_tk_widget().place(relx= 0.15, rely= 0.19, relheight= 0.65, relwidth= 0.7)
        self.canvas.draw()
        
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.frame)
        for button in self.toolbar.winfo_children():
            button.config(background='#565051')
        self.toolbar.config(background='#565051')
        self.toolbar._message_label.place(relx=0.04, rely=0.02)
        self.toolbar.place(relx= 0.26, rely= 0.85, relheight=0.19)
        self.toolbar.update()
       
        self.dists = np.array(self.centroids.geometry.apply(lambda g: self.centroids.geometry.distance(g)))
        
        #self.minDist = self.dists.mask(self.dists == 0).min().min()
        #self.maxDist = self.dists.mask(self.dists == 0).max().max() 
        
        self.minDist = np.min(self.dists[np.nonzero(self.dists)])
        self.maxDist = np.max(self.dists[np.nonzero(self.dists)])
#-----------------------------------------------------------------------------#
#--------------------------- F E A T U R E   M E N U -------------------------#
#-----------------------------------------------------------------------------#  
    def callwinFeat(self):
        self.winFeat = tk.Toplevel(self.master, bg= '#565051')
        self.winFeat.geometry("300x200+500+50")
        self.currentFeat = tk.StringVar()
        self.labelFeat = tk.Label(self.winFeat, text= "Select Attribute:", bg= '#565051').pack(pady=5)
        self.comboFeat = ttk.Combobox(self.winFeat, value= self.atr, textvariable= self.currentFeat, state="readonly")
        self.comboFeat.current(0)
        self.comboFeat.pack(pady=10)
        
        self.buttonAttrTable = tk.Button(self.winFeat, image=  self.imageAttrTable, activebackground= '#565051',command= lambda: self.showTable(), relief='flat', highlightthickness=0, bd=0)
        self.buttonAttrTable.pack(pady=15)
        tooltip.CreateToolTip(self.buttonAttrTable, "Show attribute table", 35, 30)
        
        self.buttonOkFeat = tk.Button(self.winFeat, image=  self.imageOk, activebackground= '#565051',command= lambda: self.Featclick(), bg= '#565051', relief='flat', highlightthickness=0, bd=0, width= 50)
        self.buttonOkFeat.pack(pady=10)
                    
    def Featclick(self):
        self.feat = self.currentFeat.get()
        self.ft = self.polygons[self.feat]
        self.n = len(self.ft)
        self.lfeat.set(f'{self.feat}')
        self.button3.config(state='normal')
        self.button4.config(state='normal')
        self.winFeat.destroy()
        return(self.feat)

#-----------------------------------------------------------------------------#
#---------------------- A T T R I B U T E   T A B L E ------------------------#
#-----------------------------------------------------------------------------#
    def attrTable(self):
        display = pd.options.display
        display.max_columns = 100
        display.max_rows = 100
        display.max_colwidth = 199
        display.width = None
        self.dataTable = self.polygons.loc[:, self.polygons.columns != 'geometry']
             
    def showTable(self):
        self.attrTable()
        
        self.winTable = tk.Toplevel(self.master, bg= '#565051')
        self.winTable.geometry("700x600+300+50")         
        self.winTable.title('Attribute Table')
        
        self.winFrame = tk.Frame(self.winTable, bg= '#565051')
        self.winFrame.place(relx= 0, rely=0, relheight= 1, relwidth= 1)
        
        self.pt = Table(self.winFrame, dataframe= self.dataTable, showtoolbar=True, showstatusbar=True)
        self.pt.show()  
    
#-----------------------------------------------------------------------------#        
#----------------------- O U T L I E R S   M E N U ---------------------------#
#-----------------------------------------------------------------------------#      
    def callWinOutliers(self):
        self.winOutliers = tk.Toplevel(self.master, bg= '#565051')
        self.winOutliers.geometry("300x150+500+50")
        
        self.labelOutliers = tk.Label(self.winOutliers, text= "Check for outliers:", bg= '#565051').pack(pady=10)
        self.buttonOutliers = tk.Button(self.winOutliers, image = self.imageOutliers2, activebackground= '#565051',command= lambda: self.checkOutliers(), bg= '#565051', relief='flat', highlightthickness=0, bd=0, width= 50)
        self.buttonOutliers.pack(pady=10)
        tooltip.CreateToolTip(self.buttonOutliers, "Check for outliers", 35, 30)
        
        self.buttonOkoutliers = tk.Button(self.winOutliers, image=  self.imageOk,activebackground= '#565051', command= lambda: self.winOutliers.destroy(), bg= '#565051', relief='flat', highlightthickness=0, bd=0, width= 50)
        self.buttonOkoutliers.pack(pady=10)
          
    def checkOutliers(self):       
        self.ft = pd.DataFrame(self.ft)        
        Q1 = self.ft.quantile(q=.25)
        Q3 = self.ft.quantile(q=.75)
        IQR = self.ft.apply(stats.iqr)

        self.data_clean = self.ft[~((self.ft < (Q1-1.5*IQR)) | (self.ft > (Q3+1.5*IQR))).any(axis=1)]
        self.ft = np.array(self.ft)
        self.data_clean = np.array(self.data_clean)
        self.outlist = np.setdiff1d(self.ft, self.data_clean)
        
        if np.size(self.outlist) > 0:   
            self.plotData()         
            self.responseOut = tk.messagebox.askyesno('Outliers:', f'{np.size(self.outlist)} outliers detected. Do you want to remove them?')
            
            if self.responseOut == 1:
                self.polygons = self.polygons[~self.polygons[self.feat].isin(self.outlist)]
                self.ft = self.polygons[self.feat]
                self.n = len(self.ft)
                self.plotData2()
            else:
                pass      
        else:
            self.plotData()
            tk.messagebox.showinfo('Outliers:', 'No outliers detected.')
 
#-----------------------------------------------------------------------------#
#-------------------- P L O T   W I T H   O U T L I E R S --------------------#
#-----------------------------------------------------------------------------#      
    def plotData(self): 
        self.indexes = np.arange(1,self.n+1)
        self.col = np.where(np.in1d(self.ft, self.data_clean), 'g', 'r')
                
        fig = plt.figure(figsize= (7,4), dpi=100)
        ax = fig.gca()
        ax.set_xticks(self.indexes)
        plt.title("Attribute plot",fontsize=20)
        plt.xlabel("Indexes",fontsize=14)
        plt.ylabel(f"{self.feat}",fontsize=14)      
        plt.axhline(y= self.ft.mean(), alpha=0.8, lw= 2, color='r')
                 
        plt.scatter(self.indexes, self.ft, color= self.col, marker='o', facecolors= 'none')
        plt.grid(True, axis='x', linestyle='--', alpha=0.3)       
        plt.show()

#-----------------------------------------------------------------------------#
#------------- P L O T   A F T E R   R E M O V I N G   O U T L I E R S -------#  
#-----------------------------------------------------------------------------# 
    def plotData2(self): 
        self.indexes = np.arange(1,self.n+1)
                
        fig = plt.figure(figsize= (7,4), dpi=100)
        ax = fig.gca()
        ax.set_xticks(self.indexes)
        plt.title("Attribute plot",fontsize=20)
        plt.xlabel("Indexes",fontsize=14)
        plt.ylabel(f"{self.feat}",fontsize=14)      
        plt.axhline(y= self.ft.mean(), alpha=0.8, lw= 2, color='r')
                 
        plt.scatter(self.indexes, self.ft, color= 'g', marker='o', facecolors= 'none')
        plt.grid(True, axis='x', linestyle='--', alpha=0.3)       
        plt.show()

#-----------------------------------------------------------------------------#        
#--------------------------- M E T H O D   M E N U ---------------------------#
#-----------------------------------------------------------------------------#  
    def callwinMethod(self):
        self.winMethod = tk.Toplevel(self.master, bg= '#565051')
        self.winMethod.geometry("300x150+500+50")
        self.currentMethod = tk.StringVar()
        self.labelMethod = tk.Label(self.winMethod, text= "Select Method:", bg= '#565051').pack(pady=10)
        self.methods = ['Neighbors', 'Distance', 'Inverse Distance']
        self.comboMethod = ttk.Combobox(self.winMethod, value= self.methods, textvariable= self.currentMethod, state="readonly")
        self.comboMethod.current(0)
        self.comboMethod.pack(pady=10)
        
        self.buttonOkMethod = tk.Button(self.winMethod, image=  self.imageOk,activebackground= '#565051', command= lambda: self.Methodclick(), bg= '#565051', relief='flat', highlightthickness=0, bd=0, width= 50)
        self.buttonOkMethod.pack(pady=10)
            
    def Methodclick(self):
        self.method = self.currentMethod.get()
        self.lfeat.set(f'{self.feat} / {self.method}')
        self.button5.config(state='normal')
        self.winMethod.destroy()
        
        if self.method == self.methods[1] or self.method == self.methods[2]:
            self.calcDistances()
        return(self.method)
    
#-----------------------------------------------------------------------------#   
#-------------------------- N E I G H B O R S   M E N U ----------------------#
#-----------------------------------------------------------------------------# 
    def callwinNeighbor(self):
        if self.method == self.methods[0]:
            self.winNeighbor = tk.Toplevel(self.master, bg= '#565051')
            self.winNeighbor.geometry("300x250+500+50")
            self.currentNeighbor = tk.StringVar()
            self.labelNeighMethod = tk.Label(self.winNeighbor, text= "Select Method:", bg= '#565051').pack(pady=5)
            self.neighbormethods = ['Rook', 'Queen']
            self.comboNeighbor= ttk.Combobox(self.winNeighbor, value= self.neighbormethods, textvariable= self.currentNeighbor, state="readonly")
            self.comboNeighbor.current(0)
            self.comboNeighbor.pack(pady=5)
            self.labelNeighNumber = tk.Label(self.winNeighbor, text= "Select Number of Neighbors:", bg= '#565051').pack(pady=5)
            self.neighSlider = tk.Scale(self.winNeighbor, from_= 1, to= self.n-1, orient= 'horizontal', bg= '#565051', relief='flat', bd=0)
            self.neighSlider.pack(pady=5)
            self.labelbuffer = tk.Label(self.winNeighbor, text= "Add buffer size:", bg= '#565051').pack(pady=2)
            self.entrybuffer = tk.Entry(self.winNeighbor, textvariable= self.buffer, bd=2, justify= 'c', width= 10)
            self.entrybuffer.pack(pady=2)
            
            self.buttonOkNeighbor = tk.Button(self.winNeighbor, image=  self.imageOk, activebackground= '#565051',command= lambda: self.Neighborclick(),bg= '#565051',  relief='flat', highlightthickness=0, bd=0, width= 50)
            self.buttonOkNeighbor.pack(pady=6)
            
        elif self.method == self.methods[1]:
            self.winDistance = tk.Toplevel(self.master, bg= '#565051')
            self.winDistance.geometry("300x200+500+50")
            self.currentDistance = tk.StringVar()
            self.labelDist = tk.Label(self.winDistance, text= "Select Radius Distance:", bg= '#565051').pack(pady=4)
            self.distanceSlider = tk.Scale(self.winDistance, from_= self.minDist, to= self.maxDist, orient= 'horizontal', bg= '#565051', relief='flat', bd=0)
            self.distanceSlider.pack(pady=4)
            self.labelNeighNumber = tk.Label(self.winDistance, text= "Select Number of Neighbors:", bg= '#565051').pack(pady=4)
            self.neighSlider = tk.Scale(self.winDistance, from_= 1, to= self.n-1, orient= 'horizontal', bg= '#565051', relief='flat', bd=0)
            self.neighSlider.pack(pady=4)
            
            self.buttonOkDistance = tk.Button(self.winDistance, image=  self.imageOk, activebackground= '#565051',command= lambda: self.Distclick(), bg= '#565051', relief='flat', highlightthickness=0, bd=0, width= 50)
            self.buttonOkDistance.pack(pady=4)
            
        elif self.method == self.methods[2]:
            self.winDistanceInv = tk.Toplevel(self.master, bg= '#565051')
            self.winDistanceInv.geometry("300x280+500+50")
            self.currentDistanceInv = tk.StringVar()
            self.labelDistInv = tk.Label(self.winDistanceInv, text= "Select Radius Distance:", bg= '#565051').pack(pady=4)
            self.distanceSlider = tk.Scale(self.winDistanceInv, from_= self.minDist, to= self.maxDist, orient= 'horizontal', bg= '#565051', relief='flat', bd=0)
            self.distanceSlider.pack(pady=4)
            self.labelpower = tk.Label(self.winDistanceInv, text= "Select Power:", bg= '#565051').pack(pady=4)
            self.powerSlider = tk.Scale(self.winDistanceInv, from_= 0, to= 10, orient= 'horizontal', bg= '#565051', relief='flat', bd=0)
            self.powerSlider.pack(pady=4)
            self.labelNeighNumber = tk.Label(self.winDistanceInv, text= "Select Number of Neighbors:", bg= '#565051').pack(pady=4)
            self.neighSlider = tk.Scale(self.winDistanceInv, from_= 1, to= self.n-1, orient= 'horizontal', bg= '#565051', relief='flat', bd=0)
            self.neighSlider.pack(pady=4)
            
            self.buttonOkDistanceInv = tk.Button(self.winDistanceInv, image= self.imageOk, activebackground= '#565051',command= lambda: self.DistInvclick(),bg= '#565051',  relief='flat', highlightthickness=0, bd=0, width= 50)
            self.buttonOkDistanceInv.pack(pady=4)         
    
    def Neighborclick(self):
        self.neighbormethod = self.currentNeighbor.get()
        self.neighN = self.neighSlider.get()
        self.bufferget = self.buffer.get()
        self.lfeat.set(f'{self.feat} / {self.method} / {self.neighbormethod} / N : {self.neighN} / B : {round(self.bufferget)}')
        self.button6.config(state='normal')
        self.winNeighbor.destroy()
        return(self.neighbormethod)
    
    def Distclick(self):
        self.testDist = self.distanceSlider.get()
        self.neighN = self.neighSlider.get()
        self.lfeat.set(f'{self.feat} / {self.method} / D: {round(self.testDist)} / N : {self.neighN}')
        self.button6.config(state='normal')
        self.winDistance.destroy()
        return(self.neighbormethod)
    
    def DistInvclick(self):
        self.testDist = self.distanceSlider.get()
        self.neighN = self.neighSlider.get()
        self.power = self.powerSlider.get()
        self.lfeat.set(f'{self.feat} / {self.method} / D: {round(self.testDist)} / N : {self.neighN} / p : {self.power}')
        self.button6.config(state='normal')
        self.winDistanceInv.destroy()
        return(self.neighbormethod)
    
#-----------------------------------------------------------------------------#   
#-------------------------- R E S U L T S   M E N U --------------------------#
#-----------------------------------------------------------------------------# 
    def callwinResults(self):  
        self.Permutations()
        self.winResults = tk.Toplevel(self.master, bg= '#565051')
        self.winResults.geometry("700x600+200+50")
        
        self.text = ''
        
        self.frame = tk.Frame(self.winResults, bg= '#565051')
        self.frame.place(relx= 0, rely=0, relheight= 1, relwidth= 1)
        
        self.screen = tk.Text(self.frame)
        self.screen.place(relx= 0.02, rely=0.02, relheight= 0.82, relwidth= 0.95)
        
        self.savebutton = tk.Button(self.frame, image= self.imageSave, activebackground= '#565051', bg= '#565051', relief='flat', highlightthickness=0, bd=0, command=lambda:self.savetxt())
        self.savebutton.place(relx = 0.68, rely= 0.88, relheight= 0.07, relwidth= 0.17)
        tooltip.CreateToolTip(self.savebutton, "Save results", 40, 45)
    
        self.closebutton = tk.Button(self.frame, image= self.imageBack, activebackground= '#565051', bg= '#565051', relief='flat', highlightthickness=0, bd=0, command=lambda:self.winResults.destroy())
        self.closebutton.place(relx = 0.85, rely= 0.88, relheight= 0.07, relwidth= 0.12)
        tooltip.CreateToolTip(self.closebutton, "Back", 40, 45)
    
        self.scrollbar2 = tk.Scrollbar(self.screen, jump= 1, orient= 'vertical', command= self.screen.yview)
        self.scrollbar2.pack(side= tk.RIGHT, fill= tk.Y)
         
        self.scrollbar3 = tk.Scrollbar(self.screen, jump= 1, orient= 'horizontal', command= self.screen.xview)
        self.scrollbar3.pack(side= tk.BOTTOM, fill= tk.X)
        
        self.screen.configure(wrap= 'none', yscrollcommand=self.scrollbar2.set, xscrollcommand=self.scrollbar3.set) 
        
        self.Cs = pd.DataFrame(self.Cs)
        self.Ws = pd.DataFrame(self.Ws)
        self.devp = pd.DataFrame(self.devp)
        
        if self.method == self.methods[0]:
            self.text = f"-Feature: {self.feat}\n-Method: {self.method} - {self.neighbormethod}\n-Number of neighbors: {self.neighN}\n-Buffer: {self.bufferget}\n-Moran's I: {round(self.MoransI, 5)}\n-E(I): {round(self.EI,5)}\n-μ: {round(self.MIm, 5)}\n-σ: {round(self.MIstd, 5)}\n-z: {round(self.Zscore, 3)}\n-P value: 0.001\n-Weight Table:\n{np.round(self.Cs,decimals=1)}\n-Normalized Weight Table:\n{np.round(self.Ws,decimals=1)}\n-Deviation Table:\n{np.round(self.devp,decimals=1)}\n"                    
        elif self.method == self.methods[1]:
            self.text = f"-Feature: {self.feat}\n-Method: {self.method} - {round(self.testDist)}m\n-Number of neighbors: {self.neighN}\n-Moran's I: {round(self.MoransI, 5)}\n-E(I): {round(self.EI,5)}\n-μ: {round(self.MIm, 5)}\n-σ: {round(self.MIstd, 5)}\n-z: {round(self.Zscore, 3)}\n-P value: 0.001\n-Weight Table:\n{np.round(self.Cs,decimals=1)}\n-Normalized Weight Table:\n{np.round(self.Ws,decimals=1)}\n-Deviation Table:\n{np.round(self.devp,decimals=1)}\n"                           
        elif self.method == self.methods[2]:
            self.text = f"-Feature: {self.feat}\n-Method: {self.method} - {round(self.testDist)}m\n-Power: {self.power}\n-Number of neighbors: {self.neighN}\n-Moran's I: {round(self.MoransI, 5)}\n-E(I): {round(self.EI,5)}\n-μ: {round(self.MIm, 5)}\n-σ: {round(self.MIstd, 5)}\n-z: {round(self.Zscore, 3)}\n-P value: 0.001\n-Weight Table:\n{np.round(self.Cs,decimals=1)}\n-Normalized Weight Table:\n{np.round(self.Ws,decimals=1)}\n-Deviation Table:\n{np.round(self.devp,decimals=1)}\n"                    
        
        self.screen.insert(tk.INSERT, self.text)
        self.screen.config(state= 'disabled') 
    
    def savetxt(self): 
        self.file_name = tk.filedialog.asksaveasfilename(parent= self.master)
        if self.file_name:  
            self.pos = f'-Path: {self.file_name} \n'
            with open(f'{self.file_name}.txt', 'w') as file:
                file.write(self.pos)
                file.write(self.text)
            self.Cs.to_csv(f'{self.file_name}W.csv', header=None, index=None, sep=',', mode='a')
            self.Ws.to_csv(f'{self.file_name}NW.csv', header=None, index=None, sep=',', mode='a')
            self.devp.to_csv(f'{self.file_name}D.csv', header=None, index=None, sep=',', mode='a')
            tk.messagebox.showinfo("Moran's I:", 'Result files saved successfully!')
            
#-----------------------------------------------------------------------------#
#------------------ C A L C U L A T E   M O R A N S  I------------------------#
#-----------------------------------------------------------------------------#
    def calcMoransI(self):
        if self.method == self.methods[0]:
            self.NeighborMoransI()
            self.MoransIscatterPlot()
            self.MIreps = [self.MoransI]
        
        elif self.method == self.methods[1]:
            self.DistanceMoransI()
            self.MoransIscatterPlot()
            self.MIreps = [self.MoransI]
    
        elif self.method == self.methods[2]:
            self.IDWMoransI()
            self.MoransIscatterPlot()
            self.MIreps = [self.MoransI]
        
        self.button7.config(state= 'normal')

#-----------------------------------------------------------------------------#
#---------------------- N E I G H B O R   P O L Y G O N S --------------------#
#-----------------------------------------------------------------------------#
    def NeighborMoransI(self): 
        self.ft = np.array(self.polygons[self.feat])
        self.polygons['geometry']= self.polygons.buffer(self.bufferget)
            
        #Contiguity table      
        self.xneigh = np.arange(len(self.polygons))
        if self.neighbormethod == 'Rook':
            self.polygs = np.array(self.polygons.geometry.apply(lambda g: self.polygons.geometry.intersects(g)))
        else:
            self.polygs = np.array(self.polygons.geometry.apply(lambda g: self.polygons.geometry.touches(g)))
        self.polygs = self.polygs.astype(int)
        np.fill_diagonal(self.polygs, 0)
        self.Cs = np.where( (self.polygs[self.xneigh,:] == 1) & (np.cumsum(self.polygs, axis=1)<= self.neighN), 1, 0)
        
        #Weight Table
        self.Rs = self.Cs.sum(axis=1)[:,None]
        
        #self.Ws = np.divide(self.Cs, self.Rs, out=np.zeros_like(self.Cs), where=self.Rs!=0)
  
        with np.errstate(divide='ignore', invalid='ignore'):
          self.Ws = np.true_divide(self.Cs, self.Rs)
          self.Ws[self.Ws == np.inf] = 0
          self.Ws = np.nan_to_num(self.Ws)
          
        self.MIcalc()

#-----------------------------------------------------------------------------#
#---------------------- C E N T R O I D S   D I S T A N C E ------------------#
#-----------------------------------------------------------------------------#
    def DistanceMoransI(self):     
        #centroids.to_file("centroids.shp")
        self.ft = np.array(self.polygons[self.feat])
            
        #Distances between neighbors      
        self.mins = np.argsort(self.dists, 1)[:, 0:self.neighN + 1]
     
        #Weight calculation
        self.mindists = np.take_along_axis(self.dists, self.mins, 1)
    
        self.xrange = np.arange(self.n)     
        self.Cs = np.where( (self.dists > 0) & (np.isin(self.dists[self.xrange,:], self.mindists[self.xrange,:])), 1, 0 )
        self.Rs = self.Cs.sum(axis=1)[:,None]

        #Weight Table
        #self.Cs = np.array(self.Cs).reshape((self.n, self.n))       
        #Normalize Weight Table
        #self.Ws = self.Cs / self.Cs.sum(axis=1)[:,None]            
        #self.Ws = np.divide(self.Cs, self.Rs, out=np.zeros_like(self.Cs), where=self.Rs!=0)
 
        #Normalize weight table
        with np.errstate(divide='ignore', invalid='ignore'):
           self.Ws = np.true_divide(self.Cs, self.Rs)
           self.Ws[self.Ws == np.inf] = 0
           self.Ws = np.nan_to_num(self.Ws)
           
        self.MIcalc()
            
#-----------------------------------------------------------------------------#
#-------------------- I N V E R S E   D I S T A N C E ------------------------#
#-----------------------------------------------------------------------------#
    def IDWMoransI(self):
        #centroids.to_file("centroids.shp")
        self.ft = np.array(self.polygons[self.feat])
        
        #Distances between neighbors      
        self.mins = np.argsort(self.dists, 1)[:, 0:self.neighN + 1]          
        
        #Weight calculation
        self.mindists = np.take_along_axis(self.dists, self.mins, 1)
        self.xrange = np.arange(self.n)
        
        self.idp = np.divide(1, self.dists**self.power, where=self.dists!=0)  
        self.Cs = np.where( (self.dists > 0) & (np.isin(self.dists[self.xrange,:], self.mindists[self.xrange,:])), self.idp, 0 )
       
        self.Rs = self.Cs.sum(axis=1)[:,None]
        
        #Weight Table
        #self.Cs = np.array(self.Cs).reshape((self.n, self.n))
        #self.Ws = np.divide(self.Cs, self.Rs, out=np.zeros_like(self.Cs), where=self.Rs!=0)
        
        #Normalize weight table
        with np.errstate(divide='ignore', invalid='ignore'):
           self.Ws = np.true_divide(self.Cs, self.Rs)
           self.Ws[self.Ws == np.inf] = 0
           self.Ws = np.nan_to_num(self.Ws)
        
        self.MIcalc()       
    
    def MIcalc(self):
        #Residuals calculation
        self.dev = self.ft - np.mean(self.ft)
        
        #Square residuals
        self.dev2 = self.dev**2
        
        #Sum of square residuals
        self.Sdev2 = np.sum(self.dev2)
         
        #Cartesian product of residuals
        deva, devb = np.meshgrid(self.dev, self.dev)
        self.devp = deva * devb
        self.wft = self.Ws * self.devp
        
        #Sum of cartesian product
        self.Swtf = np.sum(self.wft)
        
        #Moran's index 
        self.MoransI = self.Swtf / self.Sdev2
        #self.MoransI = round(self.MoransI, 5)
        
#-----------------------------------------------------------------------------#
#------------------ M O R A N S   I   S C A T T E R   P L O T ----------------#
#-----------------------------------------------------------------------------#
    def MoransIscatterPlot(self):
        self.zft = self.dev / np.std(self.ft)
        self.Zw = np.sum(self.Ws * self.zft, axis=1)  
        
        minzs2 = np.min(self.zft)
        maxzs2 = np.max(self.zft)
        minwzs2 = np.min(self.Zw)
        maxwzs2 = np.max(self.Zw)
        axislim = np.max(np.absolute([minzs2, maxzs2, minwzs2, maxwzs2]))
     
        plt.xlim(-1.05 * axislim, 1.05 * axislim)
        plt.ylim(-1.05 * axislim, 1.05 * axislim)
        plt.gca().set_aspect('equal', adjustable='box')     
        plt.scatter(self.zft, self.Zw, color='b', marker='o', facecolors= 'none')
        plt.title(f"Moran's I : {round(self.MoransI, 5)}",fontsize=20)
        plt.xlabel(f"{self.feat}",fontsize=14)
        plt.ylabel(f"lag {self.feat}",fontsize=14)
        plt.tick_params(axis="both",labelsize=10)
        plt.plot(np.unique(self.dev), np.poly1d(np.polyfit(self.zft, self.Zw, 1))(np.unique(self.dev)), 'r', lw= 2)
    
        plt.axvline(x= 0, ls='--', alpha=0.6, lw= 1, color= 'k')
        plt.axhline(y= 0, ls='--', alpha=0.6, lw= 1, color='k')
        plt.show()
        
#-----------------------------------------------------------------------------#    
#------------------------- P E R M U T A T I O N S ---------------------------#      
#-----------------------------------------------------------------------------#  
    def Permutations(self):  
        self.originalMI = self.MoransI
        for i in range(999):
            self.ft = np.random.permutation(self.ft)       
            self.MIcalc()
            self.MIreps.append(self.MoransI)
           
        self.MIreps = np.array(self.MIreps)    
        self.MIm = np.mean(self.MIreps)
        self.MIstd = np.std(self.MIreps)
        self.EI = -1/(len(self.ft) - 1)
        self.MIreps = np.round(self.MIreps, decimals=5)
        
        self.Zscore = (self.originalMI - self.EI) / self.MIstd
        self.pvalue = 0.001
        
        fig = plt.figure(figsize= (7,4), dpi=100)
        ax = fig.gca()
   
        plt.axvline(x= self.MIreps[0], lw= 2, color='r')
        plt.hist(self.MIreps, density= True, color= '#438D80', alpha= 0.8)
        xt = plt.xticks()[0]  
        xmin, xmax = min(xt), max(xt)  
        lnspc = np.linspace(xmin, xmax, len(self.MIreps))
        #m, s = stats.norm.fit(self.MIreps)
        pdf_g = stats.norm.pdf(lnspc, self.MIm, self.MIstd)
        plt.plot(lnspc, pdf_g, label="Norm", color= '#CD7F32')            
        
        plt.xlabel("Moran's I")
        plt.ylabel('Frequencies')
        plt.title(r'$\mathrm{Histogram\ of\ Permutations:}\ \mu=%.3f,\ \sigma=%.3f$' %(self.MIm, self.MIstd))
        
        plt.show()

#-----------------------------------------------------------------------------#
#--------------------- R E S E T   A P P L I C A T I O N ---------------------#
#-----------------------------------------------------------------------------#
    def ResetApp(self):
        self.frame.destroy()
        self.frame = tk.Frame(self.master, bg= '#565051')
        self.frame.place(relheight=1, relwidth=1)
        
        self.lfeat = tk.StringVar(value='')
        self.labelfeat = tk.Label(self.frame, textvariable= self.lfeat, bg= '#565051')
        self.labelfeat.place(relx=0, rely=0.12, relwidth= 1, relheight= 0.05)
        self.feat = ''
        self.method = ''
        self.neighbormethod = ''
        self.neighN = 0
        self.buffer = tk.DoubleVar(value= 0)
        self.Dist = 0
        self.power = 0
        self.MoransI = 0
        self.outlcheck = False
        self.perm = 0
        
        self.button = tk.Button(self.frame, image= self.imageFolder, command=lambda: self.readshp(), bg= '#565051', activebackground= '#565051', relief='flat', highlightthickness=0, bd=0)
        self.button.place(relx=0.015, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button, "Open file", 40, 45)
        
        self.button2 = tk.Button(self.frame, image= self.imageAttrTable, command=lambda: self.callwinFeat(),  state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button2.place(relx=0.135, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button2, "Select attribute", 40, 45)
        
        self.button3 = tk.Button(self.frame, image= self.imageOutliers, command=lambda: self.callWinOutliers(), state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button3.place(relx=0.255, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button3, "Check for outliers", 40, 45)
        
        self.button4 = tk.Button(self.frame, image= self.imageMethod,  command=lambda: self.callwinMethod(), state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button4.place(relx=0.375, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button4, "Select method", 40, 45)
        
        self.button5 = tk.Button(self.frame, image= self.imageDistance,  command=lambda: self.callwinNeighbor(), state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button5.place(relx=0.495, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button5, "Select number of neighbors - distance", 40, 45)
       
        self.button6 = tk.Button(self.frame, image= self.imageCalc, command=lambda: self.calcMoransI(), state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button6.place(relx=0.615, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button6, "Calculate Moran's I", 40, 45)
        
        self.button7 = tk.Button(self.frame, image= self.imageResults,  command=lambda: self.callwinResults(), state= 'disabled', bg= '#565051', activebackground= '#565051',relief='flat', highlightthickness=0, bd=0)
        self.button7.place(relx=0.735, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button7, "Show results", 40, 45)
        
        self.button8 = tk.Button(self.frame, image= self.imageRefresh, state= 'disabled', bg= '#565051', command=lambda: self.ResetApp(), activebackground= '#565051', relief='flat', highlightthickness=0, bd=0)
        self.button8.place(relx=0.855, rely=0, relwidth= 0.12, relheight= 0.1)
        tooltip.CreateToolTip(self.button8, "Reset application", 40, 45)   
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#

def main():
    root = tk.Tk()
    app = MainWindow(root)
    root.mainloop()

if __name__ == '__main__':
    main()
    
#********************************** E N D ************************************#