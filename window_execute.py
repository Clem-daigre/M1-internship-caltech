import numpy as np
import matplotlib.pyplot as plt
from LOO_Utilitie_functions import *

class signtest() :
    """ Automatically determine the location of possible change points and their statistical significance """

    def __init__ (self,infile,raterange,R) :
        ### initialization ###

        self.infile=infile
        self.raterange=raterange
        self.R=R # number of fake catalogs simulated #4000
        return
    
    def read(self) :
        ### Read data and compute local seismicity rate data set ###

        file=open(self.infile,'r')
        line=file.readline()
        self.depth,self.lat,self.lon,self.magn,self.time=[[float(x) for x in s.split(',')] for s in line.split ('][')]
        file.close() 
        splitime=eqrate(self.raterange,self.time)[0] 
        self.splitime=splitime
        splitrate=eqrate(self.raterange,self.time)[1]
        self.splitrate=splitrate
        print(self.time[-1])
        return

    def proba(self) :
        ### compute overall PDF and identify pics ###

        proba=np.zeros((len(self.splitrate))) # probability density function since 2005
        count=np.zeros((len(self.splitrate))) # number of PDF values summed for each data point

        length=10 # shifting of 10 points
        if len(self.splitrate)%length==0 :
            delimiters=[i*length for i in range(len(self.splitrate)//length+1)] # keep indexes every 10 points
        else :
            delimiters=[i*length for i in range(len(self.splitrate)//length)]+[len(self.splitrate)] # add last data point if not divisble by 10
        
        for j in range(0,len(delimiters)-4) :
            count_data = self.splitrate[delimiters[j]:delimiters[j+4]] # window of 40 points
            pblocal = cpsearch(count_data,self.R) # pdf over the window
            proba[delimiters[j]:delimiters[j+4]] += pblocal # add it to overall probability
            count[delimiters[j]:delimiters[j+4]]+=1 # add it to counter
        proba=[proba[i]/count[i] for i in range (len(proba))] # average the pdf values
        self.proba=proba # keep overall pdf
        
        cp=[0] # cp will be the change point indexes list

        #keep maximums over 10 points intervals if its pdf is above 0.02
        m=np.argmax(self.proba[1:10])+1 #first interval #cp can't be 0
        if self.proba[m]>=0.02 :
            cp.append(m)
        #iterate 10 point windows
        i=1
        while (i+1)*10<len(self.proba) :
            m=np.argmax(self.proba[i*10:(i+1)*10])+i*10 #identify a maximum
            if self.proba[m]>=0.02 :#test it is a true pic
                cp.append(m)
            i+=1
        m=np.argmax(self.proba[i*10:])+i*10 #last interval
        if self.proba[m]>=0.02 :
           cp.append(m)
        
        #if 2 cp are too close, we keep the maximum of them
        j=0
        while j <len(cp)-1 :
            if np.abs(cp[j]-cp[j+1])<=5 : # if more close than 5 points
                if max(self.proba[cp[j]],self.proba[cp[j+1]])==self.proba[cp[j]] :#keep max
                    cp.pop(j+1)
                else :
                    cp.pop(j)
            j+=1
        cp.append(len(self.splitime)-1)#add last data point as pic value
        self.cp=cp
        return

    def sign(self) :
        ### test the significance of the PDF pics, keep final change points set and return values  ###

        cp=self.cp
        i=0 # first test of pics significance
        while i+2!=len(cp) :
            count_data=self.splitrate[cp[i]:cp[i+2]]
            diff,var,cpfound=lpd(count_data,400) # lpd value of the change from last pic to next pic
            if diff-var>0 : # if change significant
                cp[i+1]=cp[i]+cpfound # keep the change point location found
                i+=1
            elif diff-var<=0 : # if change insignificant
                cp.pop(i+1) # delete pic
        
        #second check and determination if one sigma or two sigma certitude
        self.sign1=[] #one sigma
        self.sign2=[] # two sigma
        unsign=[]
        cpnew=[0] # new list of cp after the check
        i=0
        while i+2!=len(cp) :
            count_data=self.splitrate[cp[i]:cp[i+2]]
            diff,var,cpfound=lpd(count_data,400)
            if diff-2*var>0 :# 95 pourcent sure
                self.sign2.append(cp[i]+cpfound)
                cpnew.append(cp[i]+cpfound)#keep location
                cp[i+1]=cp[i]+cpfound#change change pint i+1 value
                i+=1
            elif diff-var>0 :# 68 pourcent sure
                self.sign1.append(cp[i]+cpfound)
                cpnew.append(cp[i]+cpfound)
                cp[i+1]=cp[i]+cpfound
                i+=1
            else : 
                unsign.append(cp[i]+cpfound)
                cp.pop(i+1) #delete change point i+1

        cpnew.append(len(self.splitrate))

        if unsign!=[] : # if change point deleted through this process need to recheck
            self.sign1=[]# one sigma
            self.sign2=[] # two sigma
            cp=cpnew 
            cpnew=[0]# new list of cp after the check
            i=0
            while i+2!=len(cp) :
                count_data=self.splitrate[cp[i]:cp[i+2]]
                diff,var,cpfound=lpd(count_data,400)
                
                if diff-2*var>0 :#95 pourcent sure
                    self.sign2.append(cp[i]+cpfound)
                    cpnew.append(cp[i]+cpfound)
                    cp[i+1]=cp[i]+cpfound
                    i+=1
                elif diff-var>0 :#68 pourcent sure
                    self.sign1.append(cp[i]+cpfound)
                    cpnew.append(cp[i]+cpfound)
                    cp[i+1]=cp[i]+cpfound
                    i+=1
                else : 
                    cp.pop(i+1)

        cpvalues=[0] # convert rate point index to time value for further analysis
        for c in self.sign1 :
            cpvalues.append((c+1)*self.raterange-self.raterange/2) # include pointnumber 0 # middle of the time interval
        for c in self.sign2 :
            cpvalues.append((c+1)*self.raterange-self.raterange/2) # include pointnumber 0 # middle of the time interval
        cpvalues.append(self.time[-1])
        print('time of occurence of the changes in days since 2005', sorted(cpvalues))
        self.cpvalues=sorted(cpvalues)
        return
            
    def figure(self) :
        ### plot results figure ###

        plt.figure()
        cum=[i+1 for i in range(len(self.time))]
        plt.plot(self.time, cum, 'navy', linewidth=1.5)
        plt.xlabel('Time',size=15)
        locs=[0,1000,2000,3000,4000,4800]
        labels=['2005-1','2007-9','2010-6','2013-3','2015-12','2018-2']
        plt.xticks(locs,labels,fontsize=12)
        plt.ylabel('Cumulative number of events',size=15)
        plt.yticks(fontsize=12)
        i=0
        if len(self.sign2)!=0 :
            plt.plot([self.splitime[self.sign2[i]],self.splitime[self.sign2[i]]],[0,len(cum)],color='crimson',linewidth=1,label='95%')
            for i in range(1,len(self.sign2)) :
                plt.plot([self.splitime[self.sign2[i]],self.splitime[self.sign2[i]]],[0,len(cum)],color='crimson',linewidth=1)
        if len(self.sign1)!=0 :
            plt.plot([self.splitime[self.sign1[i]],self.splitime[self.sign1[i]]],[0,len(cum)],color='k',linewidth=1,label='68%')
            for i in range(1,len(self.sign1)) :
                plt.plot([self.splitime[self.sign1[i]],self.splitime[self.sign1[i]]],[0,len(cum)],color='k',linewidth=1)
        if len(self.sign2)!=0 or len(self.sign1)!=0 :
            plt.legend(bbox_to_anchor=(-0.05, 1.05), loc=2, borderaxespad=0.,fontsize=12)
        plt.savefig('./RESULTS/'+str(self.infile)[18:-4]+str(self.cpvalues)+'.png')
        plt.show()
        return