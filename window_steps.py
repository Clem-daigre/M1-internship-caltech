import numpy as np
import matplotlib.pyplot as plt
from LOO_Utilitie_functions import *

#########################           UTILITY FUNCTIONS        ###################################
def evaluate_beta(startpt,cp,endpt,splitrate) :
    DT1=len(splitrate[startpt:cp])*20
    DT2=len(splitrate[cp:endpt])*20
    N1=np.sum(splitrate[startpt:cp])
    N2=np.sum(splitrate[cp:endpt])
    D=np.abs(N1*DT2/DT1)
    return round((N2-D)/np.sqrt(D),2)

def deltabeta(startpt,cp,endpt,splitrate) :
    Dt1=len(splitrate[startpt:cp])*20
    Dt2=len(splitrate[cp:endpt])*20
    N1=np.sum(splitrate[startpt:cp])
    N2=np.sum(splitrate[cp:endpt])
    c=scipy.stats.chi2.isf(0.975, 2*N1, loc=0, scale=1)
    d=scipy.stats.chi2.isf(0.025, (2*N1+2), loc=0, scale=1)
    lbdainf=c/(2*Dt1)
    lbdasup=d/(2*Dt1)
    betainf=(N2-Dt2*lbdainf)/np.sqrt(Dt2*lbdainf)
    betasup=(N2-Dt2*lbdasup)/np.sqrt(Dt2*lbdasup)
    return np.abs(betasup-betainf)/2

#########################           CLASS       ########################################
class signteststeps() :
    """ Automatically determine the location of possible change points and their statistical significance. Print and plot all intermediate results."""

    def __init__ (self,infile,raterange,R) :
        ### initialization ###

        self.infile=infile
        self.raterange=raterange
        self.R=R # number of fake catalogs simulated#4000
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
        return

    def proba(self) :
        ### compute overall PDF ###

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
        print('overall PDF',proba)
        return

    def cpid(self) :
        ### identify pics in the overall PDF ###

        cp=[0]# cp will be the change point indexes list

        #keep maximums over 10 points intervals if its pdf is above 0.02
        #first interval
        m=np.argmax(self.proba[1:10])+1 #cp can't be 0
        if self.proba[m]>=0.02 :
            cp.append(m)
        #iterate 10 point windows
        i=1
        while (i+1)*10<len(self.proba) :
            m=np.argmax(self.proba[i*10:(i+1)*10])+i*10 #identify a maximum
            if self.proba[m]>=0.02 :#test it is a true pic
                cp.append(m)
            i+=1
        #last interval
        m=np.argmax(self.proba[i*10:])+i*10
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
        #plot
        cum=[i+1 for i in range(len(self.time))]
        fig, ax1 = plt.subplots()
        ax1.plot(self.time, cum, 'navy', linewidth=1)
        ax1.set_xlabel('Date',size=12)
        locs=[0,1000,2000,3000,4000,4800]
        labels=['2005-1','2007-9','2010-6','2013-3','2015-12','2018-2']
        plt.xticks(locs,labels,fontsize=12)
        ax1.set_ylabel('Cumulative nomber', color='navy',size=12)
        for tl in ax1.get_yticklabels():
            tl.set_color('navy')
        ax2 = ax1.twinx()
        ax2.plot(self.splitime,self.proba,color="#8F2727")
        ax2.set_ylim([0, 0.7])
        ax2.set_ylabel('Probability Density Function',color="#8F2727",size=12)
        for tl in ax2.get_yticklabels():
            tl.set_color("#8F2727")
        for i in range(len(cp)) :
            ax2.plot([self.splitime[cp[i]],self.splitime[cp[i]]],[0,0.7],color='k',linewidth=0.3)
        plt.show()
        print('index of pics identified', self.cp)
        return

    def sign(self) :
        ### test the significiance of the pics ###

        cp=self.cp
        i=0
        while i+2!=len(cp) :
            count_data=self.splitrate[cp[i]:cp[i+2]]
            diff,var,cpfound=lpd(count_data,400) # lpd value of the change from last pic to next pic
            if diff-var>0 : # if change significant
                cp[i+1]=cp[i]+cpfound # keep the change point location found
                i+=1
            elif diff-var<=0 : # if change insignificant
                cp.pop(i+1) # delete pic
        plt.figure()
        cum=[i+1 for i in range(len(self.time))]
        plt.plot(self.time, cum, 'navy', linewidth=1)
        plt.xlabel('Date',size=12)
        locs=[0,1000,2000,3000,4000,4800]
        labels=['2005-1','2007-9','2010-6','2013-3','2015-12','2018-2']
        plt.xticks(locs,labels,fontsize=12)
        plt.ylabel('Cumulative number', color='navy',size=12)
        for i in range(len(cp)) :
            plt.plot([self.splitime[cp[i]],self.splitime[cp[i]]],[0,len(cum)],color='k',linewidth=0.3)
        plt.show()
        self.cp=cp
        print('significant change points first round', cp)
        return

    def signbis(self) :
        ### second check and determination if one sigma or two sigma certitude + result figure ###

        cp=self.cp
        self.sign1=[]
        self.sign2=[]
        diffs=[]
        stds=[]
        unsign=[]
        cpnew=[0]
        i=0
        while i+2!=len(cp) :
            count_data=self.splitrate[cp[i]:cp[i+2]]
            diff,var,cpfound=lpd(count_data,400)

            if diff-2*var>0 :#95 pourcent sure
                self.sign2.append(cp[i]+cpfound)
                cpnew.append(cp[i]+cpfound)
                cp[i+1]=cp[i]+cpfound
                i+=1
                diffs.append(diff)
                stds.append(var)
            elif diff-var>0 :#68 pourcent sure
                self.sign1.append(cp[i]+cpfound)
                cpnew.append(cp[i]+cpfound)
                cp[i+1]=cp[i]+cpfound
                i+=1
                diffs.append(diff)
                stds.append(var)
            else : 
                unsign.append(cp[i]+cpfound)
                cp.pop(i+1)
        cpnew.append(len(self.splitrate))

        if unsign!=[] : #if change point deleted through this process need to recheck
            self.sign1=[]
            self.sign2=[]
            cp=cpnew
            cpnew=[0]
            diffs=[]
            stds=[]
            i=0
            while i+2!=len(cp) :
                count_data=self.splitrate[cp[i]:cp[i+2]]
                diff,var,cpfound=lpd(count_data,400)
                
                if diff-2*var>0 :#95 pourcent sure
                    self.sign2.append(cp[i]+cpfound)
                    cpnew.append(cp[i]+cpfound)
                    cp[i+1]=cp[i]+cpfound
                    i+=1
                    diffs.append(diff)
                    stds.append(var)
                elif diff-var>0 :#68 pourcent sure
                    self.sign1.append(cp[i]+cpfound)
                    cpnew.append(cp[i]+cpfound)
                    cp[i+1]=cp[i]+cpfound
                    i+=1
                    diffs.append(diff)
                    stds.append(var)
                else : 
                    cp.pop(i+1)
        beta=[]
        dbeta=[]
        i=0
        while i+2 !=len(cpnew) :
            beta.append(evaluate_beta(cpnew[i],cpnew[i+1],cpnew[i+2],self.splitrate))
            dbeta.append(deltabeta(cpnew[i],cpnew[i+1],cpnew[i+2],self.splitrate))
            i+=1
        print('lists of diff in ldp values and its standard deviation of the change points in order',diffs,stds)
        print('lists of location of one sigma, 2 sigma change points',self.sign1,self.sign2)
        print('lists of beta and delta beta values of the change points',beta,dbeta)
        
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
        plt.show()
        return 


    #def checkcp(self) :
        #count_data=self.splitrate[178:235]
          #  print(cp[i],cp[i+2],count_data)
        #diff,var=lpdbis(count_data,400)
        #lpd(count_data,400,self.splitime
       # cp=cpsearch(count_data,self.R)
       # evaluate_beta(131,214,236,self.splitrate)
       # print(diff,var,cp)
       # return