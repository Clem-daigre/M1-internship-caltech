import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.interpolate import griddata
from obspy.geodetics import base
import scipy
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

#########################           UTILITY FUNCTIONS        ###################################
def getKey(item) :
    return item[3]

def evaluate_beta(changetime, cluster,Dt2,Dt1) :
    ### Compute beta value and its 95 % confidence interval lower and upper values###

    if len(cluster)==0 :#cluster empty
        return 10,10,10

    if cluster[-1]<=changetime :#no events after the change time
        return 10,10,10

    else :
        N1=0
        i=0
        while cluster[i] < changetime :
            N1+=1
            i+=1
        N2 = len(cluster[i+1:])
        D = np.abs(N1*Dt2/Dt1)
        beta=(N2-D) / np.sqrt(D)
        c=scipy.stats.chi2.isf(0.975, 2*N1, loc=0, scale=1) # qui square distribution degree 2*N1
        d=scipy.stats.chi2.isf(0.025, (2*N1+2), loc=0, scale=1) # qui square distribution degree 2*N1+2
        lbdainf=c/(2*Dt1)
        lbdasup=d/(2*Dt1)
        betasup=(N2-Dt2*lbdainf)/np.sqrt(Dt2*lbdainf)
        betainf=(N2-Dt2*lbdasup)/np.sqrt(Dt2*lbdasup)
        return beta,betainf,betasup

#########################           CLASS       ########################################

class map() :
    """ Draw a map of the significance of the change over Tonga zone, the change point being known. """

    def __init__ (self,infile,changetime,minnb,spacemin,lonrange,latrange,spacemax,latbottom,lattop,lonleft,lonright) :
        #### Initialisation ###

        self.infile=infile #data file
        self.changetime=changetime # known change point
        self.spacemin=spacemin # in meters, minimum radius to select the events
        self.lonrange=lonrange #in km spacing between grid points
        self.latrange=latrange #in km spacing between grid points
        self.spacemax=spacemax # in meters, maximum radius to select the events
        self.latbottom=latbottom # map edges
        self.lattop=lattop
        self.lonleft=lonleft
        self.lonright=lonright
        self.nbmin=minnb # minimum number of events to compute beta
        return

    def priliminary (self) :
        ### Read data and create a regular grid of points coordinates  ###
        
        file=open(self.infile,'r')# read infile
        line=file.readline()
        self.depth,self.lat,self.lon,self.magn,self.time=[[float(x) for x in s.split(',')] for s in line.split ('][')]
        file.close()

        self.Dt2=self.time[-1]-self.changetime+1 # number of days in the second period
        self.Dt1=self.changetime # number of days in the first period
        
        if self.lonleft > 0 : # transform all positive longitude values to negative ones
            self.lonleft = -360+self.lonleft # 174 E becomes -186 W
        for i in range (len(self.lon)) :
            if self.lon[i] > 0 :
                self.lon[i] = -360+self.lon[i]
        
        self.x=[] # longitudes of points
        self.y=[] # latitudes of points
        latcurr = self.lattop - self.latrange/111
        while latcurr > self.latbottom :
            loncurr = self.lonleft + self.lonrange/(111) #start at left top corner
            while loncurr < self.lonright : #all lons
                self.x.append (loncurr)
                self.y.append (latcurr)
                loncurr = loncurr + self.lonrange/(111)
            latcurr = latcurr - self.latrange/111
        print(self.x,self.y)
        return

    def compute(self) :
        ### Compute change status value for each point of the grid ###

        BETAS=[] # list of status values for the points of the grid
        for j in range(len(self.x)) :# iterate the grid

            info=[] # compute dist of events from current point
            for i in range(len(self.lon)) :
                info.append([self.lat[i],self.lon[i],self.time[i],base.gps2dist_azimuth(self.y[j],self.x[j],self.lat[i],self.lon[i])[0]])
            info=sorted(info,key=getKey) # sort events with increasing distance from the current point

            cluster = [] # time of events selected to compute beta value of the current point
            i=0
            while info[i][3]<self.spacemin : # keep events within 50 km radius
                cluster.append(info[i][2])
                i=i+1
           
            if len(cluster)<self.nbmin : # if not enough events

                spmin=self.spacemin+10000 # extend area
                while len(cluster)<self.nbmin and spmin<self.spacemax:#increase radius until either there are enough events or the max radius is reached
                    while info[i][3]<spmin :
                        cluster.append(info[i][2])
                        i+=1
                    spmin+=10000

                if spmin>=self.spacemax :# maximum radius reached before minimum nb reached
                    BETAS.append(np.nan)

                else : # minimum number reached
                    beta=evaluate_beta(self.changetime,sorted(cluster),self.Dt2,self.Dt1)[0]
                    print(beta)
                    if beta>=2 :
                        betainf=evaluate_beta(self.changetime, sorted(cluster),self.Dt2,self.Dt1)[1]
                        BETAS.append(betainf) # significant increase
                        print(betainf)
                    elif -2<beta<2 :
                        BETAS.append(0) # insignificant change

                    elif -2>=beta :
                        betasup=evaluate_beta(self.changetime, sorted(cluster),self.Dt2,self.Dt1)[2]
                        BETAS.append(betasup) # significant decrease
                    
                
            else : # enough events in the minimum radius
                beta=evaluate_beta(self.changetime,sorted(cluster),self.Dt2,self.Dt1)[0]
                print(beta)
                if beta>=2 :
                    betainf=evaluate_beta(self.changetime, sorted(cluster),self.Dt2,self.Dt1)[1]
                    BETAS.append(betainf) # significant increase
                    print(betainf)
                elif -2<beta<2 :
                    BETAS.append(0) # insignificant change

                elif -2>=beta :
                    betasup=evaluate_beta(self.changetime, sorted(cluster),self.Dt2,self.Dt1)[2]
                    BETAS.append(betasup) # significant decrease

        # write results on a file
        file=open('./USED_DATA_FILES/'+str(self.infile)[18:-4]+'beta.txt','w') # removing the first '['character
        file.write(str(self.x)[1:])
        file.write(str(self.y))
        file.write(str(BETAS)[:-1]) # removing the last ']' character
        file.close()
        return

    def drawmap(self) :
        ### Draw map of the change status ###

        file=open('./USED_DATA_FILES/'+str(self.infile)[18:-4]+'beta.txt','r')
        line=file.readline()
        self.x,self.y,self.z=[[float(x) for x in s.split(',')] for s in line.split ('][')]
        file.close
        
        fig=plt.figure(figsize=(6,8))
        if self.lonleft*self.lonright>=0 :
            lon_0=self.lonleft+(self.lonright-self.lonleft)/2
        else :
            lon_0=-360+(360+self.lonright-self.lonleft)/2+self.lonleft
        lat_0 = self.latbottom+ (self.lattop-self.latbottom)/2
        m = Basemap (llcrnrlon=self.lonleft,llcrnrlat=self.latbottom,urcrnrlon=self.lonright,urcrnrlat=self.lattop,projection='lcc',resolution='c',lat_0=lat_0,lon_0=lon_0)
        parallels = np.arange (self.latbottom,self.lattop,4)
        m.drawparallels (parallels,labels=[1,0,0,0],fontsize=15,linewidth=0)
        meridians = np.arange (174,180,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        meridians = np.arange (-178,-172,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        im=plt.imread('TONGA.png')
        m.imshow (np.flipud (im))
    
        yellow = np.array( [1,         1,        0.6 ,       1,       ]) # create colormap
        orange = np.array([1,         0.49803922, 0,         1,        ])
        green = ([0.4 ,       0.65098039, 0.11764706, 1        ])
        newcmp=ListedColormap([green,yellow,orange])

        ngridy=1000 # interpolate data on a regular grid
        ngridx=2000
        self.x,self.y=m(self.x,self.y)
        self.z=np.array(self.z)
        ptl,shit=m(self.lonleft,self.lattop)
        ptr,shit=m(self.lonright,self.lattop)
        shit,ptt=m(self.lonleft,self.lattop)
        shit,ptb=m(self.lonleft,self.latbottom)
        xi=np.linspace(ptl-10000,ptr+10000,ngridx)
        yi=np.linspace(ptb-10000,ptt+10000,ngridy)
        triang = tri.Triangulation(self.x, self.y)
        interpolator = tri.LinearTriInterpolator(triang, self.z)
        Xi, Yi = np.meshgrid(xi, yi)
        zi= interpolator.__call__(Xi,Yi)
        cs=m.contourf(Xi,Yi,zi,levels=[-6,-2,2,6],cmap=newcmp) # draw contour filled
        cbar=m.colorbar(cs,location='right',pad=0.05,ticks=[]) # colorbar
        cbar.ax.set_yticklabels(['Significant decrease', 'No change', 'Significant increase']) #label
        plt.savefig('./RESULTS/'+str(self.infile)[18:-4]+'map.png') # keep map
        plt.show()
        return