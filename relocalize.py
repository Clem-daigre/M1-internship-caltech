import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,cm


def select_area_depth (latlist,lonlist,depthlist,latt,latb,lonl,lonr,depthmin,depthmax) :
    ### Select events from latlist, lonlist and depthlist that are within a given depth interval and rectangle area ###

    index_selected=[]                                                                                                                 
    for i in range (len(latlist)) :
        if latb<=latlist[i]<=latt and lonl<=lonlist[i]<=lonr and depthmin<=depthlist[i]<=depthmax :
            index_selected.append (i)
    return index_selected

class relocation () :
    """ Read a file, assign unknown depth events to a category and return the modified data on an output file. Represent maps as an example."""
    
    def __init__ (self,infile,outfile,startyr,startmt,startday,lattop,latbottom,lonleft,lonright) :
        ### initialisation ###
        
        self.infile=infile
        self.outfile=outfile
        self.startyr=startyr
        self.startmt=startmt
        self.startday=startday
        self.lattop=lattop # edges of the big area studied
        self.latbottom=latbottom
        self.lonleft=lonleft
        self.lonright=lonright
        return


    def read (self) :
        ### Read data, seperate events with known depths in the 3 categories ###
        
        file=open (self.infile,'r')# read infile
        line=file.readline ()
        self.depth,self.lat,self.lon,self.magn,self.time=[[float(x) for x in s.split (',')] for s in line.split ('][')]
        file.close ()
        if self.lonleft > 0 : # transform all positive longitude values to negative ones
            self.lonleft = -360+self.lonleft # 174 E becomes -186 W
        for i in range (len(self.lon)) :
            if self.lon[i] > 0 :
                self.lon[i] = -360+self.lon[i]

        self.latint, self.lonint = [], [] # known int events
        self.latsh, self.lonsh = [], [] # known shallow events
        self.latdeep, self.londeep = [], [] # known deep events
        self.latzerobefore,self.lonzerobefore = [], [] # unknown depth events before relocalisation
        for i in range(len(self.depth)) :
            if self.depth[i]>=300 :
                self.latdeep.append(self.lat[i])
                self.londeep.append(self.lon[i])
            elif 0.01<=self.depth[i]<=70 :
                self.latsh.append(self.lat[i])
                self.lonsh.append(self.lon[i])
            elif 300>=self.depth[i]>=70 :
                self.latint.append(self.lat[i])
                self.lonint.append(self.lon[i])
            elif self.depth[i]==0.0 :
                self.latzerobefore.append(self.lat[i])
                self.lonzerobefore.append(self.lon[i])
        print('Pourcentage of unknowns events before', len(self.latzerobefore)/len(self.depth))
        return

    def mapintro (self) :
        ###  Plot the map of Tonga settings for the introduction ###

        fig = plt.figure (figsize= (6,8))
        lon_0 = self.lonleft+ (self.lonright-self.lonleft)/2
        lat_0 = self.latbottom+ (self.lattop-self.latbottom)/2
        m = Basemap (llcrnrlon=self.lonleft,llcrnrlat=self.latbottom,urcrnrlon=self.lonright,urcrnrlat=self.lattop,projection='lcc',resolution='c',lat_0=lat_0,lon_0=lon_0)
        parallels = np.arange (self.latbottom,self.lattop,4)
        m.drawparallels (parallels,labels=[1,0,0,0],fontsize=15,linewidth=0)
        meridians = np.arange (174,180,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        meridians = np.arange (-178,-172,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        im=plt.imread ('TONGA.png') # arcgis image
        m.imshow (np.flipud (im))
        data=np.loadtxt('pacific_boundary.txt',dtype='str',delimiter=' ',unpack=True)#plate boundary
        latp=[]
        lonp=[]
        for l in data[1][:] :
            latp.append(float(l))
        for l in data[2][:] :
            lonp.append(float(l))
        lonp,latp=m(lonp,latp)
        m.plot(lonp,latp,color='k',linewidth=2)
        
        lonsh,latsh = m (self.lonsh,self.latsh) # known as shallow
        m.scatter(lonsh,latsh,s=0.3,label='shallow events',color='gold')
        lonint,latint = m (self.lonint,self.latint) # known as intermediate
        m.scatter(lonint,latint,s=0.3,label='intermediate-deep events',color='darkorange')
        londeep,latdeep = m (self.londeep,self.latdeep) # known as deep
        m.scatter(londeep,latdeep,s=0.3,label='deep-focus events',color='red')
        
        lon82,lat82 = m(-178.153,-18.113) # 8.2 mw
        lon79,lat79 = m(-180.65,-18.474) # 7.9 mw
        m.scatter(lon82,lat82,s=200,marker='*',color='lime')
        m.scatter(lon79,lat79,s=200,marker='*',color='darkturquoise')
        plt.legend(loc='lower right',fontsize=13,markerscale=7)
        plt.show ()
        return

    def deep_to_unknowns (self) : 
        ### Assign all unknown depths events to deep focus category ###

        self.depthreloc = self.depth # depths will be modified on this list
        unknowns_ev = select_area_depth (self.lat,self.lon,self.depth,self.lattop,self.latbottom,self.lonleft,self.lonright,0,0) # indexes of the unknown events in the whole area
        for i in unknowns_ev :
            self.depthreloc[i] = 500.0
        return

    def reloc (self):
        ### Determine unknown events depth category from their horizontal location###
        
        self.lat0int, self.lon0int = [], [] # relocalized as int events  ## keep track of events for further mapping
        self.lat0sh, self.lon0sh = [], [] # relocalized as shallow events
        self.lat0deep, self.lon0deep = [], [] # relocalized as deep events
        self.depthreloc = self.depth # depths will be modified on this list

        latt=self.lattop # edges of the square to be analized, latb=latbottom, latt=lattop,lonr=lonright,lonl=lonleft
        latb=latt-1 # one degree latitude rectangles

        while latb >= self.latbottom : # iterating the latitudes from top to bottom
            lonl = self.lonleft
            lonr = lonl+0.5 # 0.5 degree longitude rectangles
            while lonr <= self.lonright : # iterating the longitudes from left to right
                lengthdeep=len(select_area_depth(self.lat,self.lon,self.depth,latt,latb,lonl,lonr,300,3000)) # number of deep events in the current rectangle
                lengthint=len(select_area_depth(self.lat,self.lon,self.depth,latt,latb,lonl,lonr,70,300)) # number of intermediate events in the current rectangle
                lengthsh=len(select_area_depth(self.lat,self.lon,self.depth,latt,latb,lonl,lonr,0.1,70)) # number of shallow events in the current rectangle#begins at 0.01 km because unknown depths put as zero
                unknowns_indexes = select_area_depth (self.lat,self.lon,self.depth,latt,latb,lonl,lonr,0,0) # indexes of the unknown events in the current rectangle
                
                if np.max([lengthdeep,lengthint,lengthsh])>=1 and np.max([lengthdeep,lengthint,lengthsh])==lengthdeep:
                    for i in unknowns_indexes :
                        self.depthreloc[i]=500.0
                        self.lat0deep.append(self.lat[i])
                        self.lon0deep.append(self.lon[i])
                elif np.max([lengthdeep,lengthint,lengthsh])>=1 and np.max([lengthdeep,lengthint,lengthsh])==lengthint:
                    for i in unknowns_indexes :
                        self.depthreloc[i]=185.0
                        self.lat0int.append(self.lat[i])
                        self.lon0int.append(self.lon[i])
                elif np.max([lengthdeep,lengthint,lengthsh])>=1 and np.max([lengthdeep,lengthint,lengthsh])==lengthsh :
                    for i in unknowns_indexes :
                        self.depthreloc[i]=35.0
                        self.lat0sh.append(self.lat[i])
                        self.lon0sh.append(self.lon[i])
                
                lonl+=0.5 #iterating
                lonr+=0.5
            latb=latb-1
            latt=latt-1 

        unreloc = select_area_depth (self.lat,self.lon,self.depth,self.lattop,self.latbottom,self.lonleft,self.lonright,0,0)#unknowns events in the whole area after relocation
        self.lonzeroafter,self.latzeroafter = [], [] #unrelocalized events
        for i in unreloc :
            self.latzeroafter.append (self.lat[i])
            self.lonzeroafter.append (self.lon[i])
        print ('Pourcentage of unknowns events after', len (unreloc)/len (self.depthreloc))
        return
    

    def out (self) : 
        ### Write data with depths modified on outfile ### 
        
        out=open (self.outfile, 'w')
        out.write (str(self.depthreloc)[1:])
        out.write (str(self.lat))
        out.write (str(self.lon))
        out.write (str(self.magn))
        out.write (str(self.time)[:-1])
        out.close ()
        return


    def maps(self) :
        ### Draw maps of events to check the relocation process was accurate ###

        # defining base of maps and converting events into maps coordinates
        lon_0 = self.lonleft+ (self.lonright-self.lonleft)/2
        lat_0 = self.latbottom+ (self.lattop-self.latbottom)/2
        m = Basemap (llcrnrlon=self.lonleft,llcrnrlat=self.latbottom,urcrnrlon=self.lonright,urcrnrlat=self.lattop,projection='lcc',resolution='c',lat_0=lat_0,lon_0=lon_0)
        londeep,latdeep = m (self.londeep,self.latdeep) # known as deep
        lonsh,latsh = m (self.lonsh,self.latsh) # known as shallow
        lonint,latint = m (self.lonint,self.latint) # known as intermediate
        lon0deep,lat0deep = m (self.lon0deep,self.lat0deep) # relocalized as deep
        lon0sh,lat0sh = m (self.lon0sh,self.lat0sh) # relocalized as shallow
        lon0int,lat0int = m (self.lon0int,self.lat0int) # relocalized as intermediate
        lonzerobefore,latzerobefore = m (self.lonzerobefore,self.latzerobefore) # unknown events before relocalization
        lonzeroafter,latzeroafter = m (self.lonzeroafter,self.latzeroafter) # unrelocalized events
        im=plt.imread ('TONGA.png') # arcgis image

        # map 1 before relocation
        fig = plt.figure (figsize= (6,8))
        m = Basemap (llcrnrlon=self.lonleft,llcrnrlat=self.latbottom,urcrnrlon=self.lonright,urcrnrlat=self.lattop,projection='lcc',resolution='c',lat_0=lat_0,lon_0=lon_0)
        parallels = np.arange (self.latbottom,self.lattop,4)
        m.drawparallels (parallels,labels=[1,0,0,0],fontsize=15,linewidth=0)
        meridians = np.arange (174,180,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        meridians = np.arange (-178,-172,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        m.imshow (np.flipud (im))
        m.scatter (londeep,latdeep,s=3,label='known as deep',color='red')
        m.scatter (lonsh,latsh,s=3,label='known as shallow',color='darkorange')
        m.scatter (lonint,latint,s=3,label='known as intermediate ',color='lime')
        m.scatter (lonzerobefore,latzerobefore,s=3,label='unknown depth',color='magenta')
        plt.legend (loc='lower right',fontsize=12,markerscale=4)
        plt.show ()

        # map2 deep after relocation
        fig = plt.figure(figsize=(6,8))
        m = Basemap (llcrnrlon=self.lonleft,llcrnrlat=self.latbottom,urcrnrlon=self.lonright,urcrnrlat=self.lattop,projection='lcc',resolution='c',lat_0=lat_0,lon_0=lon_0)
        parallels = np.arange (self.latbottom,self.lattop,4)
        m.drawparallels (parallels,labels=[1,0,0,0],fontsize=15,linewidth=0)
        meridians = np.arange (174,180,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        meridians = np.arange (-178,-172,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        m.imshow (np.flipud (im))
        m.scatter (londeep,latdeep,s=3,label='known deep events',color='darkorange')
        m.scatter (lon0deep,lat0deep,s=3,label='relocalized as deep',color='lime')
        plt.legend (loc='lower right',fontsize=12,markerscale=4)
        plt.show  ()

        # map3 intermediate after relocation
        fig = plt.figure (figsize= (6,8))
        m = Basemap (llcrnrlon=self.lonleft,llcrnrlat=self.latbottom,urcrnrlon=self.lonright,urcrnrlat=self.lattop,projection='lcc',resolution='c',lat_0=lat_0,lon_0=lon_0)
        parallels = np.arange (self.latbottom,self.lattop,4)
        m.drawparallels (parallels,labels=[1,0,0,0],fontsize=15,linewidth=0)
        meridians = np.arange (174,180,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        meridians = np.arange (-178,-172,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        m.imshow (np.flipud (im))
        m.scatter (lonsh,latsh,s=3,label='known shallow events',color='darkorange')
        m.scatter (lon0sh,lat0sh,s=3,label='relocalized as shallow',color='lime')
        plt.legend (loc='lower right',fontsize=12,markerscale=4)
        plt.show ()

        # map4 shallow after relocation
        fig = plt.figure (figsize= (6,8))
        m = Basemap (llcrnrlon=self.lonleft,llcrnrlat=self.latbottom,urcrnrlon=self.lonright,urcrnrlat=self.lattop,projection='lcc',resolution='c',lat_0=lat_0,lon_0=lon_0)
        parallels = np.arange (self.latbottom,self.lattop,4)
        m.drawparallels (parallels,labels=[1,0,0,0],fontsize=15,linewidth=0)
        meridians = np.arange (174,180,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        meridians = np.arange (-178,-172,4)
        m.drawmeridians (meridians,labels=[0,0,0,1],fontsize=15,linewidth=0)
        m.imshow (np.flipud (im))
        m.scatter (lonint,latint,s=3,label='known intermediate events',color='darkorange')
        m.scatter (lon0int,lat0int,s=3,label='relocalized as intermediate',color='lime')
        m.scatter (lonzeroafter,latzeroafter,s=3,label='leftover unknown depths events',color='magenta')
        plt.legend (loc='lower right',fontsize=12,markerscale=4)
        plt.show ()
        return