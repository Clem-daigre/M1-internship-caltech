import numpy as np 

def julian_date(y,m,d,hour,minute,sec) : 
    ### computes the julian date from a calendar date ###
    
	a = (14-m)/12
	ny = y+4800-a
	nm = m+12*a-3
	jd = d+(153*nm+2)/5 + 365*ny + (ny/4) -ny/100 +ny/400- 32045
	fh=hour + minute/60 + sec/3600
	jd = jd-0.5+fh/24
	return jd

class extract() : 
    """ Input is a raw catalogue text file, with one line = one event informations, separated by a str delimiter. Catalog style must be given (ISC,NEIC or CMT). Events must be ranked with increasing time from start date.
    The class read the input file, compute the time in days of events since the start time (2005-01-01) thanks to the julian date, and write extracted data on the output file.
    Output is lists of depths of events,latitudes of events,...such that elements of index i in all lists are corresponding to the same event."""

    def __init__(self, infile, depthinfile,catalogue,startyr,startmt,startday,starth,startmin,startsec,outfile):
        ### initialisation ###
    
        self.catalogue=catalogue
        self.infile=infile
        self.depthinfile=depthinfile
        self.startyr=startyr
        self.startmt=startmt
        self.startday=startday
        self.starth=starth
        self.startmin=startmin
        self.startsec=startsec
        self.outfile=outfile
        return
    
    def extract_data(self) :
        ### extract data from lines in catalogs to lists of events info ###

        if self.catalogue=='NEIC' : # NEIC format
            latitude,longitude,depth,magn=np.loadtxt(self.infile,dtype='float',delimiter=',', usecols=[1,2,3,4],unpack=True) # extract data
            self.lat,self.lon,self.depth,self.magn=[],[],[],[]#we want lists of values separated by comas
            for l in latitude :
                self.lat.append(l)
            for l in longitude :
                self.lon.append(l)
            for d in depth :
                self.depth.append(d)
            for m in magn :
                self.magn.append(m)
            date = np.loadtxt(self.infile, dtype='str', delimiter=',', usecols=[0],unpack=True)
            yr,month,day=[],[],[] # extract time
            h,minute,sec=[],[],[]
            for d in date:
                yr.append(int(d[0:4]))
                month.append(int(d[5:7]))
                day.append(int(d[8:10]))
                h.append(int(d[11:13]))
                minute.append(int(d[14:16]))
                sec.append(float(d[17:22]))
            start_jd = julian_date(self.startyr,self.startmt,self.startday,self.starth,self.startmin,self.startsec)#julian date of the start time
            self.time = [julian_date(yr[i],month[i],day[i],h[i],minute[i],sec[i])- start_jd for i in range(len(yr))]#julian date of the event minus julian date of the start time
        

        elif self.catalogue=='ISC' : # ISC format
            latitude,longitude,depth,magn=np.loadtxt(self.infile,dtype='float',delimiter=',', usecols=[4,5,6,10],unpack=True)#extract data
            self.lat,self.lon,self.depth,self.magn=[],[],[],[]#we want lists of values separated by comas
            for l in latitude :
                self.lat.append(l)
            for l in longitude :
                self.lon.append(l)
            for d in depth :
                self.depth.append(d)
            for m in magn :
                self.magn.append(m)
            date, hour = np.loadtxt(self.infile, dtype='str', delimiter=',', usecols=[2,3],unpack=True)#extract time
            yr,month,day=[],[],[]
            h,minute,sec=[],[],[]
            for d in date:
                yr.append(int(d[0:4]))
                month.append(int(d[5:7]))
                day.append(int(d[8:10]))
            for ho in hour :
                h.append(int(ho[0:2]))
                minute.append(int(ho[3:5]))
                sec.append(float(ho[6:]))
            start_jd = julian_date(self.startyr,self.startmt,self.startday,self.starth,self.startmin,self.startsec)#julian date of the start time
            self.time = [julian_date(yr[i],month[i],day[i],h[i],minute[i],sec[i])- start_jd for i in range(len(yr))]#julian date of the event minus julian date of the start time
        
        
        elif self.catalogue=='CMT' : # CMT format
            time,coor,magn=np.loadtxt(self.infile,dtype='str', delimiter=' | ',unpack=True)
            self.magn=[float(m[0:4]) for m in magn]# extract magn
            self.lat=[]#extract coord
            self.lon=[]
            for c in coor :
                self.lat.append(float(c[0:7]))
                self.lon.append(float(c[9:17]))
            yr,month,day,h,minute,sec=[],[],[],[],[],[] # extract tine
            for t in time :
                yr.append(int(t[0:4]))
                month.append(int(t[5:7]))
                day.append(int(t[8:10]))
                h.append(int(t[11:13]))
                minute.append(int(t[14:16]))
                sec.append(float(t[17:22]))
            start_jd = julian_date(self.startyr,self.startmt,self.startday,self.starth,self.startmin,self.startsec)#julian date of the start time
            self.time = [julian_date(yr[i],month[i],day[i],h[i],minute[i],sec[i])- start_jd for i in range(len(yr))]#julian date of the event minus julian date of the start time        
            depth=np.loadtxt(self.depthinfile,dtype='str',delimiter=' ', usecols=[2],unpack=True)#extract depth
            self.depth=[int(d) for d in depth]


        return


    def output(self) :
        ### write extracted data on the outputfile ###

        out=open(self.outfile, 'w')
        out.write(str(self.depth)[1:])#removing the first '['character to easily read the file with a '][' split
        out.write(str(self.lat))
        out.write(str(self.lon))
        out.write(str(self.magn))
        out.write(str(self.time)[:-1])#removing the last ']' character
        return
 
