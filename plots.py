import matplotlib.pyplot as plt
import numpy as np
from obspy.geodetics import base
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
import matplotlib.tri as tri
from scipy.interpolate import griddata

""" Various useful plotting functions """

def plot_cpres(cpvalues,reloc,dpunknowns,donly) :
    ### Plot the 3 cumulative curves deep+unknowns, deep only and deep+reloc with the cp values found. ###

    file=open(reloc,'r')
    line=file.readline()
    timereloc=[float(x) for x in line.split ('][')[4].split(',')]
    file.close()
    file=open(dpunknowns,'r')
    line=file.readline()
    timedpunknowns=[float(x) for x in line.split ('][')[4].split(',')]
    file.close()
    file=open(donly,'r')
    line=file.readline()
    timedonly=[float(x) for x in line.split ('][')[4].split(',')]
    file.close()
    plt.figure()
    plt.plot(timedonly, [i+1 for i in range(len(timedonly))], 'navy', linewidth=1.25,label='Known deep events')
    plt.plot(timereloc, [i+1 for i in range(len(timereloc))], 'deepskyblue', linewidth=1.25,label='Known deep+relocalized as deep events')
    plt.plot(timedpunknowns, [i+1 for i in range(len(timedpunknowns))], 'crimson', linewidth=1.25,label='Known deep+all unknown depth events')
    plt.xlabel('Date',size=15)
    locs=[0,1000,2000,3000,4000,4800]
    labels=['2005-1','2007-9','2010-6','2013-3','2015-12','2018-2']
    plt.xticks(locs,labels,fontsize=13)
    plt.yticks(fontsize=13)
    plt.ylabel('Cumulative number of events',size=15)
    for i in range(0,2) :
        plt.plot([cpvalues[i],cpvalues[i]],[0,len(timedpunknowns)],color='k',linewidth=1,linestyle='dashed')
    for i in range(2,len(cpvalues)) :
        plt.plot([cpvalues[i],cpvalues[i]],[0,len(timedpunknowns)],color='k',linewidth=1)
    plt.legend(bbox_to_anchor=(0.02, 1.12), loc=2, borderaxespad=0.,fontsize=11)
    plt.savefig('./RESULTS/sumupchanges.png')
    plt.show()
    return

def plot_trig(cpvalues,reloc,surround) :
    ### Plot the 2 figures to test the triggering possibility of the changes.###

    file=open(reloc,'r')
    line=file.readline()
    timereloc=[float(x) for x in line.split ('][')[4].split(',')]
    file.close()
    file=open(surround,'r')
    line=file.readline()
    depthsurround,latsurround,lonsurround,magnsurround,timesurround=[[float(x) for x in s.split(',')] for s in line.split ('][')]
    file.close()
    plt.figure()
    plt.plot(timereloc, [i+1 for i in range(len(timereloc))], 'navy', linewidth=2)
    plt.xlabel('Date',size=15)
    locs=[0,1000,2000,3000,4000,4800]
    labels=['2005-1','2007-9','2010-6','2013-3','2015-12','2018-2']
    plt.xticks(locs,labels,fontsize=13)
    plt.yticks(fontsize=13)
    plt.ylabel('Cumulative number of events',size=15)
    i=0
    plt.plot([cpvalues[i],cpvalues[i]],[0,len(timereloc)],color='r',linewidth=1.5,label='rate changes')
    plt.plot([timesurround[i],timesurround[i]],[0,len(timereloc)],color='k',linewidth=1,label='surrounding events')
    for i in range(1,len(cpvalues)) :
        plt.plot([cpvalues[i],cpvalues[i]],[0,len(timereloc)],color='r',linewidth=1.5)
    for i in range(1,len(timesurround)) :
        plt.plot([timesurround[i],timesurround[i]],[0,len(timereloc)],color='k',linewidth=1)
    plt.legend(bbox_to_anchor=(0.02, 1.12), loc=2, borderaxespad=0.,fontsize=11)
    plt.savefig('./RESULTS/testtrig.png')
    plt.show()

    dist=[base.gps2dist_azimuth(-17.5,180,latsurround[i],lonsurround[i])[0]/1000 for i in range(len(latsurround))]
    plt.figure()
    plt.scatter(dist,magnsurround,color='k')
    for i in (7,9) :
        plt.scatter(dist[i],magnsurround[i],color='r')
    plt.ylabel('Magnitude Mw',fontsize=15)
    plt.yticks(fontsize=13)
    plt.xticks(fontsize=13)
    plt.xlim(xmin=100,xmax=1000)
    plt.ylim(ymin=6.4,ymax=7.4)
    plt.xlabel('Distance to the center of the increase area (km)',fontsize=15)
    plt.grid()
    plt.savefig('./RESULTS/testtrig2.png')
    plt.show()
    return

def betamaps () :
    ### Plot the contour of the area where each change is significant on the same map. ###

    #first map
    fig=plt.figure(figsize=(6,8))
    im=plt.imread('TONGA.png')
    m = Basemap(llcrnrlon=174,llcrnrlat=-38,urcrnrlon=-170,urcrnrlat=-14,projection='lcc',resolution='c',lat_0=-26,lon_0=-178)
    parallels=np.arange(-38,-14,4)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
    meridians = np.arange(174,180,4)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)
    meridians = np.arange(-178,-172,4)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)
    m.imshow(np.flipud(im))
    ngridy=1000
    ngridx=2000

    file=open('./USED_DATA_FILES/isc<4change5beta.txt','r')
    line=file.readline()
    x,y,z=[[float(x) for x in s.split(',')] for s in line.split ('][')]
    file.close
    x,y=m(x,y)
    for i in range(len(z)):
        if str(z[i])=='nan' :
            z[i]=0     
    z=np.array(z)
    ptl,shit=m(174,-14)
    ptr,shit=m(-170,-14)
    shit,ptt=m(174,-14)
    shit,ptb=m(174,-38)
    xi=np.linspace(ptl-10000,ptr+10000,ngridx)
    yi=np.linspace(ptb-10000,ptt+10000,ngridy)
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi= interpolator.__call__(Xi,Yi)
    cs=m.contour(Xi,Yi,zi,colors='yellow',levels=[-2],linestyles='solid',linewidths=[1])
    plt.plot(0,0,color='yellow',label='2017/12/14')

    file=open('./USED_DATA_FILES/isc<4change4beta.txt','r')
    line=file.readline()
    x,y,z=[[float(x) for x in s.split(',')] for s in line.split ('][')]
    file.close
    x,y=m(x,y)
    for i in range(len(z)):
        if str(z[i])=='nan' :
            z[i]=0     
    z=np.array(z)
    ptl,shit=m(174,-14)
    ptr,shit=m(-170,-14)
    shit,ptt=m(174,-14)
    shit,ptb=m(174,-38)
    xi=np.linspace(ptl-10000,ptr+10000,ngridx)
    yi=np.linspace(ptb-10000,ptt+10000,ngridy)
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi= interpolator.__call__(Xi,Yi)
    cs=m.contour(Xi,Yi,zi,colors='limegreen',levels=[2],linestyles='solid',linewidths=[1])
    plt.plot(0,0,color='limegreen',label='2016/09/30')


    file=open('./USED_DATA_FILES/isc<4change3beta.txt','r')
    line=file.readline()
    x,y,z=[[float(x) for x in s.split(',')] for s in line.split ('][')]
    file.close
    x,y=m(x,y)
    for i in range(len(z)):
        if str(z[i])=='nan' :
            z[i]=0     
    z=np.array(z)
    ptl,shit=m(174,-14)
    ptr,shit=m(-170,-14)
    shit,ptt=m(174,-14)
    shit,ptb=m(174,-38)
    xi=np.linspace(ptl-10000,ptr+10000,ngridx)
    yi=np.linspace(ptb-10000,ptt+10000,ngridy)
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi= interpolator.__call__(Xi,Yi)
    cs=m.contour(Xi,Yi,zi,colors='magenta',levels=[2],linestyles='solid',linewidths=[1])
    plt.plot(0,0,color='magenta',label='2014/10/11')
    plt.legend(loc='lower right',fontsize=15)
    plt.savefig('./RESULTS/sumupbetamaps2.png')
    plt.show()

    #map2
    fig=plt.figure(figsize=(6,8))
    im=plt.imread('TONGA.png')
    m = Basemap(llcrnrlon=174,llcrnrlat=-38,urcrnrlon=-170,urcrnrlat=-14,projection='lcc',resolution='c',lat_0=-26,lon_0=-178)
    parallels=np.arange(-38,-14,4)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
    meridians = np.arange(174,180,4)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)
    meridians = np.arange(-178,-172,4)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)
    m.imshow(np.flipud(im))
    ngridy=1000
    ngridx=2000
       
    file=open('./USED_DATA_FILES/isc<4change2beta.txt','r')
    line=file.readline()
    x,y,z=[[float(x) for x in s.split(',')] for s in line.split ('][')]
    file.close
    x,y=m(x,y)
    for i in range(len(z)):
        if str(z[i])=='nan' :
            z[i]=0     
    z=np.array(z)
    ptl,shit=m(174,-14)
    ptr,shit=m(-170,-14)
    shit,ptt=m(174,-14)
    shit,ptb=m(174,-38)
    xi=np.linspace(ptl-10000,ptr+10000,ngridx)
    yi=np.linspace(ptb-10000,ptt+10000,ngridy)
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi= interpolator.__call__(Xi,Yi)
    cs=m.contour(Xi,Yi,zi,colors='cyan',levels=[-2],linestyles='solid',linewidths=[1])
    plt.plot(0,0,color='cyan',label='2012/03/15')
        
    file=open('./USED_DATA_FILES/isc<4change1beta.txt','r')
    line=file.readline()
    x,y,z=[[float(x) for x in s.split(',')] for s in line.split ('][')]
    file.close
    x,y=m(x,y)
    for i in range(len(z)):
        if str(z[i])=='nan' :
            z[i]=0     
    z=np.array(z)
    ptl,shit=m(174,-14)
    ptr,shit=m(-170,-14)
    shit,ptt=m(174,-14)
    shit,ptb=m(174,-38)
    xi=np.linspace(ptl-10000,ptr+10000,ngridx)
    yi=np.linspace(ptb-10000,ptt+10000,ngridy)
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi= interpolator.__call__(Xi,Yi)
    cs=m.contour(Xi,Yi,zi,colors='red',levels=[2],linestyles='solid',linewidths=[1])
    plt.plot(0,0,color='red',label='2008/08/03')
    plt.legend(loc='lower right',fontsize=15)
    plt.savefig('./RESULTS/sumupbetamaps1.png')
    plt.show()
    return