from extract import *
from relocalize import *
from window_execute import *
from window_steps import * # lpd process with plotting and values printed at each step
from betamap import *
from sort_data import *
from plots import * # plot main results figures for presentations

########################################              EXTRACT DATA           #################################################################
extract=False

if extract==True :
    print('Extracting the data from raw catalogs')
    tonga = extract('RAW_DATA_FILES/ISC_RAWDATA.txt',None,'ISC',2005,1,1,00,00,00,'USED_DATA_FILES/ISC_fulldata.txt')
    tonga.extract_data()
    tonga.output()
    tonga = extract('RAW_DATA_FILES/NEIC_RAWDATA.txt',None,'NEIC',2005,1,1,00,00,00,'USED_DATA_FILES/NEIC_fulldata.txt')
    tonga.extract_data()
    tonga.output()
    tonga = extract('RAW_DATA_FILES/CMTSOLUTION_RAWDATA.txt','RAW_DATA_FILES/PSMECA_RAWDEPTHS.txt','CMT',2005,1,1,00,00,00,'USED_DATA_FILES/SURROUNDING_EQ.txt')
    tonga.extract_data()
    tonga.output()


########################################             RELOCALIZE  DATA           #################################################################
reloc=False

if reloc==True :
    print ('Relocation of all the data')
    r = relocation('./USED_DATA_FILES/ISC_fulldata.txt', './USED_DATA_FILES/ISC_fullrelocalizeddata.txt',2005,1,1,-14,-38,174,-170)
    r.read()
    r.reloc()
    r.out()
    print('All unknowns events assigned to deep')
    r = relocation('./USED_DATA_FILES/ISC_fulldata.txt', './USED_DATA_FILES/ISC_fulldata_unknownstodeep.txt',2005,1,1,-14,-38,174,-170)
    r.read()
    r.deep_to_unknowns()
    r.out()
    print('Introduction figure') # introduction figure with events from 2016/08/18 to 2018/08/18
    sortdata('./USED_DATA_FILES/ISC_fulldata.txt', './USED_DATA_FILES/sorted.txt',-38,-14,174,-170,4247,4977,0,3000,0,10) # infile,outfile,latbottom,lattop,lonleft,lonright,startime,endtime,depthmin,depthmax,magnmin,magnmax
    r = relocation('./USED_DATA_FILES/sorted.txt', './USED_DATA_FILES/ISC_fullrelocalizeddata.txt',2005,1,1,-14,-38,174,-170)
    r.read()
    r.mapintro()
    print('Relocation example for figures') # relocation example for maps from 2005/01/01 to 2007/01/01
    sortdata('./USED_DATA_FILES/ISC_fulldata.txt', './USED_DATA_FILES/sorted.txt',-38,-14,174,-170,0,730,0,3000,0,10) # infile,outfile,latbottom,lattop,lonleft,lonright,startime,endtime,depthmin,depthmax,magnmin,magnmax
    r = relocation('./USED_DATA_FILES/sorted.txt', './USED_DATA_FILES/ISC_fullrelocalizeddata.txt',2005,1,1,-14,-38,174,-170)
    r.read()
    r.reloc()
    r.maps()


########################################             TEST Mw ABOVE 4              #################################################################

##########   ISC DATA   ###############
iscabove4=False

if iscabove4== True :
    print ('ISC deep events in the doublet area Mw>4')
    sortdata('./USED_DATA_FILES/ISC_fullrelocalizeddata.txt', './USED_DATA_FILES/reloc>4isc.txt',-21,-16,178,-179,0,7000,300,3000,4,10) #doublet area above mw 4, deep events plus relocalized as deep
    test=signtest('./USED_DATA_FILES/reloc>4isc.txt',40,4000)# 40 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    test.figure()
    sortdata('./USED_DATA_FILES/ISC_fulldata_unknownstodeep.txt', './USED_DATA_FILES/deepplusunknowns>4isc.txt',-21,-16,178,-179,0,7000,300,3000,4,10) #doublet area above mw 4, deep events plus all unknowns
    test=signtest('./USED_DATA_FILES/deepplusunknowns>4isc.txt',40,4000)# 40 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    test.figure()
    sortdata('./USED_DATA_FILES/ISC_fulldata.txt', './USED_DATA_FILES/deeponly>4isc.txt',-21,-16,178,-179,0,7000,300,3000,4,10) #doublet area above mw 4, known deep events only
    test=signtest('./USED_DATA_FILES/deeponly>4isc.txt',40,4000)# 40 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    test.figure()

########   NEIC DATA  ############
neic=False

if neic==True :
    print('NEIC deep events')
    sortdata('./USED_DATA_FILES/NEIC_fulldata.txt', './USED_DATA_FILES/neicdoublet.txt',-21,-16,178,-179,0,7000,300,3000,4,10) #doublet area, deep events
    test=signtest('./USED_DATA_FILES/neicdoublet.txt',60,4000)# 60 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    neic_cp_values = test.cpvalues
    test.figure()
    sortdata('./USED_DATA_FILES/NEIC_fulldata.txt', './USED_DATA_FILES/neicsidearea.txt',-21,-16,-178,-175,0,7000,300,3000,4,10) #nearby area, deep events
    test=signtest('./USED_DATA_FILES/neicsidearea.txt',20,4000)# 20 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    test.figure()

    print('NEIC intermediate and shallow events')
    sortdata('./USED_DATA_FILES/NEIC_fulldata.txt', './USED_DATA_FILES/neicinterm.txt',-21,-16,-178,-175,0,7000,70,300,4,10) #int and shallow area, intermediate events
    test=signtest('./USED_DATA_FILES/neicinterm.txt',40,4000)# 40 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    test.figure()
    sortdata('./USED_DATA_FILES/NEIC_fulldata.txt', './USED_DATA_FILES/neicshallow.txt',-21,-16,-178,-175,0,7000,0,70,4,10) #int and shallow area, shallow events
    test=signtest('./USED_DATA_FILES/neicshallow.txt',60,4000)# 60 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    test.figure()
    #neic_cp_values=[0,1530.0,3150.0,4920.583124652971]#with 0 value and last time value
    print('Beta maps of doublet area deep changes')
    for i in range(len(neic_cp_values)-2) :
        sortdata('./USED_DATA_FILES/NEIC_fulldata.txt','./USED_DATA_FILES/neicchange'+str(i+1)+'.txt',-38,-14,174,-170,neic_cp_values[i],neic_cp_values[i+2],300,3000,0,10 )# data set from previous change to next change
        mapping=map('./USED_DATA_FILES/neicchange'+str(i+1)+'.txt',neic_cp_values[i+1]-neic_cp_values[i],100,50000,50,50,200000,-38,-14,174,-170)#beta maps of change point i+1 #infile,changetime,minnb,spacemin,lonrange,latrange,spacemax,latbottom,lattop,lonleft,lonright
        #minnb of 100 events for the 2 neic changes
        mapping.priliminary()
        mapping.compute()
        mapping.drawmap()

   
########################################             TEST Mw BELOW 4              #################################################################
iscbelow4=False

if iscbelow4==True:
    print ('ISC deep events in the doublet area Mw < 4')
    sortdata('./USED_DATA_FILES/ISC_fullrelocalizeddata.txt', './USED_DATA_FILES/reloc<4isc.txt',-21,-16,178,-179,0,7000,300,3000,0,4) #doublet area below mw 4, deep events plus relocalized as deep
    test=signtest('./USED_DATA_FILES/reloc<4isc.txt',20,4000)# 20 days interval chosen by testing
    test.read()
    test.proba()
    isc_cp_values = test.sign()
    test.figure()
    sortdata('./USED_DATA_FILES/ISC_fulldata_unknownstodeep.txt', './USED_DATA_FILES/deepplusunknowns<4isc.txt',-21,-16,178,-179,0,7000,300,3000,0,4) #doublet area below mw 4, deep events plus all unknowns
    test=signtest('./USED_DATA_FILES/deepplusunknowns<4isc.txt',20,4000)# 20 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    test.figure()
    sortdata('./USED_DATA_FILES/ISC_fulldata.txt', './USED_DATA_FILES/deeponly<4isc.txt',-21,-16,178,-179,0,7000,300,3000,0,4) #doublet area below mw 4, known deep events only
    test=signtest('./USED_DATA_FILES/deeponly<4isc.txt',20,4000)# 20 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    test.figure()

    print('ISC intermediate and shallow events Mw<4')
    sortdata('./USED_DATA_FILES/ISC_fullrelocalizeddata.txt', './USED_DATA_FILES/interm<4isc.txt',-21,-16,-178,-175,0,7000,70,300,0,4) #int and shallow area, intermediate events, mw<4
    test=signtest('./USED_DATA_FILES/interm<4isc.txt',40,4000)# 40 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    test.figure()
    sortdata('./USED_DATA_FILES/ISC_fullrelocalizeddata.txt', './USED_DATA_FILES/shallow<4isc.txt',-21,-16,-178,-175,0,7000,0.01,70,0,4)  #int and shallow area, shallow events, mw<4
    test=signtest('./USED_DATA_FILES/shallow<4isc.txt',40,4000)# 40 days interval chosen by testing
    test.read()
    test.proba()
    test.sign()
    test.figure()   
    print('Beta maps of doublet area deep changes')
    #isc_cp_values=[0,1310,2630,3570,4290,4730,6000]#with 0 value and last time value
    for i in range(len(isc_cp_values)-2) :
        sortdata('./USED_DATA_FILES/ISC_fullrelocalizeddata.txt','./USED_DATA_FILES/isc<4change'+str(i+1)+'.txt',-38,-14,174,-170,isc_cp_values[i],isc_cp_values[i+2],300,3000,0,4)
        mapping=map('./USED_DATA_FILES/isc<4change'+str(i+1)+'.txt',isc_cp_values[i+1]-isc_cp_values[i],150,50000,50,50,200000,-38,-14,174,-170)#minnb default 150, can be changed
        mapping.priliminary()
        mapping.compute()
        mapping.drawmap()


########################################            CHECK SURROUNDING EVENTS (PLOTS FUNCTION)         #################################################################
isc_cp_values=[1310,2630,3570,4290,4730]#without 0 and last time value

trig=False
if trig==True :
    plot_trig(isc_cp_values,'./USED_DATA_FILES/reloc<4isc.txt','./USED_DATA_FILES/SURROUNDING_EQ.txt')

########################################            MAIN RESULTS FIGURES (PLOTS FUNCTION)           #################################################################
res=False
if res== True :
    plot_cpres(isc_cp_values,'./USED_DATA_FILES/reloc<4isc.txt','./USED_DATA_FILES/deepplusunknowns<4isc.txt','./USED_DATA_FILES/deeponly<4isc.txt')
    betamaps () 