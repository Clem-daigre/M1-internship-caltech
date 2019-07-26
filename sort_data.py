
def sortdata(infile,output,latbottom,lattop,lonleft,lonright,starttime,endtime,depthmin,depthmax,magnmin,magnmax) :
    """Read infile, sort data and write on outputfile selected data """

    #read
    file=open(infile,'r')
    line=file.readline()
    depth,lat,lon,magn,time=[[float(x) for x in s.split(',')] for s in line.split ('][')]
    file.close()

    #sort
    area_index=[]
    if lonright*lonleft>=0 :#both positive or both negative                                                                                                                   
        for i in range(len(lat)) :
            if latbottom<=lat[i]<=lattop and lonleft<=lon[i]<=lonright and depthmin<=depth[i]<=depthmax :
                area_index.append(i)
    elif lonright*lonleft<=0 :
        for i in range(len(lat)) :
            if latbottom<=lat[i]<=lattop and (lonleft<=lon[i] or lon[i]<=lonright) :
                area_index.append(i)

    selected_index=[]
    for i in area_index :
        if depthmin<=depth[i]<=depthmax and starttime <=time[i]<=endtime and magnmin<=magn[i]<=magnmax:
                selected_index.append(i)

    #output
    out=open(output, 'w')
    out.write(str([depth[i] for i in selected_index])[1:])#removing the first '['character to easily read the file with a '][' split
    out.write(str([lat[i] for i in selected_index]))
    out.write(str([lon[i] for i in selected_index]))
    out.write(str([magn[i] for i in selected_index]))
    out.write(str([time[i]- starttime for i in selected_index])[:-1])#removing the last ']' character
    return