from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import os
import pickle
from os import listdir
from os.path import isfile, join
from shapely.ops import cascaded_union
import numpy as np
from datetime import datetime,timedelta
from pytz import timezone
#import folium
from shapely.ops import cascaded_union


#numero settori 4302
#numero poligoni 93697
#settore nullo 4303

#con soglia a quotecut
#numero settori 4302
#numero poligoni 74250
#settore nullo 4303




NPoints = 500000
QUOTECUT = 250
QUOTECUT_sector = 280
current_path = '.'
start_datetime = "2017-09-01 00:00:00"
end_datetime = "2017-09-01 23:59:59"
confpath = "./config/"
delay_map = ""#./input/tw.txt"
file_conteggi_per_capacity = ""#"./input/conteggi_settori_nuovo_13_350.txt"
file_capacity_tw = "capacita_secondi.out"

filename_poligoni_settori = "settori_trieste_apertura_chiusura_secondi.out"
so6_folder = "./"
#filename_so6_input = "multidelay_"+os.getcwd().split("/")[-1].replace("v","")+".so6"#TriesteM1.so6"
#filename_so6_input = "multidelay_"+os.getcwd().split("/")[-1].replace("v","")+".so6"
filename_so6_input = "20170901_m1.so6"
print(filename_so6_input)
century = '20' # for 2017 etc.. use '19' for 1997 and other years
bound_file = current_path+"/config/boundary/"
temp_nvp = current_path+"/config/sectors_temp_nvp.dat"
shock_tmp = current_path+"/config/temp_nvp.dat"
capacity_file = current_path+"/config/sector_capacities.dat"
#delay_file = current_path+"/config/delay.dat"


lat_max = 82.0
lat_min = 19.0
lon_max = 46.625
lon_min = -30.0



if not os.path.exists(confpath):
    os.makedirs(confpath)
if not os.path.exists(confpath+"boundary/"):
    os.makedirs(confpath+"boundary/")

ritardi = dict()
if delay_map!="":
    with open(delay_map) as fin:
        for line in fin:
            ll = line.replace("\n","").replace("\t","").split(" ")
            fid = int(ll[0])
            delay = int(ll[1])
            if not fid in ritardi:
                ritardi[fid]=delay
npoligons = 0
capacity = dict()
if file_conteggi_per_capacity!='':
    with open(file_conteggi_per_capacity) as fin:
        for line in fin:
            ll = line.replace("\n","").split(" ")
            sector = ll[2]
            cap = int(ll[3])
            if not sector in capacity:
                capacity[sector] = cap
            else:
                if cap> capacity[sector]:
                    capacity[sector] = cap
else:
    if file_capacity_tw!='':
        with open(file_capacity_tw) as fin:
            for line in fin:
                if line[0]=='#':
                    continue
                ll = line.replace("\n","").split("\t")
                sector = ll[0]
                npoligons += int(ll[1])
                start = int(ll[2])
                #start_datetime = "2017-09-01 00:00:00"
                start  = int(datetime.timestamp((datetime.strptime(start_datetime, '%Y-%m-%d %H:%M:%S') + timedelta(seconds=int(ll[2]))).replace(tzinfo=timezone('GMT'))))
                stop = int(ll[3])
                stop = int(datetime.timestamp((datetime.strptime(start_datetime, '%Y-%m-%d %H:%M:%S') + timedelta(seconds=int(ll[3]))).replace(tzinfo=timezone('GMT'))))
                cap = int(ll[4])
                if not sector in capacity:
                    capacity[sector] = {"capacity":dict(),"bounds":dict()}
                if not (start,stop) in capacity[sector]['capacity']:
                    capacity[sector]['capacity'][(start,stop)] = cap
    else:
        print("Manca il file con le capacità")
        exit(0)

ecac = None
with open(filename_poligoni_settori) as fin:
    for line in fin:
        if line[0]=='#':
            continue
        ll = line.replace("\n","").split("\t")
        #print(ll)
        sector = ll[0]
        #[start,stop] = [int(x) for x in ll[1].split(", ")]
        #print(ll[1])
        [start, stop] = [int(datetime.timestamp((datetime.strptime(start_datetime, '%Y-%m-%d %H:%M:%S') + timedelta(seconds=int(x))).replace(tzinfo=timezone('GMT')))) for x in ll[1].split(", ")]
        #print([start,stop])
        sub = ll[2]
        low,high = int(ll[3]),int(ll[4])
        if high >= QUOTECUT:
            punti = ll[5].replace("POLYGON ((","").replace(")","").replace(",","").split(" ")
            if sector in capacity:
                if not (low,high) in capacity[sector]["bounds"]:
                    capacity[sector]["bounds"][(low,high)] = dict()
                if not sub in capacity[sector]["bounds"][(low,high)]:
                    capacity[sector]["bounds"][(low,high)][sub] = {"points":[],"Polygon":None}
                    i=0
                    while(i<len(punti)):
                        lat = float(punti[i])
                        lon = float(punti[i+1])
                        capacity[sector]["bounds"][(low,high)][sub]["points"].append((lat,lon))
                        i += 2
                    capacity[sector]["bounds"][(low,high)][sub]["Polygon"] = Polygon(capacity[sector]["bounds"][(low,high)][sub]["points"])
                    if ecac==None:
                        ecac = capacity[sector]["bounds"][(low,high)][sub]["Polygon"]
                    else:
                        ecac = cascaded_union([ecac,capacity[sector]["bounds"][(low,high)][sub]["Polygon"]])
            else:
                print("Errore!!! controllare i file dei settori")
                exit(0)

#m = folium.Map(location=[45.5236, -122.6750])
#folium.PolyLine(ecac.exterior.coords, color="green", weight=2.5, opacity=1).add_to(m)
#folium.Circle([40.2000000000000028422,31.8230000000000003979],color="red").add_to(m)
#print(ecac.contains(Point(40.2000000000000028422,31.8230000000000003979)))
#m.save("./ecac.html")

with open("./config/bound_latlon.dat","w") as fout:
    for p in ecac.exterior.coords:
        fout.write(str(p[0])+"\t"+str(p[1])+"\n")
mappa_settore_numero = dict()
mappa_numero_settore = dict()
s = 0
import os
npoly_tot = 0
for sector in capacity:
    #print(sector)
    #print(capacity[sector])
    for (start,stop) in capacity[sector]['capacity']:
        npoly = 0
        for (l,h) in capacity[sector]["bounds"]:
            for sub in capacity[sector]["bounds"][(l,h)]:
                npoly += 1
        if npoly>0:
            cap = capacity[sector]['capacity'][(start,stop)]
            if not (sector,start,stop,cap) in mappa_settore_numero:
                s += 1
                mappa_settore_numero[(sector,start,stop,cap)] = s
                mappa_numero_settore[s] = (sector,start,stop,cap)
                with open(confpath+"boundary/"+str(s)+"_bound_latlon.dat","w") as fout:
                    fout.write(str(npoly)+"\n")
            #print("poligoni",npoly)
            with open(confpath+"boundary/"+str(s)+"_bound_latlon.dat","a") as fout:
                for (l,h) in capacity[sector]["bounds"]:
                    for sub in capacity[sector]["bounds"][(l,h)]:
                        ppp = ""
                        for p in capacity[sector]["bounds"][(l,h)][sub]['points']:
                            ppp +=str(p[0])+","+str(p[1])+"\t"
                        ppp = ppp[:-1]+"\n"
                        fout.write(str(l)+"\t"+str(h)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(len(capacity[sector]["bounds"][(l,h)][sub]['points']))+"\t"+ppp)
                        npoly_tot += 1
    #print(mappa_settore_numero)
    #print(mappa_numero_settore)
    #input("?")
with open(capacity_file,"w") as fout:
    fout.write("#Sector\tCapacity\n")
    for ss in sorted(mappa_numero_settore):
        #mappa_settore_numero[(sector,start,stop,cap)] = s
        fout.write(str(ss)+"\t"+str(mappa_numero_settore[ss][3])+"\n")
settore_nullo = s + 1
mappa_numero_settore[settore_nullo] = ("ECAC",0,0,9999)
pickle.dump(mappa_numero_settore,open("mappa_numero_settore.pp","wb"))
pickle.dump(mappa_settore_numero,open("mappa_settore_numero.pp","wb"))
print("numero settori",s)
print("numero poligoni",npoly_tot)
print("settore nullo",settore_nullo)
voli = dict()

new_fid = -1
with open(filename_so6_input.replace(".so6",".ids"),"w") as fout:
    fout.write("#so6 abm\n")
    with open(so6_folder+filename_so6_input) as f:
        old_fid = None
        for row in f:
            campi = row.strip("\n").strip("\r").split(" ")
            #segment = campi[0]
            segorigin = campi[0].split("_")[0]
            segdest = campi[0].split("_")[1]
            origin = campi[1]
            dest = campi[2]
            aircraft = campi[3]
            begintime = campi[4]
            endtime = campi[5]
            flbegin = float(campi[6])
            flend = float(campi[7])
            status = campi[8]
            callsign = campi[9]
            datestart = campi[10]
            datestop = campi[11]
            latitudebegin = "%.3f" % (float(campi[12])/60.)
            latitudebegin = float(latitudebegin)
            longitudebegin = "%.3f" % (float(campi[13])/60.)
            longitudebegin = float(longitudebegin)
            latitudeend = "%.3f" % (float(campi[14])/60.)
            latitudeend = float(latitudeend)
            longitudeend = "%.3f" % (float(campi[15])/60.)
            longitudeend = float(longitudeend)
            fid = campi[16]
            seq = int(campi[17])
            seg_len = campi[18]
            parity = campi[19]
            p1 = (latitudebegin,longitudebegin)
            p2 = (latitudeend,longitudeend)

            if old_fid!=fid:
                new_fid += 1
                old_fid = fid

            if not new_fid in voli:
                voli[new_fid] = {"fid":fid}
                fout.write(str(fid)+" "+str(new_fid)+"\n")

            if not seq in voli[new_fid]:
                if seq==1:
                    #{"Origin":segorigin,"Aircraft":aircraft,
                    voli[new_fid][0] = {"Time":begintime,"Date":datestart,"latitude":float(latitudebegin),"longitude":float(longitudebegin),"quota":flbegin}
                voli[new_fid][seq] = {"Time":endtime,"Date":datestop,"latitude":float(latitudeend),"longitude":float(longitudeend),"quota":flend}
print(len(voli),"voli caricati!")

#with open(confpath+"bound_latlon.dat","w") as fout:
#    for point in ecac:
#        fout.write(str(point[0])+"\t"+str(point[1])+"\n")

usati = dict()
punti_usati = dict()
Nflight = 0
print("so6")
with open(filename_so6_input+".abm","w") as fout:
    lines = ''
    for fid in voli:
        del voli[fid]['fid']
        sequences = list(sorted(voli[fid].keys()))
        line = ''
        pp = 0
        lat_lon = dict()

        for seq in sorted(sequences):
            if voli[fid][seq]["quota"]>QUOTECUT:
                #print(voli[fid][seq])
                #check for points inside ecac area (segment level)
                if not ecac.contains(Point(voli[fid][seq]['latitude'],voli[fid][seq]['longitude'])):
                    if seq+1 in sequences:
                        if not ecac.contains(Point(voli[fid][seq+1]['latitude'],voli[fid][seq+1]['longitude'])):
                            continue
                    else:
                        if seq-1 in sequences:
                            if not ecac.contains(Point(voli[fid][seq-1]['latitude'],voli[fid][seq-1]['longitude'])):
                                continue
                        else:
                            continue
                # here the segment is inside ecac
                # check for duplicate points on the route
                if (voli[fid][seq]['latitude'],voli[fid][seq]['longitude']) in lat_lon:
                    continue
                lat_lon[(voli[fid][seq]['latitude'],voli[fid][seq]['longitude'])] = None
                # updating "punti usati" list to avoid using duplicates point on random navps
                if not (voli[fid][seq]['latitude'],voli[fid][seq]['longitude']) in punti_usati:
                    punti_usati[(voli[fid][seq]['latitude'],voli[fid][seq]['longitude'])] = None

                #Add line to input
                pp +=1
                line += str(voli[fid][seq]['latitude'])+","+str(voli[fid][seq]['longitude'])
                line += ","+str(float(voli[fid][seq]['quota']))+","
                line += century+voli[fid][seq]['Date'][:2]+"-"+voli[fid][seq]['Date'][2:4]+"-"+voli[fid][seq]['Date'][4:]
                line += " "+voli[fid][seq]['Time'][:2]+":"+voli[fid][seq]["Time"][2:4]+":"+voli[fid][seq]["Time"][4:]
                point = Point(float(voli[fid][seq]['latitude']),float(voli[fid][seq]['longitude']))
                quota = voli[fid][seq]['quota']
                st = settore_nullo
                #print(voli[fid][seq]["Date"],voli[fid][seq]["Time"])
                if not voli[fid][seq]["Date"][4:] == start_datetime[8:10]:
                    st = 0
                else:
                    if ecac.contains(point):
                        #secondi = int(voli[fid][seq]["Time"][:2])*3600+int(voli[fid][seq]["Time"][2:4])*60+int(voli[fid][seq]["Time"][4:])
                        secondi = int(datetime.timestamp(datetime.strptime(voli[fid][seq]["Date"]+" "+voli[fid][seq]["Time"],"%y%m%d %H%M%S").replace(tzinfo=timezone('GMT'))))
                        if voli[fid][seq]['quota']>=QUOTECUT_sector:
                            for sector in capacity:
                                if st != settore_nullo:
                                    break
                                for (start,stop) in capacity[sector]['capacity']:
                                    if st != settore_nullo:
                                        break
                                    if start <= secondi  and secondi <stop:
                                        for (l,h) in capacity[sector]['bounds']:
                                            if st !=settore_nullo:
                                                break
                                            if l <=quota and quota < h:
                                                for sub in capacity[sector]['bounds'][(l,h)]:
                                                    if capacity[sector]['bounds'][(l,h)][sub]["Polygon"].contains(point):
                                                        st = mappa_settore_numero[sector,start,stop,capacity[sector]['capacity'][(start,stop)]]
                                                        #print(st)
                                                        break
                        else:
                            st=0
                    else:
                        st = 0

                line += ","+str(st)+"\t"
            #print(line)
        #print(line)
        #input("?")
        line = str(fid)+"\t"+str(pp)+"\t"+line+"\n"
        #print(line)
        if (pp>1):
            #fout.write(line)
            Nflight += 1
            lines += line
            if not fid in usati:
                usati[fid] = None
            #new_fid += 1
    fout.write(str(Nflight)+"\tNflight\n")
    fout.write(lines)

#pickle.dump(voli,open("voli_nuovo.pp","wb"))
nsim = 50
max_ang = 0.2745
extr_ang = 0.4745
direct_thr = 0.21275862069
x_capacity = 0.672413793103
rer_active = 1
ls=1200
as_=1.
max_t = 1200
xdelay =0
pdelay = 0
use_delay = 0
t_w = 45
t_d = 90
t_i = 10
t_r = 0.4
shortest_path = 1
d_thr = 10000
noise_d_thr = 10000
geom = 1
sig_V = 0
laplacian_vel = 0
Nm_shock = 0
radius = 18500
shock_f_lvl_min = 240
shock_f_lvl_max = 300
lifetime = 3
tmp_from_file = 1
if current_path == '':
    current_path = os.getcwd()






def old_delay():
    print("flight id usati",len(usati))
    with open(delay_file,"w") as fout:
        fout.write("#FlightID\tDelay\n")
        for fid in ritardi:
            if str(fid) in usati:
                fout.write(str(fid)+"\t"+ritardi[str(fid)]+"\n")
                continue
            if int(fid) in usati:
                fout.write(str(fid)+"\t"+str(ritardi[int(fid)])+"\n")

with open("./config/config.cfg","w") as fout:
    fout.write("# Configuration File. Attention each value has to be followed by the sharp with the label of variable\n\n")
    fout.write("# Number of simulation performed by the ABM\n")
    fout.write(str(nsim)+"\t#nsim\n\n")
    fout.write("# Maximum Angle of deviation from original trajectory in rerouting (rad)\n")
    fout.write("# and Extreame angle for deviation (rad)\n")
    fout.write(str(max_ang)+"\t#max_ang\n")
    fout.write(str(extr_ang)+"\t#extr_ang\n\n")
    fout.write("# Percentage of possibility to have a direct\n")
    fout.write(str(direct_thr)+"\t#direct_thr\n\n")
    fout.write("#A moltiplicative factor for capacity\n")
    fout.write(str(x_capacity)+"\t#x_capacity\n\n")
    fout.write("#To Activate the rerouting module (Boolean)\n")
    fout.write(str(rer_active)+"\t#rer_active\n\n")
    fout.write("# Minimum Improvement of a direct (meters)\n")
    fout.write(str(ls)+"\t#ls\n\n")
    fout.write("# Sensitivity Angle for direct (deg)\n")
    fout.write(str(as_)+"\t#as\n\n")
    fout.write("# Maximum reroute time for direct\n")
    fout.write(str(max_t)+"\t#max_t\n\n")
    fout.write("# Maximum amount of delay on departure (sec)\n")
    fout.write(str(xdelay)+"\t#xdelay\n\n")
    fout.write("# Percentage of flight with xdelay\n")
    fout.write(str(pdelay)+"\t#pdelay\n\n")
    fout.write("# Use external delay file: 1 Yes, 0 No\n")
    fout.write(str(use_delay)+"\t#use_delay\n\n")

    fout.write("# Number of elementary time increment in a time-step\n")
    fout.write(str(t_w)+"\t#t_w\n\n")
    fout.write("# Number of elementary time increment for direct\n")
    fout.write(str(t_d)+"\t#t_d\n\n")
    fout.write("# Size of a time incremet (sec)\n")
    fout.write(str(t_i)+"\t#t_i\n\n")
    fout.write("# Fraction of t_w after which the alghorithm is updated\n")
    fout.write(str(t_r)+"\t#t_r\n\n")
    fout.write("#Boolean 1) shortest path 0) minimum deviation (rerouting)\n")
    fout.write(str(shortest_path)+"\t#shortest_path\n\n")
    fout.write("#threshold value of the safety distance between aircraft (meters)\n")
    fout.write(str(d_thr)+"\t#d_thr\n\n")
    fout.write("#threshold value of the safety event at 15m (meters)\n")
    fout.write(str(noise_d_thr)+"\t#noise_d_thr\n\n")
    fout.write("#Boolean 1) Peter-Gall projection 2) Spheric Geometry\n")
    fout.write(str(geom)+"\t#geom\n\n")
    fout.write("#Width of distribution of noise on velocity. Needs to be between -1 and 1 (not included).\n")
    fout.write(str(sig_V)+"\t#sig_V\n\n")
    fout.write("# Boolean to have a laplacian variation of velocity\n")
    fout.write(str(laplacian_vel)+"\t#laplacian_vel\n\n")
    fout.write("# Average number of shock per time-step per flight level; (Unstable)\n")
    fout.write(str(Nm_shock)+"\t#Nm_shock\n\n")
    fout.write("# Radius of the shock (meters);  (Unstable)\n")
    fout.write(str(radius)+"\t#radius\n\n")
    fout.write("# Minimum and Maximum flight level for shocks;  (Unstable)\n")
    fout.write(str(shock_f_lvl_min)+"\t#shock_f_lvl_min\n")
    fout.write(str(shock_f_lvl_max)+"\t#shock_f_lvl_max\n\n")
    fout.write("# Average lifetime of a shock ( t_w*t_r*t_i unity );  (Unstable)\n")
    fout.write(str(lifetime)+"\t#lifetime\n\n")
    fout.write("# Boolean. If 1, new temporary navpoints are read from the disk. Otherwise they are generated. Remark: should always be set to 1! TODO: remove this.\n")
    fout.write(str(tmp_from_file)+"\t#tmp_from_file\n\n")
    fout.write("# Stating and Ending Datetime of the Simulation Year-Mounth-Day Hour:minute:second\n")
    fout.write(start_datetime+"\t#start_datetime\n")
    fout.write(end_datetime+"\t#end_datetime\n")
    fout.write("# Directories for boundaries. Should be fixed!!\n")
    fout.write(bound_file+"\t#bound_file\n\n")
    fout.write("# file for temporary points. Should be fixed!!\n")
    fout.write(temp_nvp+"\t#temp_nvp\n\n")
    fout.write("# file for shocks. Can be equal of temp_nvp but should be fixed!!\n")
    fout.write(shock_tmp+"\t#shock_tmp\n\n")
    fout.write("#file for capacity Should be fixed!!\n")
    fout.write(capacity_file+"\t#capacity_file\n\n")
    fout.write(str(npoly_tot+2)+"\t#npol\n\n")
    #fout.write("#file for maximum delay of each flight in minutes\n")
    #fout.write(delay_file+"\t#delay_file\n\n")



navp_file = 'temp_nvp.dat'
points = 0
if os.path.exists("./input/"+navp_file):
    with open("./input/"+navp_file) as fin:
        with open("./config/sectors_"+navp_file,"w") as fout:
            for line in fin:
                ll = line.replace("\n","").split("\t")
                lat = "%.3f" % float(ll[0])
                lat = float(lat)
                lon = "%.3f" % float(ll[1])
                lon = float(lon)
                if not (lat,lon) in punti_usati:
                    punti_usati[(lat,lon)] = None
                    p = Point(lat,lon)
                    if ecac.contains(p):
                        fout.write(str(lat)+"\t"+str(lon)+"\n")
                        points += 1
                if points==NPoints:
                    break
print("punti")
with open("./config/sectors_"+navp_file,"a") as fout:
    while (points<NPoints):
        #print(points,NPoints)
        lat = np.random.uniform(lat_min,lat_max)
        lat = "%.3f" % lat
        lat = float(lat)
        lon = np.random.uniform(lon_min,lon_max)
        lon = "%.3f" % lon
        lon = float(lon)
        if not (lat,lon) in punti_usati:
            punti_usati[(lat,lon)] = None
            p = Point(lat,lon)
            if ecac.contains(p):
                fout.write(str(lat)+"\t"+str(lon)+"\n")
                points += 1
