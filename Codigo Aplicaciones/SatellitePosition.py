import pandas as pd
from skyfield.api import EarthSatellite,Topos, load
import math
from skyfield.constants import (ERAD)
import time,sys
import datetime
from multiprocessing import Pool,cpu_count,freeze_support
import tqdm
import pandas
import argparse

def fillValues(P_SATID,P_CITY,P_UTCCALCDATE,P_TLEEEPOCH,P_UTC0R,P_UTC0S,P_UTC15R,P_UTC15S,P_UTC30R,P_UTC30S,P_UTCMAX,
               P_MAXELEV,P_SUNELEVAT0R,P_SUNELEVAT0S,P_SUNELEVAT15R,P_SUNELEVAT15S,P_SUNELEVAT30R,P_SUNELEVAT30S,P_SUNELEVATMAX,P_ISSUNLIT0R,
               P_ISSUNLIT0S,P_ISSUNLIT15R,P_ISSUNLIT15S,P_ISSUNLIT30R,P_ISSUNLIT30S,P_ISSUNLITMAX,P_UTCSHADOW,P_ELEVATSHADOW,
               P_MAG0R,P_MAG0S,P_MAG15R,P_MAG15S,P_MAG30R,P_MAG30S,P_MAGMAX,P_MAGPRESHADOW,P_AZ0R,P_AZ15R,P_AZ30R,P_AZMAX,P_AZ30S,P_AZ15S,P_AZ0S):

    NewObject = []
    NewObject.append(P_SATID)
    NewObject.append(P_CITY)
    NewObject.append(str(P_UTCCALCDATE.utc_iso(' ')).replace("Z",""))
    NewObject.append(P_TLEEEPOCH)
    NewObject.append(str(P_UTC0R.utc_iso (' ')).replace("Z","") if P_UTC0R  is not None else None)
    NewObject.append(str(P_UTC0S.utc_iso (' ')).replace("Z","") if P_UTC0S  is not None else None)
    NewObject.append(str(P_UTC15R.utc_iso(' ')).replace("Z","") if P_UTC15R is not None else None)
    NewObject.append(str(P_UTC15S.utc_iso(' ')).replace("Z","") if P_UTC15S is not None else None)
    NewObject.append(str(P_UTC30R.utc_iso(' ')).replace("Z","") if P_UTC30R is not None else None)
    NewObject.append(str(P_UTC30S.utc_iso(' ')).replace("Z","") if P_UTC30S is not None else None)
    NewObject.append(str(P_UTCMAX.utc_iso(' ')).replace("Z","") if P_UTCMAX is not None else None)
    NewObject.append(P_MAXELEV)
    NewObject.append(P_SUNELEVAT0R)
    NewObject.append(P_SUNELEVAT0S)
    NewObject.append(P_SUNELEVAT15R)
    NewObject.append(P_SUNELEVAT15S)
    NewObject.append(P_SUNELEVAT30R)
    NewObject.append(P_SUNELEVAT30S)
    NewObject.append(P_SUNELEVATMAX)
    NewObject.append(P_ISSUNLIT0R)
    NewObject.append(P_ISSUNLIT0S)
    NewObject.append(P_ISSUNLIT15R)
    NewObject.append(P_ISSUNLIT15S)
    NewObject.append(P_ISSUNLIT30R)
    NewObject.append(P_ISSUNLIT30S)
    NewObject.append(P_ISSUNLITMAX)
    NewObject.append(str(P_UTCSHADOW.utc_iso(' ')).replace("Z","") if P_UTCSHADOW is not None else None)
    NewObject.append(P_ELEVATSHADOW)
    NewObject.append(P_MAG0R)
    NewObject.append(P_MAG0S)
    NewObject.append(P_MAG15R)
    NewObject.append(P_MAG15S)
    NewObject.append(P_MAG30R)
    NewObject.append(P_MAG30S)
    NewObject.append(P_MAGMAX)
    NewObject.append(P_MAGPRESHADOW)
    NewObject.append(P_AZ0R)    
    NewObject.append(P_AZ15R)
    NewObject.append(P_AZ30R)
    NewObject.append(P_AZMAX)
    NewObject.append(P_AZ30S)
    NewObject.append(P_AZ15S)
    NewObject.append(P_AZ0S)

    return(NewObject)

def calculateMag(t,SatToposDiff,SunEarth,EarthLocSun,stdmag):
    topocentric = SatToposDiff.at(t)
    alt, az, distance = topocentric.altaz()
    SunEarthDistance = SunEarth.at(t).distance().km
    LocSunAngle = EarthLocSun 
    
    a = SunEarthDistance - ERAD/1000 #distance sun from observer (Km)
    b = distance.km # distance to ISS from observer (Km)
    angle_c = topocentric.separation_from(LocSunAngle).degrees * 0.0174533
    c = math.sqrt( math.pow(a,2) + math.pow(b,2) - 2*a*b*math.cos( angle_c) )
    angle_a = math.acos((math.pow(b,2) + math.pow( c,2) - math.pow(a,2)) / (2 * b * c)) 
    phase_angle = angle_a * 57.2958 

    distanceToSatellite =  topocentric.distance().km #This is in KM
    phaseAngleDegrees = phase_angle #Angle from sun->satellite->observer
    pa = phaseAngleDegrees * 0.0174533 #Convert the phase angle to radians
    intrinsicMagnitude = stdmag #also called standard magnitude
    term_1 = intrinsicMagnitude
    term_2 = 5.0 * math.log10(distanceToSatellite / 1000) # or 2.5*LOG((distanceToSatellite/1000)^2) 
    arg = math.sin(pa) + (math.pi - pa) * math.cos(pa)
    term_3 = -2.5 * math.log10(arg)
    apparentMagnitude = term_1 + term_2 + term_3
    return apparentMagnitude


def calculatePasses(params):
    try:
        ephem = load(params[13])
        sun = ephem['sun']
        earth = ephem['earth']
        CityElevation = int(params[15])
        Location = Topos(params[12].split(',')[0],params[12].split(',')[1],elevation_m=CityElevation)
        
        tle0=params[0]
        tle1=params[1]
        tle2=params[2]
        NORADID=params[3]
        tleEpoch=params[4]
        City=params[5]
        YearFrom=params[6]
        MonthFrom=params[7]
        DayFrom=params[8]
        YearTo=params[9]
        MonthTo=params[10]
        DayTo=params[11]
        stdmag = params[14]
        ts = load.timescale()    
        E_TLEEEPOCH         = tleEpoch
        E_SATID             = NORADID
        E_CITY              = City
        E_UTCCALCDATE       = ts.now()
        E_TLEEEPOCH         = tleEpoch
        satellite           = EarthSatellite(tle1, tle2, tle0, ts)
        difference          = satellite - Location
        EarthLoc            = (earth + Location)
        SunEarth            = (sun - earth)
        obj =[]
        t0 = ts.utc(YearFrom, MonthFrom, DayFrom)
        t1 = ts.utc(YearTo, MonthTo, DayTo)
    
        lastevent = -1   
        times, events = satellite.find_events(Location, t0, t1, altitude_degrees=0)
    
        if len(events)==0:
            return obj
        
        if events[len(events)-1]!=2:
            t1 = ts.utc(YearTo, MonthTo, DayTo,0,40,0)
            times, events = satellite.find_events(Location, t0, t1, altitude_degrees=0)
        
        for ti, event in zip(times, events):
            if lastevent==-1 and event>0:
                continue
            else:
                lastevent=event
            if event==0:
                E_UTC0R             =None
                E_UTC0S             =None
                E_UTC15R            =None
                E_UTC15S            =None
                E_UTC30R            =None
                E_UTC30S            =None
                E_UTCMAX            =None
                E_MAXELEV           =None
                E_SUNELEVAT0R       =None
                E_SUNELEVAT0S       =None
                E_SUNELEVAT15R      =None
                E_SUNELEVAT15S      =None
                E_SUNELEVAT30R      =None
                E_SUNELEVAT30S      =None
                E_SUNELEVATMAX      =None
                E_ISSUNLIT0R        =None
                E_ISSUNLIT0S        =None
                E_ISSUNLIT15R       =None
                E_ISSUNLIT15S       =None
                E_ISSUNLIT30R       =None
                E_ISSUNLIT30S       =None
                E_ISSUNLITMAX       =None
                E_UTCSHADOW         =None
                E_ELEVATSHADOW      =None
                E_MAG0R             =None
                E_MAG0S             =None
                E_MAG15R            =None
                E_MAG15S            =None
                E_MAG30R            =None
                E_MAG30S            =None
                E_MAGMAX            =None
                E_MAGPRESHADOW      =None
                E_AZ0R              =None
                E_AZ15R             =None
                E_AZ30R             =None
                E_AZMAX             =None
                E_AZ30S             =None
                E_AZ15S             =None
                E_AZ0S              =None
                
                E_UTC0R = ti
                topocentric = difference.at(E_UTC0R)
                alt, az, distance = topocentric.altaz()
                E_AZ0R = az.degrees
                E_ISSUNLIT0R = satellite.at(E_UTC0R).is_sunlit(ephem)
                EarthLocSun = EarthLoc.at(E_UTC0R).observe(sun)
                altSun,azSun,distanceSun = EarthLocSun.apparent().altaz()
                E_SUNELEVAT0R = altSun.degrees
                if E_ISSUNLIT0R:
                    E_MAG0R = calculateMag(E_UTC0R,difference,SunEarth,EarthLocSun,stdmag)
            if event ==1:
                E_UTCMAX = ti
                topocentric = difference.at(E_UTCMAX)
                alt, az, distance = topocentric.altaz()
                E_MAXELEV = alt.degrees
                E_AZMAX = az.degrees
                E_ISSUNLITMAX = satellite.at(E_UTCMAX).is_sunlit(ephem)
                EarthLocSun =  EarthLoc.at(E_UTCMAX).observe(sun)
                altSun,azSun,distanceSun = EarthLocSun.apparent().altaz()
                E_SUNELEVATMAX = altSun.degrees
                if E_ISSUNLITMAX:
                    E_MAGMAX = calculateMag(E_UTCMAX,difference,SunEarth,EarthLocSun,stdmag)
            if event ==2:
                E_UTC0S = ti
                topocentric = difference.at(E_UTC0S)
                alt, az, distance = topocentric.altaz()
                E_AZ0S = az.degrees
                E_ISSUNLIT0S = satellite.at(E_UTC0S).is_sunlit(ephem)
                EarthLocSun = EarthLoc.at(E_UTC0S).observe(sun)
                altSun,azSun,distanceSun = EarthLocSun.apparent().altaz()
                E_SUNELEVAT0S = altSun.degrees
                if E_ISSUNLIT0S:
                    E_MAG0S = calculateMag(E_UTC0S,difference,SunEarth,EarthLocSun,stdmag)
                if E_MAXELEV >15:
                    times15, events15 = satellite.find_events(Location, E_UTC0R, E_UTC0S, altitude_degrees=15)
                    for ti15, event15 in zip(times15, events15): 
                        if event15==0: 
                            E_UTC15R = ti15
                            topocentric = difference.at(E_UTC15R)
                            alt, az, distance = topocentric.altaz()
                            E_AZ15R = az.degrees
                            E_ISSUNLIT15R = satellite.at(E_UTC15R).is_sunlit(ephem)
                            EarthLocSun=EarthLoc.at(E_UTC15R).observe(sun)
                            altSun,azSun,distanceSun = EarthLocSun.apparent().altaz()
                            E_SUNELEVAT15R = altSun.degrees
                            if E_ISSUNLIT15R:
                                E_MAG15R = calculateMag(E_UTC15R,difference,SunEarth,EarthLocSun,stdmag)
                        if event15==2:
                            E_UTC15S = ti15
                            topocentric = difference.at(E_UTC15S)
                            alt, az, distance = topocentric.altaz()
                            E_AZ15S = az.degrees
                            E_ISSUNLIT15S = satellite.at(E_UTC15S).is_sunlit(ephem)
                            EarthLocSun=EarthLoc.at(E_UTC15S).observe(sun)
                            altSun,azSun,distanceSun = EarthLocSun.apparent().altaz()
                            E_SUNELEVAT15S = altSun.degrees
                            if E_ISSUNLIT15S:
                                E_MAG15S = calculateMag(E_UTC15S,difference,SunEarth,EarthLocSun,stdmag)
                if E_MAXELEV >30:
                    times30, events30 = satellite.find_events(Location, E_UTC15R, E_UTC15S, altitude_degrees=30)
                    for ti30, event30 in zip(times30, events30): 
                        if event30 == 0:
                            E_UTC30R = ti30
                            topocentric = difference.at(E_UTC30R)
                            alt, az, distance = topocentric.altaz()
                            E_AZ30R = az.degrees
                            E_ISSUNLIT30R = satellite.at(E_UTC30R).is_sunlit(ephem)
                            EarthLocSun = EarthLoc.at(E_UTC30R).observe(sun)
                            altSun,azSun,distanceSun = EarthLocSun.apparent().altaz()
                            E_SUNELEVAT30R = altSun.degrees
                            if E_ISSUNLIT30R:
                                E_MAG30R = calculateMag(E_UTC30R,difference,SunEarth,EarthLocSun,stdmag)
                        if event30 == 2:
                            E_UTC30S = ti30
                            topocentric = difference.at(E_UTC30S)
                            alt, az, distance = topocentric.altaz()
                            E_AZ30S = az.degrees
                            E_ISSUNLIT30S = satellite.at(E_UTC30S).is_sunlit(ephem)
                            EarthLocSun = EarthLoc.at(E_UTC30S).observe(sun)
                            altSun,azSun,distanceSun = EarthLocSun.apparent().altaz()
                            E_SUNELEVAT30S = altSun.degrees
                            if E_ISSUNLIT30S:
                                E_MAG30S = calculateMag(E_UTC30S,difference,SunEarth,EarthLocSun,stdmag)
                
                SHADOWDIRECTION =0
                if (E_ISSUNLIT0R and (not E_ISSUNLITMAX or not E_ISSUNLIT0S)):
                    SHADOWDIRECTION = 1
                
                if (not E_ISSUNLIT0R and (E_ISSUNLITMAX or E_ISSUNLIT0S)):
                    SHADOWDIRECTION = -1
                
                if SHADOWDIRECTION !=0:
                    SatRiseTime = datetime.datetime.strptime(E_UTC0R.utc_iso(' '),'%Y-%m-%d %H:%M:%SZ')
                    SatSetTime = datetime.datetime.strptime(E_UTC0S.utc_iso(' '),'%Y-%m-%d %H:%M:%SZ')
                    SatTimeDiff = (SatSetTime-SatRiseTime).total_seconds()
                    TimeIncrease = 1
    
                    if SHADOWDIRECTION>0:
                        MinLitTime = SatRiseTime
                        MaxLitTime = SatSetTime
                        TTC1= SatRiseTime
                    if SHADOWDIRECTION<0:
                        MinLitTime = SatRiseTime
                        MaxLitTime = SatSetTime
                        TTC1= SatSetTime
                    
                    while TimeIncrease>0:
                        TimeIncrease = math.ceil(SatTimeDiff/2)
                        TTC1 = TTC1 + datetime.timedelta(seconds= SHADOWDIRECTION * TimeIncrease)
                        TTC2 = TTC1 + datetime.timedelta(seconds= SHADOWDIRECTION)
                        tsunlit1=ts.utc(TTC1.year,TTC1.month,TTC1.day,TTC1.hour,TTC1.minute,TTC1.second) 
                        tsunlit2=ts.utc(TTC2.year,TTC2.month,TTC2.day,TTC2.hour,TTC2.minute,TTC2.second) 
                        ShadowSunLit1 = satellite.at(tsunlit1).is_sunlit(ephem)
                        ShadowSunLit2 = satellite.at(tsunlit2).is_sunlit(ephem)
                        if ShadowSunLit1 != ShadowSunLit2:
                            TimeIncrease = -1
                        else:
                            if ShadowSunLit1 and SHADOWDIRECTION > 0:
                                MinLitTime = TTC1
                                
                            elif not ShadowSunLit1 and SHADOWDIRECTION > 0:
                                MaxLitTime = TTC1
                                SHADOWDIRECTION = -1
                                
                            elif ShadowSunLit1 and SHADOWDIRECTION <0:
                                MinLitTime = TTC1
                                SHADOWDIRECTION=1
                                
                            elif not ShadowSunLit1 and SHADOWDIRECTION <0:
                                MinLitTime=TTC1
                             
                            SatTimeDiff = (MaxLitTime-MinLitTime).total_seconds()
                            if SatTimeDiff <=1:
                                TimeIncrease=-1
                    
                    topocentric = difference.at(tsunlit1)
                    alt, az, distance = topocentric.altaz()                    
                    E_ELEVATSHADOW=alt.degrees
                    E_UTCSHADOW = tsunlit1
                    E_MAGPRESHADOW = calculateMag(E_UTCSHADOW,difference,SunEarth,EarthLocSun,stdmag)
    
                obj.append(fillValues( P_SATID=E_SATID,P_CITY=E_CITY,P_UTCCALCDATE=E_UTCCALCDATE,P_TLEEEPOCH=E_TLEEEPOCH,P_UTC0R=E_UTC0R,
                            P_UTC0S=E_UTC0S,P_UTC15R=E_UTC15R,P_UTC15S=E_UTC15S,P_UTC30R=E_UTC30R,P_UTC30S=E_UTC30S,P_UTCMAX=E_UTCMAX,P_MAXELEV=E_MAXELEV,
                            P_SUNELEVAT0R=E_SUNELEVAT0R,P_SUNELEVAT0S=E_SUNELEVAT0S,P_SUNELEVAT15R=E_SUNELEVAT15R,P_SUNELEVAT15S=E_SUNELEVAT15S,P_SUNELEVAT30R=E_SUNELEVAT30R,
                            P_SUNELEVAT30S=E_SUNELEVAT30S,P_SUNELEVATMAX=E_SUNELEVATMAX,P_ISSUNLIT0R=E_ISSUNLIT0R,P_ISSUNLIT0S=E_ISSUNLIT0S,P_ISSUNLIT15R=E_ISSUNLIT15R,
                            P_ISSUNLIT15S=E_ISSUNLIT15S,P_ISSUNLIT30R=E_ISSUNLIT30R,P_ISSUNLIT30S=E_ISSUNLIT30S,P_ISSUNLITMAX=E_ISSUNLITMAX,P_UTCSHADOW=E_UTCSHADOW,
                            P_ELEVATSHADOW=E_ELEVATSHADOW,P_MAG0R=E_MAG0R,P_MAG0S=E_MAG0S,P_MAG15R=E_MAG15R,P_MAG15S=E_MAG15S,P_MAG30R=E_MAG30R,P_MAG30S=E_MAG30S,
                            P_MAGMAX=E_MAGMAX,P_MAGPRESHADOW=E_MAGPRESHADOW,P_AZ0R=E_AZ0R,P_AZ15R=E_AZ15R,P_AZ30R=E_AZ30R,P_AZMAX=E_AZMAX,P_AZ30S=E_AZ30S,P_AZ15S=E_AZ15S,P_AZ0S=E_AZ0S)
                         )
        return obj
    except Exception as e:
        if hasattr(e, 'message'):
            return "Error {0} processing values: {1}".format(e.message,",".join(str(x) for x in params))
        else:
            return "Error {0} processing values: {1}".format(str(e),",".join(str(x) for x in params))

if __name__ == '__main__':
    freeze_support()
    parser=argparse.ArgumentParser()
    parser.add_argument('--tle', help='TLE csv file (must contain EPOCH,NORAD_CAT_ID,TLE_LINE0,TLE_LINE1,TLE_LINE2,STDMAG [optional]). Default: tle.csv',default = 'tle.csv')
    parser.add_argument('--output', help='output file path (will be overwritten). Default: satevents.csv',default ='satevents.csv')
    parser.add_argument('--errors', help='output file path for errors (will be overwritten). Default sateventserrors.txt',default ='sateventserrors.txt')
    parser.add_argument('--ephem', help='ephem file location. Default: de421.bsp',default ='de421.bsp')
    parser.add_argument('--CityName', help='Name of the city. Default None',default =None)
    parser.add_argument('--CityCoordinates', help='Name of the city. Default None',default =None)
    parser.add_argument('--DayFrom', help='Format YYYY-MM-DD',default =None)
    parser.add_argument('--DayTo', help='Format YYYY-MM-DD',default =None)
    parser.add_argument('--CityElevation', help='Meters above sea level. Default 0',default =0)

    args=parser.parse_args()
    TleFile = args.tle
    OutputFile = args.output
    CityName = args.CityName
    CityCoordinates = args.CityCoordinates
    DayFromParam = args.DayFrom
    DayToParam = args.DayTo
    CityElevation = args.CityElevation
    ErrorFile = args.errors
    if CityName is None or CityCoordinates is None or DayFromParam is None or DayToParam is None:
        print(parser.format_help())
        sys.exit()

    EventArray=[]
    tlecount=0
    TLEFound = 0
    start_time = time.time()
    paramToFunction = []
    AllTLE = pd.read_csv(TleFile,header=0)
    tlecount = 0
    TLEFound = len(AllTLE.index)
    
    YearFrom = int(DayFromParam.split('-')[0])
    MonthFrom = int(DayFromParam.split('-')[1])
    DayFrom = int(DayFromParam.split('-')[2])
    YearTo = int(DayToParam.split('-')[0])
    MonthTo = int(DayToParam.split('-')[1])
    DayTo = int(DayToParam.split('-')[2])
    	
    print ('using ',cpu_count(), ' cores')
    for index,row in AllTLE.iterrows():
        Params = []
        tlecount+=1
        tle0 = row["TLE_LINE0"]
        tle1 = row["TLE_LINE1"]
        tle2 = row["TLE_LINE2"]
        tleEpoch = row["EPOCH"]
        NORADID = row["NORAD_CAT_ID"]
        if 'STDMAG' in AllTLE:
            stdmag = row["STDMAG"]
        else:
            stdmag = 5 #row["TLE_LINE2"]
        
        Params.append(tle0)
        Params.append(tle1)
        Params.append(tle2)
        Params.append(NORADID)
        Params.append(tleEpoch)
        Params.append(CityName)
        Params.append(YearFrom)
        Params.append(MonthFrom)
        Params.append(DayFrom)
        Params.append(YearTo)
        Params.append(MonthTo)
        Params.append(DayTo)
        Params.append(CityCoordinates)
        Params.append(args.ephem)
        Params.append(stdmag)
        Params.append(CityElevation)
        paramToFunction.append(Params)

    if len(paramToFunction)>0:
        pool = Pool()
        function_output = pool.imap_unordered(calculatePasses, paramToFunction)
        pool.close()
        output = []
        errors = []
        for _ in tqdm.tqdm(function_output, total=len(paramToFunction)):
            if isinstance(_, str) :
                errors.append(_)
            else:
                output.append(_)
        pool.join()
        print("--- %s seconds calculating orbits ---" % ((time.time() - start_time)))
        if errors is not None:
            if len(errors)>0:
                with open(ErrorFile, 'w') as f:
                    for item in errors:
                        f.write("%s\n" % item)
                print('The process encountered errors processing TLEs, check the error output file. Errors in TLE: ', len(errors))
        EventArray= [item for sublist in output for item in sublist]
        df = pd.DataFrame(data=EventArray,
                  columns= ["SATID",
                            "CITY",
                            "UTCCALCDATE",
                            "TLEEEPOCH",
                            "UTC0R",
                            "UTC0S",
                            "UTC15R",
                            "UTC15S",
                            "UTC30R",
                            "UTC30S",
                            "UTCMAX",
                            "MAXELEV",
                            "SUNELEVAT0R",
                            "SUNELEVAT0S",
                            "SUNELEVAT15R",
                            "SUNELEVAT15S",
                            "SUNELEVAT30R",
                            "SUNELEVAT30S",
                            "SUNELEVATMAX",
                            "ISSUNLIT0R",
                            "ISSUNLIT0S",
                            "ISSUNLIT15R",
                            "ISSUNLIT15S",
                            "ISSUNLIT30R",
                            "ISSUNLIT30S",
                            "ISSUNLITMAX",
                            "UTCSHADOW",
                            "ELEVATSHADOW",
                            "MAG0R",
                            "MAG0S",
                            "MAG15R",
                            "MAG15S",
                            "MAG30R",
                            "MAG30S",
                            "MAGMAX",
                            "MAGPRESHADOW",
                            "AZ0R",
                            "AZ15R",
                            "AZ30R",
                            "AZMAX",
                            "AZ30S",
                            "AZ15S",
                            "AZ0S",
                            ])
							
    
    date = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    df.to_csv(OutputFile,index=False,header=True)
    print("--- %s seconds ---" % ((time.time() - start_time)))
   
