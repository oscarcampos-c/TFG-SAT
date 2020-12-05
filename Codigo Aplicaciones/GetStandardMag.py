# -*- coding: utf-8 -*-

import requests
import pandas as pd
import argparse

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--outputfile', help='Output file path',default = 'Satflare.csv')
    args=parser.parse_args()
    outputfile = args.outputfile

    satlafrebase = 'Satflare'
    url = 'http://www.satflare.com/advsearch.aspx?sat=&ExcDeb=1&ExcRB=1&max=30000'
    html = requests.get(url).content
    open("{0}.html".format(satlafrebase), 'wb').write(html)
    print('list downloaded')
    raw = open("{0}.html".format(satlafrebase), "r")
    lines = raw.readlines()
    raw.close()
    startadding =False
    stopadding = False
    new = open("{0}Mod.html".format(satlafrebase), "w")
    for line in lines:
        if (startadding and not stopadding):
            new.write(line)      
        if line.find('TRACK ALL THESE OBJECTS') > 0:
            startadding = True
        if line.find('RA and DEC are J2000') > 0:
            stopadding = True  
        
    new.close()

    html_file = open("{0}Mod.html".format(satlafrebase),mode="r").read()
    df_list = pd.read_html(html_file, decimal=',',thousands=None)
    df = df_list[-1]
    headers = df.iloc[0]
    df  = pd.DataFrame(df.values[1:], columns=headers)
    df.to_csv(outputfile,index=False)

