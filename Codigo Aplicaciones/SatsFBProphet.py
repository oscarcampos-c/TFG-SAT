import pandas as pd
import argparse
from fbprophet import Prophet

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--NumMonths', help='Number of months to predict. Default: 12',default = 12)
    parser.add_argument('--NumMonthsBack', help='Number of months to look back for the prediction. Default: 60',default = 60)
    parser.add_argument('--input', help='inputfile Default: SatStatistics.csv',default ='SatStatistics.csv')
    parser.add_argument('--output', help='output file path (will be overwritten). Default: SatPredictionPH.csv',default ='SatPredictionPH.csv')
    args=parser.parse_args()
    NumMonths = args.NumMonths
    outputfile = args.output
    inputfile = args.input
    NumMonthsBack = args.NumMonthsBack

    data = pd.read_csv(inputfile,parse_dates=['MesAnno'])
    data['ds'] = data['MesAnno']
    data['y'] = data['NumSats']
    data.drop(labels=['MesAnno','NumSats'],axis=1,inplace=True)
    data =data.head(NumMonthsBack)

    m = Prophet()
    m.fit(data)

    future = m.make_future_dataframe(periods=int(NumMonths),freq='MS')
    fcst = m.predict(future)
    
    result = fcst[~fcst.ds.isin(data.ds)]
    result.drop(result.columns.difference(['ds','trend']), 1, inplace=True)
    result = pd.merge(data,
            result[['ds', 'trend']],
            on='ds', 
            how='outer')
    result.sort_values(by=['ds'],inplace=True)
    result.columns = ['MesAnno', 'NumSats','Prediction']
    result.to_csv(outputfile,index=False)
