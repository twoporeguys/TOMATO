import scipy.io as sio
import numpy as np
import pandas as pd
import os
import sys


def data_append(file_names):
    dwell = []
    meanAmp = []
    maxAmp = []
    sdAmp = []
    setNum = []
    medAmp = []
    SNR = []
    PoreNum = []
    Exp = []
    NewArea = []

    for file_name in file_names:
        mat_contents = sio.loadmat(file_name)
        data = mat_contents['eventDataPassFour']
        if int(data['PostTestingEventTotal'][0][0][0]) > 0:
            dwell.extend(np.log10(data['dwell'][0][0][0]))
            meanAmp.extend(data['meanAmp'][0][0][0])
            maxAmp.extend(data['maxAmp'][0][0][0])
            sdAmp.extend(data['sdAmp'][0][0][0])
            setNum.extend(int(data['PostTestingEventTotal'][0][0][0])*[data['setNum'][0][0][0]])
            medAmp.extend(data['medAmp'][0][0][0])
            SNR.extend(data['SNR'][0][0][0])
            PoreNum.extend(int(data['PostTestingEventTotal'][0][0][0])*[data['PoreNum'][0][0][0]])
            Exp.extend(int(data['PostTestingEventTotal'][0][0][0])*[data['Exp'][0][0][0]])
            NewArea.extend(data['NewArea'][0][0][0])

    data_extracted = np.column_stack((dwell, meanAmp, maxAmp,
                                      sdAmp, setNum, medAmp, SNR,
                                      PoreNum, Exp, NewArea))
    df = pd.DataFrame(data_extracted, columns=('dwell',
                                               'meanAmp', 'maxAmp', 'sdAmp', 'setNum',
                                               'medAmp', 'SNR', 'PoreNum', 'Exp',
                                               'NewArea'))
    return df


def load_file_names(folder_name):
    files_list = []
    for file in os.listdir(folder_name):
        if file.endswith(".mat"):
            files_list.append(os.path.join(folder_name, file))
    files_list.sort()
    return files_list


if __name__ == '__main__':
    loc = sys.argv[1]
    save_loc = sys.argv[2]
    file_names = load_file_names(loc)
    data = data_append(file_names)
    data.to_csv(save_loc+'/all_events.csv')
