import csv, os
import pandas as pd


file_path = '/Users/Kimokeo Bowden Jr/Desktop/capstone/armstrong/src/data_exp_63315-v4_task-rrep.csv'
# file_path = '/Users/danieldang/code/armstrong/src/data_exp_63315-v4_task-rrep.csv'

def getVariablesFromCSV(file_path):
    variables = {}
    df = pd.read_csv(file_path)
    # with open(file_path, "r") as f:
    #     dialect = csv.Sniffer().sniff(f.read(10240))
    #     f.seek(0)
    #     reader = csv.reader(f, dialect)
    #     header = next(reader)
    startRow, endRow = getStartAndFinishRowsOfFirstID(df)
    print(getConditions(df, startRow, endRow))
    #     print(f'startRow: {startRow} and endRow: {endRow}')
    #     stimuli = getStimuli(reader, header, startRow, endRow)
    #     print(f'stimuli: {stimuli}')
    #     printValAtStartRow(reader, header, startRow)
    #     nTrials = getTrialNumber(reader, header, endRow)
    #     nTrialsPerStim = getNTrialsPerStim(nTrials, stimuli)

def getStartAndFinishRowsOfFirstID(df):
    startRow = None
    endRow = None

    for i, line in df.iterrows():
        if line["display"] == 'trial' and startRow is None:
            startRow = i
            print(f'first time we see trial: {line["display"]} at i = {i}')
        if line["Trial Number"] == 'END TASK' and endRow is None:
            endRow = i - 1
            break
    return startRow, endRow

def getConditions(df, startRow, endRow):
    conditions = []
    for i, line in df.iterrows():
        if i < startRow or type(line['condition']) == float:
             continue
        if i == endRow: break
        condition = line['condition']
        if condition not in conditions:
            conditions.append(condition)
    return conditions


    
getVariablesFromCSV(file_path)



