import csv, os
import pandas as pd


file_path = '/Users/Kimokeo Bowden Jr/Desktop/capstone/armstrong/src/data_exp_63315-v4_task-rrep.csv'
# file_path = '/Users/danieldang/code/armstrong/src/data_exp_63315-v4_task-rrep.csv'

def getVariablesFromCSV(file_path):
    variables = {}
    df = pd.read_csv(file_path)
    startRow, endRow = getStartAndFinishRowsOfFirstID(df)
    print(getConditions(df, startRow, endRow))
    print(f'startRow: {startRow} and endRow: {endRow}')
    stimuli = getStimuli(df, startRow, endRow)
    print(f'stimuli: {stimuli}')
    printValAtStartRow(df, startRow)
    nTrials = getTrialNumber(df, endRow)
    nTrialsPerStim = getNTrialsPerStim(nTrials, stimuli)
    print(f'nTrials : {nTrials} \n nTrialsPerStim : {nTrialsPerStim}')

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

def getNTrialsPerStim(nTrials, stimuli):
    randomStimulus = list(stimuli.keys())[0]
    numStimuli = len(stimuli[randomStimulus])
    return nTrials // numStimuli

def getTrialNumber(df, endRow):
    responseCount = 0
    for i, line in df.iterrows():
        if i == endRow:
            break
        if line['Response'] == '+':
            responseCount += 1
    return responseCount - 1

def printValAtStartRow(df, startRow):
    print(f'startRow in printValAtStartRow: {startRow}')
    for i, line in df.iterrows():
        print("i: " , i , ',' , "value: " , line['display'])
        if i == startRow:
            print("val at startRow: " + line['display'])
            break


def getStimuli(df, startRow, endRow):
    conditions = getConditions(df, startRow, endRow)
    stimuli = {}

    if len(conditions) == 1:
        onlyCondition = conditions[0]
        for i, line in df.iterrows():
            if type(line['left_stim']) is float:
                continue
            if i >= startRow and i <= endRow:
                indexUpToExtension = line['left_stim'].find('.') - 1 
                value = line['left_stim'][:indexUpToExtension + 1]
                if onlyCondition not in stimuli:
                    stimuli[onlyCondition] = [value]
                else:
                    if value not in stimuli[onlyCondition]:
                        stimuli[onlyCondition].append(value)
    
    elif len(conditions) > 1:
        for i, line in df.iterrows():
            if i >= startRow and i <= endRow:
                if type(line['left_stim']) is float:
                    continue
                indexUpToExtension = line['left_stim'].find('.') - 1            # -1 assuming N is < 10
                possibleKey = line['left_stim'][:indexUpToExtension]
                if possibleKey == '': continue
                value = line['left_stim'][:indexUpToExtension + 1]
                if possibleKey not in stimuli:
                    stimuli[possibleKey] = [value]
                else:
                    if value not in stimuli[possibleKey]:
                        stimuli[possibleKey].append(value)
            elif i > endRow:
                break
    
    return stimuli


    
getVariablesFromCSV(file_path)



