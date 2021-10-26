import csv, os

file_path = '/Users/danieldang/code/armstrong/src/data_exp_63315-v4_task-rrep.csv'

def getVariablesFromCSV(file_path):
    variables = {}
    with open(file_path, "r") as f:
        dialect = csv.Sniffer().sniff(f.read(10240))
        f.seek(0)
        reader = csv.reader(f, dialect)
        header = next(reader)
        startRow, endRow = getStartAndFinishRowsOfFirstID(reader, header)
        # print(f'startRow: {startRow} and endRow: {endRow}')
        stimuli = getStimuli(reader, header, startRow, endRow)
        # print(f'stimuli: {stimuli}')
        printValAtStartRow(reader, header, startRow)
        # nTrials = getTrialNumber(reader, header, endRow)
        # nTrialsPerStim = getNTrialsPerStim(nTrials, stimuli)

def getNTrialsPerStim(nTrials, stimuli):
    randomStimulus = list(stimuli.keys())[0]
    numStimuli = len(stimuli[randomStimulus])
    return nTrials // numStimuli

def getTrialNumber(reader, header, endRow):
    responseColumn = header.index('Response')
    responseCount = 0
    for i, line in enumerate(reader):
        if i == endRow:
            break
        if line[responseColumn] == '+':
            responseCount += 1
    return responseCount - 1


def getStartAndFinishRowsOfFirstID(reader, header): #hint the return type
    displayColumn = header.index('display')
    trialNumberColumn = header.index('Trial Number')

    startRow = None
    endRow = None

    for i, line in enumerate(reader):
        if line[displayColumn] == 'trial' and startRow is None:
            startRow = i
            print(f'first time we see trial: {line[displayColumn]} at i = {i}')
        if line[trialNumberColumn] == 'END TASK' and endRow is None:
            endRow = i - 1
            break
    return startRow, endRow

def getStimuli(reader, header, startRow, endRow):
    leftStimColumn = header.index('left_stim')
    conditions = getConditions(reader, header, endRow)
    stimuli = {}

    #if number of conditions == 1:
    # if len(conditions) == 1:
    #     onlyCondition = conditions[0]
    #     for i, line in enumerate(reader):
    #         if i >= startRow and i <= endRow:
    #             indexUpToExtension = line[leftStimColumn].find('.') - 1 
    #             value = line[leftStimColumn][:indexUpToExtension + 1]
    #             if onlyCondition not in stimuli:
    #                 stimuli[onlyCondition] = [value]
    #             else:
    #                 stimuli[onlyCondition].append(value)

    #if number of conditions > 1:
    # elif len(conditions) > 1:
    for i, line in enumerate(reader):
        # print(f'i: {i}, {line[leftStimColumn]}')
        if i >= startRow and i <= endRow:
            indexUpToExtension = line[leftStimColumn].find('.') - 1            # -1 assuming N is < 10
            possibleKey = line[leftStimColumn][:indexUpToExtension]
            if possibleKey == '': continue
            value = line[leftStimColumn][:indexUpToExtension + 1]
            if possibleKey not in stimuli:
                stimuli[possibleKey] = [value]
            else:
                if value not in stimuli[possibleKey]:
                    stimuli[possibleKey].append(value)
        elif i > endRow:
            break
    return stimuli
    
def getConditions(reader, header, endRow):
    conditionsColumn = header.index('condition')
    conditions = []
    for i, line in enumerate(reader):
        if i == endRow: break
        condition = line[conditionsColumn]
        if condition not in conditions:
            conditions.append(condition)
    return conditions

def printValAtStartRow(reader, header, startRow):
    print(f'startRow in printValAtStartRow: {startRow}')
    col = header.index('display')
    for i, line in enumerate(reader):
        print(f'i: {i}, {line[col]}')
        if i == startRow:
            print(f'val at startRow: {line[col]}')
            break
    
getVariablesFromCSV(file_path)



