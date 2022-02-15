import xlrd 

#access a cell 
# folder = input('where to look?\n')
EXPERIMENT_FOLDER = input('What is the name of the folder where your data is? \
    (press enter when finished)\n')
    
wb = xlrd.open_workbook(EXPERIMENT_FOLDER + '/userInputVariables.xls')
sheet = wb.sheet_by_index(0)

def getTrialDuration():
    return sheet.cell_value(0,1)

def getStimWidth():
    return sheet.cell_value(1,1)

def getStimHeight():
    return sheet.cell_value(2,1)

def getCSVFileName():
    return sheet.cell_value(3,1)

def getTrialDataFolder():
    return sheet.cell_value(4,1)

