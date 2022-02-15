from zoneVariable import *
from variableGetters import EXPERIMENT_FOLDER
from xlwt import Workbook

def writeHeader(pictureDwellTimes, sheet):
    # TODO : write getTotalFreq()

    totalPicShownFreq = 2 #getTotalFreq()
    currColumn = 1
    for currentPicShownFreq in range(1, totalPicShownFreq + 1):
        for pictureName in pictureDwellTimes:
            if 'neutral' in pictureName:
                continue
            headerName = pictureName + '_' + str(currentPicShownFreq)
            sheet.write(0, currColumn, headerName)
            currColumn += 1

def writeParticipantIDs(pictureDwellTimes, sheet):
    sheet.write(0, 0, 'Participant ID')
    randomPictureName = list(pictureDwellTimes.keys())[0]
    currRow = 1
    for participant in pictureDwellTimes[randomPictureName]:
        sheet.write(currRow, 0, participant)
        currRow += 1

def writeAllColumns(pictureDwellTimes, sheet):
    def writeSingleColumn(pictureName, column, freq):
        row = 1
        for participant in pictureDwellTimes[pictureName]:
            cellData = pictureDwellTimes[pictureName][participant]['raw'][freq][0]
            sheet.write(row, column, cellData)
            row += 1

    totalPicShownFreq = 2 #getTotalFreq()
    currColumn = 1

    for currentPicShownFreq in range(totalPicShownFreq):
        for pictureName in pictureDwellTimes:
            if 'neutral' in pictureName:
                continue
            writeSingleColumn(pictureName, currColumn, currentPicShownFreq)
            currColumn += 1


def main():
    wb = Workbook()
    pictureDwellTimes = createDwellTimesDict() #FIX THIS LATER 1/26
    sheet = wb.add_sheet('Sheet 1')
    writeHeader(pictureDwellTimes, sheet)
    writeParticipantIDs(pictureDwellTimes, sheet)
    writeAllColumns(pictureDwellTimes, sheet)

    wb.save(EXPERIMENT_FOLDER + '/pictureDwellSummaryTimes.xls')

if __name__ == '__main__':
    main()

'''
PSEUDOCODE FOR ROW FUNCTION:


column = 2
totalPicShownFreq = len(raw)


for currentPicShownFreq in range(1, totalPicShownFreq + 1):
    for every picture:
        newName = pictureName + _ + currPicShownFreq
        put in cell
        column += 1

'''