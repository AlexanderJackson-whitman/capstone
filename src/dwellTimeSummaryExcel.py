from zoneVariable import *
from variableGetters import EXPERIMENT_FOLDER
from xlwt import Workbook

def writeHeader(pictureDwellTimes, sheet, totalFrequencyShown):
    currColumn = 1
    for currentPicShownFreq in range(1, totalFrequencyShown + 1):
        for pictureName in pictureDwellTimes:
            if 'neutral' in pictureName:
                continue
            headerName = pictureName + '_' + str(currentPicShownFreq)
            sheet.write(0, currColumn, headerName)
            currColumn += 1

def getTotalFrequencyShown(pictureDwellTimes, randomPictureName, randomParticipant):
    return len(pictureDwellTimes[randomPictureName][randomParticipant]['raw'])

def writeParticipantIDs(pictureDwellTimes, sheet, randomPictureName):
    sheet.write(0, 0, 'Participant ID')
    currRow = 1
    for participant in pictureDwellTimes[randomPictureName]:
        sheet.write(currRow, 0, participant)
        currRow += 1

def writeAllColumns(pictureDwellTimes, sheet, totalFrequencyShown):
    def writeSingleColumn(pictureName, column, freq):
        row = 1
        for participant in pictureDwellTimes[pictureName]:
            cellData = pictureDwellTimes[pictureName][participant]['raw'][freq][0]
            sheet.write(row, column, cellData)
            row += 1

    currColumn = 1

    for currentPicShownFreq in range(totalFrequencyShown):
        for pictureName in pictureDwellTimes:
            if 'neutral' in pictureName:
                continue
            writeSingleColumn(pictureName, currColumn, currentPicShownFreq)
            currColumn += 1

def getRandomPictureName(pictureDwellTimes):
    for pictureName in pictureDwellTimes:
        if 'neutral' in pictureName:
            continue
        return pictureName

def main():
    wb = Workbook()
    pictureDwellTimes = createDwellTimesDict()
    sheet = wb.add_sheet('Sheet 1')

    randomPictureName = getRandomPictureName(pictureDwellTimes)
    randomParticipant = list(pictureDwellTimes[randomPictureName].keys())[0]
    totalFrequencyShown = getTotalFrequencyShown(pictureDwellTimes, randomPictureName, randomParticipant)

    writeHeader(pictureDwellTimes, sheet, totalFrequencyShown)
    writeParticipantIDs(pictureDwellTimes, sheet, randomPictureName)
    writeAllColumns(pictureDwellTimes, sheet, totalFrequencyShown)

    wb.save(EXPERIMENT_FOLDER + '/pictureDwellSummaryTimes.xls')

if __name__ == '__main__':
    main()