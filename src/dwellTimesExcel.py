import xlwt
from zoneVariable import *
from xlwt import Workbook
  
# Workbook is created
wb = Workbook()
  
# add_sheet is used to create sheet.
pictureDwellTimes = createDwellTimesDict() #FIX THIS LATER 1/26

for picture in pictureDwellTimes:
    sheet = wb.add_sheet(picture)
    sheet.write(0, 0, 'Participant ID')
    sheet.write(0, 1, 'Total Dwell Time (ms)')
    sheet.write(0, 2, 'Average Dwell Time (ms)')
    sheet.write(0, 3, 'Raw Dwell Data (ms)')
    sheet.write(0, 4, 'Paired Stimulus')

    row = 1
    for participant in pictureDwellTimes[picture]:
        sheet.write(row, 0, participant)
        sheet.write(row, 1, pictureDwellTimes[picture][participant]['total'] )
        sheet.write(row, 2, pictureDwellTimes[picture][participant]['average'])
        for dwellTime, comparedPic in pictureDwellTimes[picture][participant]['raw']:
            sheet.write(row, 3, dwellTime)
            sheet.write(row, 4, comparedPic)
            row += 1
  
wb.save('pictureDwellTimes.xls')


