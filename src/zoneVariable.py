import variableGetters 
import gorilla

TRIAL_DURATION = variableGetters.getTrialDuration()
STIM_width = variableGetters.getStimWidth()
STIM_height = variableGetters.getStimHeight()
DATA_FILE = variableGetters.EXPERIMENT_FOLDER + '/' + variableGetters.getCSVFileName()
TRIAL_DATA_FOLDER = variableGetters.EXPERIMENT_FOLDER + '/' + variableGetters.getTrialDataFolder()

STIM_ratio = STIM_height/STIM_width

def createPicInfoDict(dimensions):
    picInfo = {}
    x, y, w, h = dimensions
    picInfo['startX'] = x
    picInfo['endX'] = x + w
    picInfo['startY'] = y
    picInfo['endY'] = y + h
    return picInfo
    
def isInPicture(x,y, picInfo):
    isInXbound = x >= picInfo['startX'] and x <= picInfo['endX']
    isInYbound = y >= picInfo['startY'] and y <= picInfo['endY']
    return isInXbound and isInYbound

def checkValidTrial(totalOnLeft, totalOnRight): 
    return (totalOnLeft + totalOnRight) >= TRIAL_DURATION/4

def calculateDwellTimes(xCoordinates, yCoordinates, times, leftStimPicInfo, rightStimPicInfo):
    numSamples = len(xCoordinates)
    totalOnLeft = 0
    totalOnRight = 0

    for i in range(1, numSamples):
        timeDifference = times[i] - times[i-1]
        if isInPicture(xCoordinates[i], yCoordinates[i], leftStimPicInfo):
            totalOnLeft += timeDifference
        elif isInPicture(xCoordinates[i], yCoordinates[i], rightStimPicInfo):
            totalOnRight += timeDifference
    
    return totalOnLeft, totalOnRight

def createAOIRect(participant, data):
    '''
    go get code from analysis script
    reach into first trial and create aoi rect
    width and height will be hardcoded...
    return aoi rect dictionary
    '''
    # Reformat zone shape messages to match old format
    # Create dict with zone coordinates
    Zone_rect = {}
    coord_list = ['x','y','w','h']
    if len(data[participant]["trials"][0]["msg"]) >3: # >3 so it won't crash on empty trials without this information in msgs
        for zone in ['Zone1','Zone2']:
            Zone_rect[zone] = []
            if zone == 'Zone1':
                Zonemsg = data[participant]["trials"][0]["msg"][6][1]
            else:
                Zonemsg = data[participant]["trials"][0]["msg"][7][1]
            Zonemsg_split = Zonemsg.split(';')
            for ji,coord in enumerate(coord_list):
                if coord == 'x':
                        temp = Zonemsg_split[2]
                elif coord == 'y':
                        temp = Zonemsg_split[3]
                elif coord == 'w':
                        temp = Zonemsg_split[4]
                elif coord == 'h':
                        temp = Zonemsg_split[5]
                
                temp_split = temp.split('=')
                Zone_rect[zone].append(int(round(float(temp_split[1]))))

    ''' at this point zone_rect = {'Zone1': [x,y,w,h], 'Zone2': [x,y,w,h]} '''
    
    aoi_rect = {}
    data[participant]["trials"][0]["aoi_coord"] = []
    if len(data[participant]["trials"][0]["msg"]) >3:
        for zone in ['Zone1','Zone2']:
            i_rect = []
            z_rect = Zone_rect[zone]
            zone_center = [(z_rect[0]+(z_rect[2]/2)),(z_rect[1]+(z_rect[3]/2))] #list with x,y for zone center
            if z_rect[2] >= STIM_width :#presented in actual size
                i_rect = [zone_center[0]-(STIM_width/2),(zone_center[1]-(STIM_height/2)),STIM_width,STIM_height] #x,y,w,h for actual size
            else: #shrunk to width of zone w aspect ratio maintained
                i_rect = [zone_center[0]-(z_rect[2]/2),zone_center[1]-((z_rect[2]*STIM_ratio)/2),z_rect[2],z_rect[2]*STIM_ratio]
            for ir in range(len(i_rect)):
                i_rect[ir] = int(round(i_rect[ir]))
            if zone == "Zone1":
                z_field = 'left_stim'
            else:
                z_field = 'right_stim'
            aoi_rect[z_field] = i_rect
            # data[participant]["trials"][0][zone] = i_rect             removed because unnecessary
    return aoi_rect

def createDwellTimesDict():
    data = gorilla.read_file(DATA_FILE, TRIAL_DATA_FOLDER, \
    custom_fields=["condition", "left_stim", "right_stim"], \
    verbose=True)

    pictureDwellTimes = {}

    for participant in data:
        aoi_rect = createAOIRect(participant, data)
        leftStimPicInfo = createPicInfoDict(aoi_rect['left_stim'])
        rightStimPicInfo = createPicInfoDict(aoi_rect['right_stim'])
        for trial in data[participant]['trials']:
            x = trial['x']
            y = trial['y']
            time = trial['time']
            totalOnLeft, totalOnRight = calculateDwellTimes(x, y, time, \
                leftStimPicInfo, rightStimPicInfo)
            isValidTrial = checkValidTrial(totalOnLeft, totalOnRight)

            pictureNameRight = trial['msg'][0][1]
            pictureNameLeft = trial['msg'][1][1]

            if pictureNameRight not in pictureDwellTimes:
                pictureDwellTimes[pictureNameRight] = {}
            if participant not in pictureDwellTimes[pictureNameRight]:
                pictureDwellTimes[pictureNameRight][participant] = \
                    {'total': 0, 'average': 0, 'raw': [], 'numDropped': 0}

            if pictureNameLeft not in pictureDwellTimes:
                pictureDwellTimes[pictureNameLeft] = {}
            if participant not in pictureDwellTimes[pictureNameLeft]:
                pictureDwellTimes[pictureNameLeft][participant] = \
                    {'total': 0, 'average': 0, 'raw': [], 'numDropped': 0}

            dwellPicLeft = pictureDwellTimes[pictureNameLeft][participant]
            dwellPicRight = pictureDwellTimes[pictureNameRight][participant]

            if not isValidTrial:
                dwellPicLeft['numDropped'] += 1
                dwellPicLeft['raw'].append(('INVALID', pictureNameRight))
                dwellPicRight['numDropped'] += 1
                dwellPicRight['raw'].append(('INVALID', pictureNameLeft))
                continue

            dwellPicLeft['total'] += totalOnLeft
            dwellPicLeft['raw'].append((totalOnLeft, pictureNameRight))
            dwellPicLeft['average'] = dwellPicLeft['total'] / len(dwellPicLeft['raw']) - dwellPicLeft['numDropped']
            pictureDwellTimes[pictureNameLeft][participant] = dwellPicLeft

            dwellPicRight['total'] += totalOnRight
            dwellPicRight['raw'].append((totalOnRight, pictureNameLeft))
            dwellPicRight['average'] = dwellPicRight['total'] / len(dwellPicRight['raw']) - dwellPicRight['numDropped']
            pictureDwellTimes[pictureNameRight][participant] = dwellPicRight
            
    return pictureDwellTimes
