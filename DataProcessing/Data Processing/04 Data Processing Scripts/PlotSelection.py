 #!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import operator
import scipy.stats as stats
import os
import csv

from statsmodels.stats.anova import AnovaRM
#import statsmodels.api as sm
#from statsmodels.formula.api import ols

from ImportFiles import ImportFiles
from FigureInfo import FigureInfo
from PlotData import PlotData

# Nicer logging:
import logging
logging.basicConfig()
log = logging.getLogger()
# Make the logs appear prettyer
formatter = logging.Formatter('%(levelname)s:%(message)s')
for handler in log.handlers:
    handler.setFormatter(formatter)

class PlotSelection(object):
    def __init__(self, basePathData, doSwarm = True, plotqq = False, doTransformation = False, doCheckInteractions = False, loglevel = logging.INFO):
        # Class settings
        log.setLevel(loglevel)
        self.outputPath = "../03 Statistical Data/"  # This is the path into which the figures will be saved
        self.analysisPath = "../02 Analysis Data/"   # This is the path into which the analysis data will be saved

        self.shapirodatapath = "../03 Statistical Data/shapiroresults.csv"  #  This file contains all results of the shapiro wilks test for easy referencing.
        self.anovadatapath = "../03 Statistical Data/anovaresults.csv"  #  This file contains all results of the shapiro wilks test for easy referencing.

        self.plotqq = plotqq
        self.doTransformation = doTransformation
        self.doCheckInteractions = doCheckInteractions

        # init classes
        self.IF = ImportFiles(basePathData)
        self.PD = PlotData(self.outputPath, doSwarm)
        self.FI = FigureInfo()

        self.labels = ["MC","HC","MN","HN"]  #  M=Monitor,H=HMD;C=Cues,N=No cues

    ########: Time plots ########:
    def timePlot(self, cond, logType, data, names = None, trialNo = [0,1]):  #  These can all be done using a single function, as the data is all of the same structure
        if isinstance(data,str):
            data = [data]
        if isinstance(trialNo,int):
            trialNo = [trialNo]

        if names == None or names == "all":
            names = self.IF.allParticipants

        fig = None # Figure Handle
        axes = [None for i in range(len(data))]

        for name in names:
            self.IF.selectParticipant(name)
            dfDict = self.IF.findLog(cond, logType, trialNo) # contains all logs of a specified type
            infoList = [self.FI.getInfo(logType, currData) for currData in data]

            for i,file in enumerate(dfDict): # each df will make a new figure
                fileSplit = file.split("_")
                repnum = fileSplit[-1].replace(".csv","")

                supTitle = fileSplit[1] + " " + fileSplit[0]

                if len(data) > 1: # each info entry will make a new subplot
                    fig,axes[0] = self.PD.timeFigure(dfDict[file],infoList[0],fig,axes[0],supTitle)
                    for l,info in enumerate(infoList[1:]): # skip the first entry
                        fig,axes[l+1] = self.PD.addTimeSubPlot(dfDict[file],info,fig,axes[l+1])
                else:
                    fig,axes[0] = self.PD.timeFigure(dfDict[file],infoList[0],fig,axes[0],supTitle)

    ########: Calculate Data and Plot Results ########:

    ### Continuous Metrics

    # Removed dross fraction
    def plotDross(self, includefails = False):
        #for each participant
        dataList = []

        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)

            #for each condition
            dataList_p = []
            for i in range(4):
                resultList = []
                logType = "obi"
                dfDict = self.IF.findLogType(i, logType) # contains all logs of a specified trial number

                drossInBin = 0
                totalParticles = 0
                for file in dfDict:
                    bin1infoList = self.FI.getInfo(logType, "bin1")
                    bin2infoList = self.FI.getInfo(logType, "bin2")
                    bathInfoList = self.FI.getInfo(logType, "bath")

                    timeList = list(dfDict[file][bin1infoList["xAxis"]])

                    if self.checkTime(timeList) or includefails:
                        bin1LastEntry = list(dfDict[file][bin1infoList["yAxis"]])[-1] # last entry
                        bin2LastEntry = list(dfDict[file][bin2infoList["yAxis"]])[-1] # last entry
                        drossInBin += bin1LastEntry + bin2LastEntry
                        totalParticles += 750  #  max(list(dfDict[file][bathInfoList["yAxis"]])) would also work

                fractionRemoved = drossInBin/totalParticles
                dataList_p.append(fractionRemoved)

            dataList.append(dataList_p)

        # Plotting & Saving Figure
        self.boxplot(dataList,'Dross Removed','Fraction removed')

        return dataList

    def plotAccuracy(self, includefails = False):
        #for each participant
        dataList = []

        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)

            #for each condition
            dataList_p = []
            for i in range(4):
                resultList = []
                logType = "obi"
                dfDict = self.IF.findLogType(i, logType) # contains all logs of a specified trial number

                drossInBath = []
                drossInBin = []
                totalParticles = 0
                for file in dfDict:
                    bin1infoList = self.FI.getInfo(logType, "bin1")
                    bin2infoList = self.FI.getInfo(logType, "bin2")
                    bathInfoList = self.FI.getInfo(logType, "bath")

                    timeList = list(dfDict[file][bin1infoList["xAxis"]])
                    if self.checkTime(timeList) or includefails:
                        bin1LastEntry = list(dfDict[file][bin1infoList["yAxis"]])[-1] # last entry
                        bin2LastEntry = list(dfDict[file][bin2infoList["yAxis"]])[-1] # last entry
                        drossInBin.append(bin1LastEntry + bin2LastEntry)
                        drossInBath.append(list(dfDict[file][bathInfoList["yAxis"]])[-1])
                        #drossMissed = 750-drossInBath-drossInBin
                        totalParticles += 750

                accuracy = sum(drossInBin)/(totalParticles-sum(drossInBath))
                dataList_p.append(accuracy)
                if accuracy < .75: log.debug("Participant {} had accuracy of ({}) in condition {}".format(name,accuracy,i))

            dataList.append(dataList_p)

        # Plotting & Saving Figure
        self.boxplot(dataList,'Accuracy dross dumping','Fraction successful')

        return dataList

    def plotAvgScoopSize(self, includefails = False):
        #for each participant
        dataList = []

        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)

            #for each condition
            dataList_p = []
            for i in range(4):
                resultList = []
                logType = "obi"
                dfDict = self.IF.findLogType(i, logType) # contains all logs of a specified trial number

                scoopSizes = []
                dumpSizes = []
                for file in dfDict:
                    bin1infoList = self.FI.getInfo(logType, "bin1")
                    bin2infoList = self.FI.getInfo(logType, "bin2")
                    bathInfoList = self.FI.getInfo(logType, "bath")
                    timeList = list(dfDict[file][bin1infoList["xAxis"]])
                    startElement = self.findStartTime(timeList)+1  #  used to cut off the leading part of the data (for some reason part of the previous dataset flows over to the next one)

                    if self.checkTime(timeList) or includefails:
                        bin1 = list(dfDict[file][bin1infoList["yAxis"]])[startElement:]
                        bin2 = list(dfDict[file][bin2infoList["yAxis"]])[startElement:]
                        bath = list(dfDict[file][bathInfoList["yAxis"]])[startElement:]

                        d_bin1 = np.abs(np.diff(bin1))  #  we dont care whether it is up (bin) or down (bath) for the purpose of this function
                        d_bin2 = np.abs(np.diff(bin2))
                        d_bath = np.abs(np.diff(bath))
                        # find the starts of each scoop motion
                        scoopStarts = []
                        dumpStarts = []
                        j = 0
                        l = 0
                        while j < len(bin1)-2:
                            if d_bath[j]>0: #  dross is taken from the bath
                                scoopStarts.append(j)
                                while d_bin1[j+l] == 0 and d_bin2[j+l] ==0 and j+l<len(bin1)-2:  #  dross is dumped in the bins
                                    l += 1
                                dumpStarts.append(j+l)
                                j += l+1
                                l = 0
                            else:
                                j += 1

                        for k in range(len(scoopStarts)):
                            scoopSizes.append(bath[scoopStarts[k]] - bath[dumpStarts[k]])
                            if k + 1 < len(scoopStarts):
                                dumpSizes.append((bin1[scoopStarts[k+1]] - bin1[dumpStarts[k]]) + (bin2[scoopStarts[k+1]] - bin2[dumpStarts[k]]))
                NScoops = len(scoopSizes)
                NDumps = len(dumpSizes)
                loss = [scoopSizes[x] - dumpSizes[x] for x in range(NDumps)]
                #print(scoopSizes)
                avgScoop = sum(scoopSizes)/NScoops
                avgDump = sum(dumpSizes)/NDumps
                avgLoss = sum(loss)/NDumps   # very minimal, not really relevant

                dataList_p.append(avgDump)
                #if dataList_p[-1] > 25: log.debug("Participant {} score of ({}) in condition {}".format(name,dataList_p[-1] ,i))

            dataList.append(dataList_p)

        # Plotting & Saving Figure
        self.boxplot(dataList,'Average scoopsize','Dross particles (count)')

        return dataList

    def plotNScoops(self, includefails = False):
        #for each participant
        dataList = []

        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)

            #for each condition
            dataList_p = []
            for i in range(4):
                resultList = []
                logType = "obi"
                dfDict = self.IF.findLogType(i, logType) # contains all logs of a specified trial number

                NScoops = 0
                for file in dfDict:
                    bin1infoList = self.FI.getInfo(logType, "bin1")
                    bin2infoList = self.FI.getInfo(logType, "bin2")
                    bathInfoList = self.FI.getInfo(logType, "bath")
                    timeList = list(dfDict[file][bin1infoList["xAxis"]])
                    startElement = self.findStartTime(timeList)+1  #  used to cut off the leading part of the data (for some reason part of the previous dataset flows over to the next one)

                    if self.checkTime(timeList) or includefails:
                        bin1 = list(dfDict[file][bin1infoList["yAxis"]])[startElement:]
                        bin2 = list(dfDict[file][bin2infoList["yAxis"]])[startElement:]
                        bath = list(dfDict[file][bathInfoList["yAxis"]])[startElement:]

                        d_bin1 = np.abs(np.diff(bin1))  #  we dont care whether it is up (bin) or down (bath) for the purpose of this function
                        d_bin2 = np.abs(np.diff(bin2))
                        d_bath = np.abs(np.diff(bath))
                        # find the starts of each scoop motion
                        scoopStarts = []
                        dumpStarts = []
                        j = 0
                        l = 0
                        while j < len(bin1)-2:
                            if d_bath[j]>0: #  dross is taken from the bath
                                scoopStarts.append(j)
                                while d_bin1[j+l] == 0 and d_bin2[j+l] ==0 and j+l<len(bin1)-2:  #  dross is dumped in the bins
                                    l += 1
                                dumpStarts.append(j+l)
                                j += l+1
                                l = 0
                            else:
                                j += 1
                        NScoops += len(scoopStarts)

                dataList_p.append(NScoops)
                #if dataList_p[-1] > 25: log.debug("Participant {} score of ({}) in condition {}".format(name,dataList_p[-1] ,i))

            dataList.append(dataList_p)

        # Plotting & Saving Figure
        self.boxplot(dataList,'Number of scoops','Number of scoops (count)')

        return dataList

    def plotCollision(self, includefails = False):
        #for each participant
        dataList = []
        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)

            #for each condition
            dataList_p = []
            for i in range(4):
                resultList = []
                logType = "RB"

                dfDict = self.IF.findLogType(i, logType) # contains all logs of a specified trial number

                result = 0
                for file in dfDict:
                    infoList0 = self.FI.getInfo(logType, "mag")

                    timeList = list(dfDict[file][infoList0["xAxis"]])
                    startElement = self.findStartTime(timeList)+1  #  used to cut off the leading part of the data (for some reason part of the previous dataset flows over to the next one)

                    if self.checkTime(timeList) or includefails:

                        forceVec = list(dfDict[file][infoList0["yAxis"]])[startElement:]
                        max_index, maxForce = max(enumerate(forceVec), key=operator.itemgetter(1))
                        trimmedtime = timeList[startElement:]

                        # Theres some weird business going on with forces locking at high values when disconnecting the grabjoint in unity
                        # The force shoots up and stays the exact same value, which doesn't happen otherwise, this is whats used to filter these out
                        # Furthermore, sometimes the logging rate exceeds the rb update rate, resulting in two identical datapoints, thats why 3 datapoints are taken instead to find these frozen forces
                        # The results of this "algorithm" were checked with the real data to make sure everything works as it should
                        if max_index + 2 < len(forceVec) - 1: #  When the max force is the final force, else this throws an error
                            while forceVec[max_index] == forceVec[max_index+1] and forceVec[max_index] == forceVec[max_index+2] or forceVec[max_index] == forceVec[max_index-1] and forceVec[max_index] == forceVec[max_index-2]:
                                log.debug("Incorrect force of mag: " + str(maxForce) + " at Time stamp: " + str(trimmedtime[max_index]) + ", rescanning for smaller force")
                                #trimmedtime =  [trimmedtime for e in forceVec if e < maxForce]
                                forceVec = [e for e in forceVec if e < maxForce]  #  Only take values < the incorrect maxvalue
                                max_index, maxForce = max(enumerate(forceVec), key=operator.itemgetter(1))
                                if max_index + 2 > len(forceVec)-1:
                                    break
                        log.debug("Found correct maxforce: "+ str(maxForce) + "N, Time stamp: " + str(trimmedtime[max_index]))  #  The timestamps wont be 100% correct here as some elements have been removed, but they will be close (so you can find the points in the dataset)
                        if maxForce > result:  # take the maximum of the two trials
                            result = maxForce

                if result > 40: result = 40  #  This was the threshold set in the simulation at which the trial would fail (due to stepsize this would sometimes surpass 40, though higher then 40 doesn't make sense, so 40 is the maximum)
                dataList_p.append(result)

            dataList.append(dataList_p)

        # Plotting & Saving Figure
        self.boxplot(dataList,'Maximum collision force','Force [N]')

        return dataList

    def plotProximity(self, includefails = False):
        #for each participant
        dataList = []
        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)

            #for each condition
            dataList_p = []
            for i in range(4):
                resultList = []
                logType = "prox"
                logTypeIF = "proximity_Box002"  #  there are also logs of the wrist of the robot, but these are not so interesting as that never got close to the steel strip.

                dfDict = self.IF.findLogType(i, logTypeIF) # contains all logs of a specified trial number

                result = 0
                for file in dfDict:
                    infoList0 = self.FI.getInfo(logType, "mag")

                    timeList = list(dfDict[file][infoList0["xAxis"]])
                    startElement = self.findStartTime(timeList)+1  #  used to cut off the leading part of the data (for some reason part of the previous dataset flows over to the next one)

                    if self.checkTime(timeList) or includefails:
                        proxVec = list(dfDict[file][infoList0["yAxis"]])[startElement:]
                        min_index, minProx = min(enumerate(proxVec), key=operator.itemgetter(1))
                        trimmedtime = timeList[startElement:]

                        #log.debug("{} - had a minimum distance to steel strip: {} m, Time stamp: {}".format(name,minProx,trimmedtime[min_index]))  #  The timestamps wont be 100% correct here as some elements have been removed, but they will be close (so you can find the points in the dataset)
                        if minProx > result:  # take the maximum of the two trials
                            result = minProx

                dataList_p.append(result)

            dataList.append(dataList_p)

        # Plotting & Saving Figure
        self.boxplot(dataList,'Minimum distance to steel strip','Distance [m]')

        return dataList

    def plotScoopFlux(self, includefails = False):
        #for each participant
        dataList = []
        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)

            #for each condition
            dataList_p = []
            for i in range(4):
                resultList = []
                logType = "ivolume"

                dfDict = self.IF.findLogType(i, logType) # contains all logs of a specified trial number

                result = 0
                for file in dfDict:
                    infoList0 = self.FI.getInfo(logType, "flux")

                    timeList = list(dfDict[file][infoList0["xAxis"]])
                    startElement = self.findStartTime(timeList)+1  #  used to cut off the leading part of the data (for some reason part of the previous dataset flows over to the next one)
                    timeList = timeList[startElement:]
                    #log.debug("Timelist[0]: {}".format(timeList[0]))

                    if self.checkTime(timeList) or includefails:
                        fluxVec = list(dfDict[file][infoList0["yAxis"]])[startElement:]
                        max_index, maxFlux = max(enumerate(fluxVec), key=operator.itemgetter(1))
                        maxFlux = max(fluxVec)

                        if maxFlux > result:  # take the maximum of the two trials
                            result = maxFlux
                            resultTime = timeList[max_index]

                #log.debug(("Found correct maxflux: {}, Time stamp: {}, person: {}").format(result,resultTime,name))
                dataList_p.append(result)

            dataList.append(dataList_p)

        # Plotting & Saving Figure
        self.boxplot(dataList,'Maximum scoop flux','Flux')

        return dataList

    def plotAvgVel(self, includefails = False):
        dataList = []

        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)

            dataList_p = []
            for i in range(4):
                resultList = []
                dfDict = self.IF.findLogType(i, "RB")

                finalTs = []
                V = 0
                for file in dfDict:  # for both of the reps
                    logNum = file.strip(".csv").split("_")[-1]
                    file = self.findLogNum(logNum,dfDict)
                    infoList = self.FI.getInfo("RB", "vel")
                    t = list(dfDict[file][infoList["xAxis"]])

                    if self.checkTime(t) or includefails:  #  check if it failed, else we skip it
                        s = self.findStartTime(t)+1  #  Start element (theres some excess data at the start due to the way the loggers were written in unity)
                        t = t[s:]  #  corrected time vector: Trim of the start data
                        dt = np.diff(t)
                        #V = 0
                        # put the three "axial" velocities in a single vector
                        speeds = []
                        for m in range(3): speeds.append(list(dfDict[file][infoList["yAxis"][m]])[s:])

                        for l in range(len(dt)):  # Integrate over time
                            speedVec = [speeds[0][l],speeds[1][l],speeds[2][l]]
                            speed = self.mag(speedVec)
                            V += speed * dt[l]

                        finalTs.append(t[-1])
                        #averageVel = V/t[-1] #  Divide over total time
                        #resultList.append(averageVel)  # put into the resulting array
                # Take find the average velocity
                averageVel = V/sum(finalTs) #  Divide over total time
                dataList_p.append(averageVel)
                #dataList_p.append(self.avgData(resultList))

            dataList.append(dataList_p)

        # Plotting & Saving Figure
        self.boxplot(dataList,'Average velocity','Average velocity [m/s]')

        return dataList

    def plotDisturbance(self, includefails = False):
        scaleByPerformance = False
        dataList = []

        if scaleByPerformance:
            drossData = self.plotDross(includefails)

        #drossResults = self.plotDross()
        for pnum, name in enumerate(self.IF.allParticipants):
            self.IF.selectParticipant(name)

            dataList_p = []
            for i in range(4):
                dfDict_iv = self.IF.findLogType(i, "ivolume")
                dfDict_rb = self.IF.findLogType(i, "RB")

                dDisturbance = 0
                finalTs = []
                for file_iv in dfDict_iv:
                    logNum = file_iv.strip(".csv").split("_")[-1]
                    file_rb = self.findLogNum(logNum,dfDict_rb)
                    infoList_iv = self.FI.getInfo("ivolume", "vol")
                    infoList_rb = self.FI.getInfo("RB", "vel")

                    # Still need to modify the time vectors so they correctly line up, this part of the code achieves just that
                    t_iv = list(dfDict_iv[file_iv][infoList_iv["xAxis"]])
                    t_rb = list(dfDict_rb[file_rb][infoList_rb["xAxis"]])
                    if self.checkTime(t_iv) or includefails:  #  check if it failed, else we skip it
                        s_iv = self.findStartTime(t_iv)+1 # element num to start at
                        s_rb = self.findStartTime(t_rb)+1
                        maxLength = min([len(t_iv)-s_iv,len(t_rb)-s_rb])  # max number of elements that can be used

                        t = t_iv[s_iv:s_iv+maxLength]  #  same as t_rb[s_rb:s_rb+maxLength]

                        dt = np.diff(t)
                        speeds = []
                        for m in range(3):
                            speeds.append(list(dfDict_rb[file_rb][infoList_rb["yAxis"][m]])[s_rb:s_rb+maxLength])
                        iv = list(dfDict_iv[file_iv][infoList_iv["yAxis"]])[s_iv:s_iv+maxLength]

                        for l in range(len(dt)):
                            speedVec = [speeds[0][l],speeds[1][l],speeds[2][l]]
                            speed = self.mag(speedVec)
                            dDisturbance += dt[l] * speed * iv[l]

                        finalTs.append(t[-1])

                averageD = dDisturbance/sum(finalTs) #  Divide over total time

                if scaleByPerformance:
                    averageD /= drossData[pnum][i]

                if averageD>.05:
                    log.debug(("{} had disturbance: {}, conditionNo: {}").format(name,averageD,i))  #  The timestamps wont be 100% correct here as some elements have been removed, but they will be close (so you can find the points in the dataset)

                dataList_p.append(averageD)
                #dataList_p = [dataList_p[l]/drossResults[i][l] for l in range(len(dataList_p))]

            dataList.append(dataList_p)

        # Plotting & Saving Figure
        self.boxplot(dataList,'Average Disturbance','Average Disturbance')
        return dataList

    def plotvdLaan(self):  #  https://www.hfes-europe.org/accept/accept.htm
        #for each participant, Not actually continuous, but data analysis works the same here
        dataLists_s = []
        dataLists_u = []

        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)

            dataList_p_s = []
            dataList_p_u = []

            for i in range(4): # for each condition
                dfDict = self.IF.findvdLaan(i)

                corrValues = dfDict["Value"]  #  range recorded 0..4, need +2..-2
                corrValues = (corrValues-2)*-1

                mirrorCorr = [2,5,7]   #  the questions that need their score mirrored
                corrValues[mirrorCorr] *= -1

                usefulness_qnums = [0,2,4,6,8]
                usefulness = sum(corrValues[usefulness_qnums])/len(usefulness_qnums)

                satisfying_qnums = [1,3,5,7]
                satisfying = sum(corrValues[satisfying_qnums])/len(satisfying_qnums)

                dataList_p_s.append(satisfying)
                dataList_p_u.append(usefulness)

            dataLists_s.append(dataList_p_s)
            dataLists_u.append(dataList_p_u)

        # Plotting
        self.vdLaanPlot(dataLists_u, dataLists_s,'van der Laan Acceptance Scale')

        self.boxplot(dataLists_u,'van der laan Usefulness','Usefulness score')
        self.boxplot(dataLists_s,'van der laan Satisfying','Satisfying score')

        return dataLists_u, dataLists_s


    ### Discrete Metrics
    def plotFailure(self):
        #for each participant
        dataList = []
        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)

            dataList_p = []
            for i in range(4):
                resultList = []

                dfDict = self.IF.findLogType(i, "RB")

                for l,file in enumerate(dfDict):
                    logType = file.split("_")[0]
                    infoList0 = self.FI.getInfo(logType, "mag")

                    timeList = list(dfDict[file][infoList0["xAxis"]])
                    forceList = list(dfDict[file][infoList0["yAxis"]])
                    if self.checkTime(timeList):
                        result = 1
                        resultList.append(result)
                    else:
                        log.info("Participant {} failed condition {}, trialno: {}, at timestamp {}, final collision force: {}".format(name,i,l,timeList[-1],forceList[-1]))
                        if forceList[-1] < 40.0:
                            self.checkCollisionsForSteelStrip(name)

                if (len(resultList)>0):
                    dataList_p.append(2-sum(resultList))
                else:
                    log.warning("ERROR ERROR 2 fails WTF")

            dataList.append(dataList_p)

        dataListSum = [0,0,0,0]
        for i in range(len(dataList)):
            for j in range(4):
                dataListSum[j] += dataList[i][j]

        #dataListSum = [(1/16)*100*dataListSum[i] for i in range(len(dataListSum))]
        # Plotting
        #self.makeBarPlot(dataListSum,'Percentage of failed trials, (out of n=16)','Percentage of failed trials')
        self.barplot(dataListSum,'Number of failed trials, (out of n=32)','Number of failed trials')

        return dataListSum

    def plotPreference(self):
        dataListSum = [0,0,0,0]

        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)
            dfDict = self.IF.findQ()
            pref = int(dfDict["Value"][1])
            log.debug(name + " Prefers condition: " + str(pref))

            dataListSum[pref] += 1
        self.barplot(dataListSum,'Preferred Condition, n = {}'.format(16),'Number of participants')

        return dataListSum

    def printFinalQ(self,Qnum):  #  A simple function to print out all the answers to 1 of the questions of the final questionnaire
        for name in self.IF.allParticipants:
            self.IF.selectParticipant(name)
            dfDict = self.IF.findQ()
            answer = dfDict["Value"][Qnum]
            log.info("Answer question {}, participant {}: {}".format(Qnum+1,name,answer))

    def checkCollisionsForSteelStrip(self, name):  # check to see if there was collision with the steel strip at any point
        self.IF.selectParticipant(name)

        #for each condition
        for i in range(4):
            logType = "RB"

            dfDict = self.IF.findLogType(i, logType) # contains all logs of a specified trial number

            for file in dfDict:
                infoList0 = self.FI.getInfo(logType, "collidingBody")

                timeList = list(dfDict[file][infoList0["xAxis"]])

                collisionVec = list(dfDict[file][infoList0["yAxis"]])

                blackList = ["Bottom Drip","Top Steel shiny"]  #  Collider names belonging to the steel strip

                for x in range(len(collisionVec)):
                    if collisionVec[x] == blackList[0] or collisionVec[x] == blackList[1]:
                        log.info("{} COLLIDED WITH STEEL STRIP! Condition No: {},Time stamp: {}".format(name,i,timeList[x]))

    ##################### PLOTTING FUNCTIONS #######################
    def barplot(self,dataList,title,ylabel):
        self.getFileName(title)
        df = self.reformatData(dataList)  #  Reformat the data to suit 2way-RM ANOVA

        self.PD.makeBarPlot(df,title,ylabel)

    def boxplot(self,dataList,title,ylabel):
        printtocsv = False
        self.getFileName(title)

        ### Check Assumptions for 2 way RM ANOVA
        rDataList = self.reformatData(dataList)  #  Reformat the data to suit 2way-RM ANOVA

        parametric = self.checkAssumptions(dataList, rDataList)

        anova_table = self.stattests(rDataList)
        if not parametric and self.doTransformation: anova_table = self.stattestsART(anova_table)
        if printtocsv: self.printAnovaToCSV(anova_table)

        self.PD.makeBoxPlot(rDataList,anova_table,title,ylabel)

    def vdLaanPlot(self,usefullness, satisfying, title = ''):
        self.getFileName(title)

        parametric_u = self.checkAssumptions(usefullness)
        parametric_s = self.checkAssumptions(satisfying)

        self.PD.makeVdLaanPlot(usefullness,satisfying,title)

    def reformatData(self,dataList):
        dataType = ""
        df = pd.DataFrame()

        if isinstance(dataList[0],int):  #  So that count metrics can be reformatted using the safe method
            dataList = [dataList]
            dataType = "CountMetrics/"

        for i, dp in enumerate(dataList):
            dftemp  = pd.DataFrame({"Sub_id":[i], "value": [dp[0]], "Display": ["Monitor"],"Cues": ["With"]})
            df = df.append(dftemp)
            dftemp  = pd.DataFrame({"Sub_id":[i], "value": [dp[1]], "Display": ["HMD"],"Cues": ["With"]})
            df = df.append(dftemp)
            dftemp  = pd.DataFrame({"Sub_id":[i], "value": [dp[2]], "Display": ["Monitor"],"Cues": ["Without"]})
            df = df.append(dftemp)
            dftemp  = pd.DataFrame({"Sub_id":[i], "value": [dp[3]], "Display": ["HMD"],"Cues": ["Without"]})
            df = df.append(dftemp)

            # Save the dataframe to a CSV, so reviewers have an easier time analyzing
            df.to_csv(self.analysisPath + dataType + self.outputfilename + ".csv", columns=["Sub_id","Display","Cues","value"],index=False)

        return df

    def checkAssumptions(self, dataList, df = []):  #  Check the assumptions for 2-way repeated measures ANOVA
        shapirocsv = False  #  Save the shapiro wilks results to a single csv, to make it easier to read out

        shapiro_ws = {}
        shapiro_ps = {}
        with open(self.filename, "w") as text_file:  #  To write the results to text file - w mode used here to overwrite the previous file
            # To take the shapiro set over all the data
            #s_w, s_pvalue = stats.shapiro(df['value'])
            #log.debug("Overall Shapiro:\n  w: {}, p: {} \n".format(s_w, s_pvalue))
            #if s_pvalue < 0.05:
                #parametric = False

            for i in range(4):
                shapiroData = [dataList[x][i] for x in range(len(dataList))]  #  Get thedata from the current condition

                w, pvalue = stats.shapiro(shapiroData)  #  Null hypothesis: data is drawn from normal distribution.
                print("Shapiro {}:\n  w: {}, p: {} \n".format(self.labels[i], w, pvalue), file=text_file)
                log.debug("Shapiro {}:\n  w: {}, p: {}".format(self.labels[i], w, pvalue))

                shapiro_ws['dataset'] = "w value: " + self.title
                shapiro_ps['dataset'] = "p value: " +  self.title
                shapiro_ws[self.labels[i]] = np.round(w,3)
                shapiro_ps[self.labels[i]] = np.round(pvalue,3)

                if pvalue < 0.05:
                    parametric = False
                    log.info("ATTENTION:\nData is NOT parametric, performing Aligned Rank Transformation".format(self.labels[i], w, pvalue))
                    print("ATTENTION:\nData is NOT parametric, perform Aligned Rank Transformation".format(self.labels[i], w, pvalue), file=text_file)

                if self.plotqq:
                    self.PD.qqPlot(shapiroData, self.title, self.labels[i], pvalue)

        if shapirocsv:
            doHeader = not os.path.isfile(self.shapirodatapath)
            with open(self.shapirodatapath, "a", newline='') as csvfile:
                fieldnames = ['dataset']
                fieldnames.extend(self.labels)
                wr = csv.DictWriter(csvfile, fieldnames=fieldnames)

                if doHeader: wr.writeheader()

                wr.writerow(shapiro_ws)
                wr.writerow(shapiro_ps)

    def stattests(self, df):
        aovrm2way = AnovaRM(df, 'value', 'Sub_id', within=['Cues', 'Display'])

        res2way = aovrm2way.fit()
        #df_res2way = pd.DataFrame({"Pr > F" : res2way["Pr > F"]})
        anova_table = res2way.anova_table

        #del anova_table["Den DF"]
        #del anova_table["Num DF"]

        # Other libraries (to check if data is calculated correctly)
        #aovrm2way_pg = pg.rm_anova(dv='value', subject='Sub_id', within=['Cues', 'Display'], data=data)
        #print(aovrm2way_pg)

        p_int = anova_table['Pr > F']['Cues:Display']
        if p_int < 0.05 and self.doCheckInteractions:
            if anova_table['Pr > F']['Display'] < 0.05:
                self.checkInteraction(df, 'Display')
            if anova_table['Pr > F']['Cues'] < 0.05:
                self.checkInteraction(df, 'Cues')

        self.printAnova(res2way)
        self.printDFstats(df)

        return anova_table

    def stattestsART(self,anova_table):
        ARTFilePath = self.analysisPath + self.outputfilename + ".art.csv"
        if os.path.isfile(ARTFilePath):
            df = pd.read_csv(ARTFilePath)

            aovrm2way_d = AnovaRM(df, 'ART(value) for Display', 'Sub_id', within=['Cues', 'Display'])
            aovrm2way_c = AnovaRM(df, 'ART(value) for Cues', 'Sub_id', within=['Cues', 'Display'])
            aovrm2way_dxc = AnovaRM(df, 'ART(value) for Display*Cues', 'Sub_id', within=['Cues', 'Display'])

            res2way_d = aovrm2way_d.fit()
            res2way_c = aovrm2way_c.fit()
            res2way_dxc = aovrm2way_dxc.fit()

            # Take the correct values
            res2way = res2way_d
            anova_table = res2way.anova_table

            anova_table['F Value']['Cues'] = res2way_c.anova_table['F Value']['Cues']  #  Display is already correct
            anova_table['Pr > F']['Cues'] = res2way_c.anova_table['Pr > F']['Cues']
            anova_table['F Value']['Cues:Display'] = res2way_dxc.anova_table['F Value']['Cues:Display']
            anova_table['Pr > F']['Cues:Display'] = res2way_dxc.anova_table['Pr > F']['Cues:Display']

            p_int = anova_table['Pr > F']['Cues:Display']
            if p_int < 0.05 and self.doCheckInteractions:
                if anova_table['Pr > F']['Display'] < 0.05:
                    self.checkInteraction(df, 'Display', 'value', False)  #  ART(value) for Display
                if anova_table['Pr > F']['Cues'] < 0.05:
                    self.checkInteraction(df, 'Cues', 'value', False)  #  'ART(value) for Cues'

            #del anova_table["Den DF"]
            #del anova_table["Num DF"]

            with open(self.filename, "a") as text_file:  #  a mode used here to append to shapiro results
                print("ANOVA Results of (ART) Transformed data:", file=text_file)

            self.printAnova(res2way)

        else:
            log.warning('ART DATA HAS NOT YET BEEN GENERATED, PLEASE DO SO USING THE ARTTool in \00 Data Storage\02 Analysis Data\ARToolExe, and run this method again')

        return anova_table

    def printAnova(self,res2way): # Print to text file (and in console)
        anova_table = res2way.anova_table

        with open(self.filename, "a") as text_file:  #  a mode used here to append to shapiro results
            print(res2way, file=text_file)
            log.info(res2way)
            print("\nFull p-values:\n", file=text_file)  #  shows more decimals
            print("p Cues: {}".format(anova_table['Pr > F']['Cues']), file=text_file)
            print("p Display: {}".format(anova_table['Pr > F']['Display']), file=text_file)
            print("p Cues:Display: {}\n".format(anova_table['Pr > F']['Cues:Display']), file=text_file)

    def printAnovaToCSV(self, anova_table):
        doHeader = not os.path.isfile(self.anovadatapath)
        with open(self.anovadatapath, "a", newline='') as csvfile:
            fieldnames = ['dataset','Cues: Pr > F','Cues: F Value','Display: Pr > F','Display: F Value','Cues:Display: Pr > F','Cues:Display: F Value']
            wr = csv.DictWriter(csvfile, fieldnames=fieldnames)
            if doHeader: wr.writeheader()

            anova_csv = {}
            anova_csv['dataset'] = self.title
            anova_csv['Cues: F Value'] = anova_table['F Value']['Cues']
            anova_csv['Cues: Pr > F'] = anova_table['Pr > F']['Cues']
            anova_csv['Display: F Value'] = anova_table['F Value']['Display']
            anova_csv['Display: Pr > F'] = anova_table['Pr > F']['Display']
            anova_csv['Cues:Display: F Value'] = anova_table['F Value']['Cues:Display']
            anova_csv['Cues:Display: Pr > F'] = anova_table['Pr > F']['Cues:Display']

            wr.writerow(anova_csv)

    def printDFstats(self,df): #  Print Mean and Standard deviation to text file
        with open(self.filename, "a") as text_file:  #  a mode used here to append to shapiro results
            df_c = df[df['Cues']=='With']
            df_n = df[df['Cues']=='Without']
            MC = df_c[df_c['Display']=='Monitor']
            HC = df_c[df_c['Display']=='HMD']
            MN = df_n[df_n['Display']=='Monitor']
            HN = df_n[df_n['Display']=='HMD']

            print("MN:\n  Mean : {}\n  Standard Deviation: {}\n".format(np.mean(MN['value']),np.std(MN['value'])), file=text_file)
            print("MC:\n  Mean : {}\n  Standard Deviation: {}\n".format(np.mean(MC['value']),np.std(MC['value'])), file=text_file)
            print("HN:\n  Mean : {}\n  Standard Deviation: {}\n".format(np.mean(HN['value']),np.std(HN['value'])), file=text_file)
            print("HC:\n  Mean : {}\n  Standard Deviation: {}\n".format(np.mean(HC['value']),np.std(HC['value'])), file=text_file)

            print("------\nWhole set:\n  Mean : {}\n  Standard Deviation: {}\n".format(np.mean(df['value']),np.std(df['value'])), file=text_file)

    def checkInteraction(self, df, mainfactor, metric = 'value',parametric = True):
        self.PD.interactionPlot(df, self.title)

        log.info("Cues*Display effect found! Investigating main effect of {}".format(mainfactor))

        ## Cues Pairwise tests:
        # Get the seperate datasets
        df_c = df[df['Cues']=='With']
        df_n = df[df['Cues']=='Without']
        MC = df_c[df_c['Display']=='Monitor']
        HC = df_c[df_c['Display']=='HMD']
        MN = df_n[df_n['Display']=='Monitor']
        HN = df_n[df_n['Display']=='HMD']

        if mainfactor == 'Cues':
            order = [MC,MN,HC,HN]
            labels= ['MC','MN','HC','HN']
        elif mainfactor == 'Display':
            order = [MC,HC,MN,HN]
            labels = ['MC','HC','MN','HN']

        if parametric:
            # perform i t-test to evaluate if both are significant, to check whether there actually is a main effect
            t1_statistic, t1_pvalue = stats.ttest_rel(order[0][metric],order[1][metric])
            t2_statistic, t2_pvalue = stats.ttest_rel(order[2][metric],order[3][metric])
            str1 = "Paired sample t-test, {}-{}:\n  t-statistic: {}, p: {} \n".format(labels[0],labels[1],t1_statistic, t1_pvalue)
            str2 = "Paired sample t-test, {}-{}:\n  t-statistic: {}, p: {} \n".format(labels[2],labels[3],t2_statistic, t2_pvalue)
        else:
            t1_statistic, t1_pvalue = stats.wilcoxon(order[0][metric],order[1][metric])
            t2_statistic, t2_pvalue = stats.wilcoxon(order[2][metric],order[3][metric])
            str1 = "Wilcoxon signed-rank test, {}-{}:\n  t-statistic: {}, p: {} \n".format(labels[0],labels[1],t1_statistic, t1_pvalue)
            str2 = "Wilcoxon signed-rank test, {}-{}:\n  t-statistic: {}, p: {} \n".format(labels[2],labels[3],t2_statistic, t2_pvalue)
        log.info(str1)
        log.info(str2)
        # Write t-test to text file
        with open(self.filename, "a") as text_file:
            print("Cues*Display effect found! Investigating main effect of {}".format(mainfactor), file=text_file)
            print(str1, file=text_file)
            print(str2, file=text_file)

    def tryTransformations(self, data):
        # Try out some transformations
        pvals = []
        wvals = []
        tfs = transformations()
        for l,tf in enumerate(tfs.allTransforms):
            try:  #  Try because sometimes it fails, e.g. due to division by 0
                data_t = [tf(x) for x in data]
                w, pvalue = stats.shapiro(data_t)
                #print("TESTDATA: Shapiro tranformation {}:\n  w: {}, p: {}".format(l+1, w, pvalue))
                pvals.append(pvalue)
                wvals.append(w)
            except:
                pvals.append(np.nan)
                wvals.append(np.nan)
        log.debug("Shapiro test on Transformed data: \n   wvalues = {} \n   pvalues = {}".format(wvals,pvals))

    def checkCorrelations(self,func1,func2):  #  Simple function to check out some correlations
        df1 = self.reformatData(func1)
        df2 = self.reformatData(func2)

        log.info(np.corrcoef(df1['value'], df2['value']))
        self.PD.newfig()
        plt.plot(df1['value'],df2['value'],'o')

    ######### Helper Functions ##########
    def getFileName(self,title):
        self.title = title
        self.outputfilename = title.replace(" ", "")
        self.filename = self.outputPath + "Statistics/" + self.outputfilename + ".txt"  # the path used for saving the figure
        log.info(">>>>>>{}<<<<<<<".format(title))
        return self.filename

    def checkTime(self,timeList):
        finalTime = timeList[-1]
        if finalTime < 295:
            #log.debug("FAILED! at time: " + str(finalTime))
            return False
        else:
            #log.debug("SUCCES!")
            return True

    def findRelData(self, dataList):
        mean = sum(dataList)/len(dataList)
        #print(dataList)
        #print(mean)
        dataList = [p/mean for p in dataList]
        #print(dataList)
        return dataList

    def findStartTime(self,timeList):
        dt = np.diff(timeList)
        negJumps = np.where(dt < 0)[0]

        if negJumps.shape[0] == 0: return 0
        else: return int(negJumps[0])

    def avgData(self,resultList):
        if (len(resultList)>0):
            avg = sum(resultList)/len(resultList)  #  Take average of all trials per participant per condition (2 max :P)
        else:
            log.warning("ERROR ERROR 2 fails WTF, what to do!")  # This never happened, luckily

        return avg

    def dfFromDict(self, usedDict):
        keyList = list(usedDict.keys())
        columns = [str(colname) for colname in keyList]
        df = pd.DataFrame (list(usedDict.values()), columns)
        return df

    def findLogNum(self, logNumString, dfDict):
        for df in dfDict:
            if df.strip(".csv").split("_")[-1] == logNumString:
                return df

    def mag(self,x):
        return math.sqrt(sum(i**2 for i in x))


class transformations(object):
    #https://www.statisticssolutions.com/transforming-data-for-normality/
    def __init__(self):
        self.allTransforms = [self.t1,self.t2,self.t3,self.t4,self.t5,self.t6,self.t7,self.t8,self.t9,self.t10,self.t11,self.t12,self.t13]

    def t1(self,x):
        return 1/np.log10(x)

    def t2(self,x):
        return 1/x

    def t3(self,x):
        return 1/(np.log10(x)**0.5)

    def t4(self,x):
        return 1/(np.log10(x)**(1/3))

    def t5(self,x):
        return 1/(x**2)

    def t6(self,x):
        return 1/(np.log10(x)**2)

    def t7(self,x):
        return np.log10(x)

    def t8(self,x):
        return x**(1/2)

    def t9(self,x):
        return x**(1/4)

    def t10(self,x):
        return x**(1/3)

    def t11(self,x):
        return np.arcsin(x)**(1/2)

    def t12(self,x):
        return np.log10(x/(1-x))

    def t13(self,x):
        return 0.5*np.log10((1+x)/(1-x))


if __name__ == "__main__":  # this part of the script is executed when this file is invoked directly (handy for testing)
    dataMainDir = "./../01 Raw Data/01 Recorded Data/"  #os.getcwd() + "/Data/"
    doSwarm = True
    plotqq = False
    doCheckInteractions = True
    doTransformation = True
    # Note that if data needs to be transformed, it needs to be done using the ARTool in \00 Data Storage\02 Analysis Data\ARToolExe
    # Run this script first, with the transformations disabled, then transform the data, and rerun this script with the transformations enabled
    loggingLevel = logging.INFO  # DEBUG, INFO
    PS = PlotSelection(dataMainDir, doSwarm, plotqq, doTransformation, doCheckInteractions, loggingLevel)
    #PS.printFinalQ(3)

    includeFails = True
    PS.plotNScoops(includeFails)
    PS.plotAvgScoopSize(includeFails)
    PS.plotvdLaan()
    PS.plotScoopFlux(includeFails)
    PS.plotProximity(includeFails)
    PS.plotDross(includeFails)
    PS.plotDisturbance(includeFails)
    PS.plotAvgVel(includeFails)
    PS.plotCollision(includeFails)
    PS.plotAccuracy(includeFails)
    PS.plotFailure()
    PS.plotPreference()

    ### Plot trial data
    #PS.timePlot(1,"obi",["bath"],[1,2],[0])
    #PS.timePlot(1,"RB",["mag"],[1,2],[0])  #  Some of the data still needs some filtering...
    #PS.timePlot(1,"proximity",["mag"],[1,2],[0])
    #PS.timePlot(1,"omni",["mag"],[1,2],[0])

    #PS.checkCorrelations(PS.plotDross(includeFails),PS.plotDisturbance(includeFails))