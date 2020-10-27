#!/usr/bin/env python
import pandas as pd
import os
from os import listdir


"""
This class finds all the data in a given filepath to the data folder, and is able
to create dictionaries of either all logs of 1 filetype or all logs of 1 trial no.
based on the participant name or number that is given.

These dictionaries consists of:
    key     = path to the particular file
    value   = dataframes (df) that is basically a dictionary containing all the
    lists of data in a column (the value), speficied by its header (the key)

The participant can be changed without having to create another instance of this class=
"""

class ImportFiles(object):
    def __init__(self, basePath, participant = 0): # find all the files in the basepath
        self.basePath = basePath # the path containing all the data divided in folders per participant

        # Get a dict of participants & their data folders
        allParticipantNames = listdir(self.basePath) # Find all the folders containing the data for each participant
        allParticipantDirs = [self.basePath + str(Name) for Name in allParticipantNames] # add the full path to the path string
        self.allParticipants = dict(zip(allParticipantNames, allParticipantDirs)) # put them in a dict
        self.selectParticipant(participant)

        # PlayerInfo is found by just going to participantPath / PlayerInfo.csv
        # FinalQuestionnaire is found by just going to participantPath / FinalQuestionnaire.csv


    def selectParticipant(self, participant):
        # Pick the given participant and set the class variables accordingly
        if (type(participant) == int): # get the requested data from the participant number (in order of the directory list found by listdir)
            self.participantPath = list(self.allParticipants.values())[participant]
        elif (type(participant) == str): # get the requested data from the participant name
            self.participantPath = self.allParticipants[participant]

        allFileNames = listdir(self.participantPath) # Find all the folders containing the data for each participant
        allFileDirs = [self.participantPath +"/"+ str(Name) for Name in allFileNames] # add the full path to the path string

        self.allFiles = dict(zip(allFileNames, allFileDirs)) # put them in a dict --> This isnt really useful as it only contains the folders with the conditions/


    def _findFiles(self, fileList, fileExtension): # find all the files of a certain extension
        files = [] # make a new list for this file extension

        for file in fileList:
            filestr = str(file)
            if filestr.endswith(fileExtension):
                files.append(file) # add the dict item to the new dict

        if len(files) == 1: # what is this used for? ...
            files = str(files[0])

        return files

    def findQ(self):
        fileList = self._findFiles(listdir(self.participantPath),"csv")  # based on the condition number (which is the name of the folder) -> Right now this is actually just a list (not a dict), but that should work the same

        if isinstance(fileList, list):
            for file in fileList:
                if (file.startswith("FinalQuestionnaire")):
                    dfDict = pd.read_csv(self.participantPath +"/"+ str(file))
        else:
            dfDict = pd.read_csv(self.participantPath +"/"+ fileList)

        return dfDict

    def findvdLaan(self, condNo):
        condPath = self.participantPath + "/" + str(condNo)
        fileList = self._findFiles(listdir(condPath),"csv")  # based on the condition number (which is the name of the folder) -> Right now this is actually just a list (not a dict), but that should work the same

        if isinstance(fileList, list):
            for file in fileList:
                if (file.startswith("Questionnaire")):
                    dfDict = pd.read_csv(condPath +"/"+ str(file))
        return dfDict

    def findLog(self, condNo, logType, logNum): # put data of all trials for 1 type in a dict & return it
        if isinstance(logNum, int):
            logNum = [logNum]

        # possible logTypes: ivolume, obi, omni, proximity, RB, Questionnaire
        condPath = self.participantPath + "/" + str(condNo)
        fileList = self._findFiles(listdir(condPath),"csv")  # based on the condition number (which is the name of the folder) -> Right now this is actually just a list (not a dict), but that should work the same

        dfDict = {}

        if isinstance(fileList, list):
            for file in fileList:
                # Check the log number
                correctLogNum = False
                for no in logNum:
                    if file.strip(".csv").endswith(str(no)):
                        correctLogNum = True

                if file.startswith(logType) and correctLogNum:
                    dfDict[file] = pd.read_csv(condPath +"/"+ str(file))

        return dfDict

    def findLogType(self, condNo, logType): # put data of all trials for 1 type in a dict & return it
        # possible logTypes: ivolume, obi, omni, proximity, RB, Questionnaire
        condPath = self.participantPath + "/" + str(condNo)
        fileList = self._findFiles(listdir(condPath),"csv")  # based on the condition number (which is the name of the folder) -> Right now this is actually just a list (not a dict), but that should work the same

        dfDict = {}

        if isinstance(fileList, list):
            for file in fileList:
                if (file.startswith(logType)):
                    dfDict[file] = pd.read_csv(condPath +"/"+ str(file))

        return dfDict

    def findLogNum(self, condNo, logNum): # put data of all logtypes of 1 trial in a dict & return it
        # logNum is an int of the trial number for the log
        condPath = self.participantPath + "/" + str(condNo)
        fileList = self._findFiles(listdir(condPath),"csv")  # based on the condition number (which is the name of the folder)

        dfDict = {}

        if isinstance(fileList, list):
            for file in fileList:
                if (file.strip(".csv").endswith(str(logNum))):
                    dfDict[file] = pd.read_csv(condPath +"/"+ str(file))

        return dfDict




if __name__ == "__main__":  # this part of the script is executed when this file is invoked directly (handy for testing)
    dataMainDir = "./Data/"  #os.getcwd() + "/Data/"
    IF = ImportFiles(dataMainDir, 0)

    typePandas = IF.findLogType(0,"proximity")
    #print([typePandas[item].head() for item in typePandas])
    itemList = list(typePandas.values())[0]
    print(itemList["Time Stamp"])

    numPandas = IF.findLogNum(0,0)
    #print([typePandas[item].head() for item in typePandas])