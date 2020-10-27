#!/usr/bin/env python

'''
This class is able to return a dictionary containing information useful for
creating plots from the datasets. This information is:
        'title': "Figure Title",
        'xLabel': "Label on x-axis",
        'yLabel': "Label on y-axis",
        'xAxis': "The label of the COLUMN in the CSV files containing the
        data to be plotted on the x-axis",
        'yAxis': "The label of the COLUMN in the CSV files containing the
        data to be plotted on the y-axis"
'''

class FigureInfo(object):
    def __init__(self):
        pass

    def getInfo(self,logType,plotSpecification='empty'):
        switcher = {
            'ivolume': self.ivolume,
            'obi': self.obi,
            'omni': self.omni,
            'prox': self.proximity,
            'RB': self.RB
        }
        if logType in switcher:
            func = switcher.get(logType, lambda: "invalid logType name: {}".format(logType))
            if plotSpecification == 'empty':
                self.info = func()
                return func()
            else:
                self.info = func(plotSpecification)
                return func(plotSpecification)

        else:
            print("invalid logType name: {}".format(logType))

    def ivolume(self, plotType = 'vol'):
        plotInfoDict = {}

        PID_vol = {'title': "Sumberged volume scoop",
                    'xLabel': "Time [s]",
                    'yLabel': "Volume Fraction",
                    'xAxis': "Time Stamp",
                    'yAxis': "Intersecting Volume"
        }
        PID_flux = {'title': "Volumetric flux scoop",
                    'xLabel': "Time [s]",
                    'yLabel': "Volumetric flux [1/s]",
                    'xAxis': "Time Stamp",
                    'yAxis': "Volume Flux"
        }
        if plotType == 'vol':
            plotInfoDict = PID_vol
        elif plotType == 'flux':
            plotInfoDict = PID_flux
        elif plotType == 'all':
            plotInfoDict = {
                'vol' : PID_vol,
                'flux' : PID_flux
                }
        else:
            print("plotType not found")
        return plotInfoDict

    def obi(self, plotType = 'bath'):
        plotInfoDict = {}
        PID_bath = {'title': "Particles in Bath",
                    'xLabel': "Time [s]",
                    'yLabel': "Number of partles",
                    'xAxis': "Time Stamp",
                    'yAxis': "Particles in Bath"
        }
        PID_bin1 = {'title': "Particles in Bin 1",
                    'xLabel': "Time [s]",
                    'yLabel': "Number of partles",
                    'xAxis': "Time Stamp",
                    'yAxis': "Particles in Bin 1"
        }
        PID_bin2 = {'title': "Particles in Bin 2",
                    'xLabel': "Time [s]",
                    'yLabel': "Number of partles",
                    'xAxis': "Time Stamp",
                    'yAxis': "Particles in Bin 2"
        }
        if plotType == 'bath':
            plotInfoDict = PID_bath
        elif plotType == 'bin1':
            plotInfoDict = PID_bin1
        elif plotType == 'bin2':
            plotInfoDict = PID_bin2
        elif plotType == 'all':
            plotInfoDict = {
                'bath' : PID_bath,
                'bin1' : PID_bin1,
                'bin2' : PID_bin2
                }
        else:
            print("plotType not found")
        return plotInfoDict

    '''def questionnaire(self, plotType = 'bath'):
        plotInfoDict = {}
        PID_bath = {'title': "Particles in Bath",
                    'xLabel': "Time [s]",
                    'yLabel': "Number of partles",
                    'xAxis': "Time Stamp",
                    'yAxis': "Particles in Bath"
        }
        PID_bin1 = {'title': "Particles in Bin 1",
                    'xLabel': "Time [s]",
                    'yLabel': "Number of partles",
                    'xAxis': "Time Stamp",
                    'yAxis': "Particles in Bin 1"
        }
        PID_bin2 = {'title': "Particles in Bin 2",
                    'xLabel': "Time [s]",
                    'yLabel': "Number of partles",
                    'xAxis': "Time Stamp",
                    'yAxis': "Particles in Bin 2"
        }
        if plotType == 'v':
            plotInfoDict = PID_bath
        elif plotType == 'bin1':
            plotInfoDict = PID_bin1
        elif plotType == 'bin2':
            plotInfoDict = PID_bin2
        elif plotType == 'all':
            plotInfoDict = {
                'bath' : PID_bath,
                'bin1' : PID_bin1,
                'bin2' : PID_bin2
                }
        else:
            print("plotType not found")
        return plotInfoDict'''

    def omni(self, plotType = 'mag'): # could add the anchorPosition etc. not sure how useful it is though
        plotInfoDict = {}
        PID_force = {'title': "Omni Stylus Force",
                    'xLabel': "Time [s]",
                    'yLabel': "Force [N]",
                    'xAxis': "Time Stamp",
                    'yAxis': ["Stylus Force X","Stylus Force Y","Stylus Force Z"],
                    'legend': ["Force X","Force Y","Force Z"]
        }
        PID_mag = {'title': "Omni Stylus Force Magnitude",
                    'xLabel': "Time [s]",
                    'yLabel': "Force [N]",
                    'xAxis': "Time Stamp",
                    'yAxis': "Stylus Force Magnitude"
        }
        PID_pos = {'title': "Omni Stylus Position",
                    'xLabel': "Time [s]",
                    'yLabel': "Position [mm]",
                    'xAxis': "Time Stamp",
                    'yAxis': ["stylusPositionRaw X","stylusPositionRaw Y","stylusPositionRaw Z"],
                    'legend': ["X position","Y position","Z position"]
        }
        PID_vel = {'title': "Omni Stylus Velocity",
                    'xLabel': "Time [s]",
                    'yLabel': "Velocity [mm/s]",
                    'xAxis': "Time Stamp",
                    'yAxis': ["stylusVelocityRaw X","stylusVelocityRaw Y","stylusVelocityRaw Z"],
                    'legend': ["X velocity","Y velocity","Z velocity"]
        }
        PID_rot = {'title': "Omni Stylus Rotation",
                    'xLabel': "Time [s]",
                    'yLabel': "Rotation [deg]",
                    'xAxis': "Time Stamp",
                    'yAxis': ["stylusRotationWorld X","stylusRotationWorld Y","stylusRotationWorld Z"],
                    'legend': ["X rotation","Y rotation","Z rotation"]
        }
        if plotType == 'force':
            plotInfoDict = PID_force
        elif plotType == 'mag':
            plotInfoDict = PID_mag
        elif plotType == 'pos':
            plotInfoDict = PID_pos
        elif plotType == 'vel':
            plotInfoDict = PID_vel
        elif plotType == 'rot':
            plotInfoDict = PID_rot
        elif plotType == 'all':
            plotInfoDict = {
                'force' : PID_force,
                'pos' : PID_pos,
                'vel' : PID_vel,
                #'rot' : PID_rot
                }
        else:
            print("plotType not found")
        return plotInfoDict

    def proximity(self, plotType = 'mag'): # could add the anchorPosition etc. not sure how useful it is though
        plotInfoDict = {}
        PID_mag = {'title': "Distance to environment",
                    'xLabel': "Time [s]",
                    'yLabel': "Distance [m]",
                    'xAxis': "Time Stamp",
                    'yAxis': "Magnitude"
        }
        PID_pos = {'title': "Distance to environment",
                    'xLabel': "Time [s]",
                    'yLabel': "Distance [m]",
                    'xAxis': "Time Stamp",
                    'yAxis': ["Distance X","Distance Y","Distance Z"],
                    'legend': ["X distance","Y distance","Z distance"]
        }
        if plotType == 'mag':
            plotInfoDict = PID_mag
        elif plotType == 'pos':
            plotInfoDict = PID_pos
        elif plotType == 'all':
            plotInfoDict = {
                'mag' : PID_mag,
                'pos' : PID_pos
                }
        else:
            print("plotType not found")
        return plotInfoDict

    def RB(self, plotType = 'force'): # could do the titles based on the object name but since scoop is the only object I log... might add later if i would decide to log other RBs
        plotInfoDict = {}
        PID_collidingBody = {'title': "Colliding Body",
                    'xLabel': "Time [s]",
                    'yLabel': "Colliding Body",
                    'xAxis': "Time Stamp",
                    'yAxis': "Colliding Body",
        }
        PID_forceMag = {'title': "Scoop Collision Force Magnitude",
                    'xLabel': "Time [s]",
                    'yLabel': "Scoop Collision Force [N]",
                    'xAxis': "Time Stamp",
                    'yAxis': "Impact Force Magnitude",
        }
        PID_force = {'title': "Scoop Collision Force",
                    'xLabel': "Time [s]",
                    'yLabel': "Scoop Collision Force [N]",
                    'xAxis': "Time Stamp",
                    'yAxis': ["Impact Force X","Impact Force Y","Impact Force Z"],
                    'legend': ["Collision Force X","Collision Force Y","Collision Force Z"]
        }
        PID_pos = {'title': "Scoop Position",
                    'xLabel': "Time [s]",
                    'yLabel': "Position [mm]",
                    'xAxis': "Time Stamp",
                    'yAxis': ["Position X","Position Y","Position Z"],
                    'legend': ["X position","Y position","Z position"]
        }
        PID_vel = {'title': "Scoop Velocity",
                    'xLabel': "Time [s]",
                    'yLabel': "Velocity [mm/s]",
                    'xAxis': "Time Stamp",
                    'yAxis': ["LinVelocity X","LinVelocity Y","LinVelocity Z"],
                    'legend': ["X velocity","Y velocity","Z velocity"]
        }
        PID_rot = {'title': "Scoop Rotation",
                    'xLabel': "Time [s]",
                    'yLabel': "Rotation [deg]",
                    'xAxis': "Time Stamp",
                    'yAxis': ["Rotation X","Rotation Y","Rotation Z"],
                    'legend': ["X rotation","Y rotation","Z rotation"]
        }
        PID_angvel = {'title': "Scoop Angular Velocity",
                    'xLabel': "Time [s]",
                    'yLabel': "Angular Velocity [deg/s]", # <<<<??? units
                    'xAxis': "Time Stamp",
                    'yAxis': ["AngVelocity X","AngVelocity Y","AngVelocity Z"],
                    'legend': ["X Angular Velocity","Y Angular Velocity","Z Angular Velocity"]
        }
        if plotType == 'mag':
            plotInfoDict = PID_forceMag
        elif plotType == 'force':
            plotInfoDict = PID_force
        elif plotType == 'pos':
            plotInfoDict = PID_pos
        elif plotType == 'vel':
            plotInfoDict = PID_vel
        elif plotType == 'rot':
            plotInfoDict = PID_rot
        elif plotType == 'angvel':
            plotInfoDict = PID_angvel
        elif plotType == 'collidingBody':
            plotInfoDict = PID_collidingBody
        elif plotType == 'all':
            plotInfoDict = {
                'mag' : PID_forceMag,
                #'force' : PID_force,
                'pos' : PID_pos,
                #'vel' : PID_vel,
                #'rot' : PID_rot,
                #'angvel' : PID_angvel
                }
        else:
            print("plotType not found")

        return plotInfoDict