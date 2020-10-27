#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import pingouin as pg

from matplotlib.ticker import MaxNLocator
from statsmodels.graphics.factorplots import interaction_plot

import logging
log = logging.getLogger()
'''
This class outputs figures based on the data that is provided to it.
It serves more as a tool that can be used for the making of figures containing
plots, based on the df (dataframes) that are given as input.
'''

class PlotData(object):
    def __init__(self, outputPath, doSwarm):
        # Settings:
        self.titlefontsize = 15
        self.save = True
        self.trimTime = True                    #  Set this to false to leave on the leading bit of the data (overflow from previous repetition)

        # Class Params
        self.doSwarm = doSwarm
        self.outputPath = outputPath            #  Path where the figures will be saved to
        self.labels = ["MC","HC","MN","HN"]     #  labels used for plotting, M=Monitor,H=HMD;C=Cues,N=No cues
        self.figCount = 1                       #  To keep figures from overlapping

        # Prettier Colors  :   https://matplotlib.org/3.1.1/gallery/style_sheets/style_sheets_reference.html
        plt.style.use('ggplot')
        self.colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # get the colors of this theme as I think they look nicer
        plt.style.use('default')

    ########### Time Plots ###########

    def timePlot(self, ax, df, dfXColumn, dfYColumns, legend = []):
        # dfXColumn is the header of the column we want to plot on the X axis
        # dfYColumns is a list of headers of all the columns we want to plot on the y axis

        if isinstance(dfYColumns,str): # if its not a list, make it a list so it doesnt error out
            dfYColumns = [dfYColumns]

        xValues = df[dfXColumn]

        # Trim of the lead data
        if self.trimTime: startElement = self.findStartTime(xValues)+1  #  used to cut off the leading part of the data (for some reason part of the previous dataset flows over to the next one)
        else: startElement = 0

        plt.sca(ax)  # set the current axes

        for dfYColumn in dfYColumns:
            yValues = df[dfYColumn]
            plt.plot(xValues[startElement:], yValues[startElement:], label = dfYColumn)

        return plt

    def timeFigure(self,df,figureInfo, fig=False,ax=False, supTitle = ""):
        if fig==None and ax==None:
            fig,ax = plt.subplots()

        # get the info
        title = figureInfo.get('title')
        legend = figureInfo.get('legend')
        X = figureInfo.get('xAxis')
        Ys = figureInfo.get('yAxis')
        xlabel = figureInfo.get('xLabel')
        ylabel = figureInfo.get('yLabel')

        # set in figure
        ax.set_title(title)
        self.timePlot(ax, df, X, Ys)

        if len(supTitle) > 0:
            fig.suptitle(supTitle)

        currhandles, currlabels = ax.get_legend_handles_labels()
        if isinstance(legend, list):
            ax.legend(loc='best', handles = currhandles, labels = legend)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid()
        #ax.autoscale(enable=True)
        return fig,ax

    def addTimeSubPlot(self,df,figureInfo,fig,ax = None):
        if ax == None:
            ax = fig.axes[0]

            geometry=ax.get_geometry()
            geometry1 = []
            geometry2 = []
            nrows = []
            ncols = []

            tupleAdd = lambda t1,t2: tuple(np.array(t1)+np.array(t2))
            geometry1 = tupleAdd((1,0,0),geometry)
            nrows =geometry1[0]
            ncols =geometry1[1]

            geometry2 = tupleAdd((0,0,1),geometry1)

            ax.change_geometry(*geometry1)
            ax2 = fig.add_subplot(*geometry2)
        else:
            ax2 = ax

        # set on gridspec
        fig,ax2 = self.timeFigure(df,figureInfo,fig,ax2)
        #allaxes = fig.get_axes()
        #print(len(allaxes))
        plt.subplots_adjust(top = 2, right = 2, hspace = 0.5)

        return fig,ax2

    ########### Bar Plots/Boxplots ###########
    def makeBarPlot(self, df, title = '', ylabel = '', save = True):
        # Set class variables
        self.title = title

        ### Initialize
        self.newfig() # Start a new figure

        ax = sns.barplot(x="Cues",y='value',hue="Display",data=df, palette="Set1")
        fig = ax.get_figure()

        # Figure Formatting
        ax.set_title(title, fontsize = self.titlefontsize)
        #ax.set_xlabel('Condition')
        ax.set_ylabel(ylabel)
        #ax.yaxis.set_major_locator(MaxNLocator(integer=True))  #  only show ints on y-axis

        if self.save:
            self.saveFig(title, fig)

        return fig , ax

    def makeBoxPlot(self, rDataList, anova_table, title = '', ylabel = ''):
        # Settings
        a = 0.05  #  criteria p < a for significance.

        self.maxvalue = rDataList['value'].max()

        # Set class variables
        self.title = title

        fig = self.newfig() # Start a new figure

        if self.doSwarm: #  Make the scatter plots
            swarmplt = sns.swarmplot(x="Cues", y="value",hue="Display", size=6, marker="D",linewidth=1,edgecolors="black", dodge=True, data=rDataList, palette="Set1")

        ax = sns.boxplot(x="Cues", y="value", hue="Display",dodge=True, data = rDataList, palette="Set1")

        ax,legendBbox = self.formatBoxPlot(ax,fig)
        self.annotateboxplot(fig, anova_table, legendBbox, a)  #  Annotate the significant p values into the plot

        # Figure Formatting
        ax.set_title(title, fontsize = self.titlefontsize)
        #fig.subplots_adjust(top=0.80)
        ax.set_ylabel(ylabel)

        if self.save:
            self.saveFig(title, fig)

        return fig , ax

    def formatBoxPlot(self, ax, fig = None):
        if fig == None: fig = self.fig
        ### Format Plot
        #ax.set_aspect(self.get_aspect(ax)*1.5) # Set aspect ratio so it looks cleaner
        # Set Correct Legend
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend(handles[0:2],["Monitor","HMD"], loc="upper left", bbox_to_anchor=(1.1, 0.6), borderaxespad=0.)  # Set custom legend
        legendBbox = leg.get_tightbbox(fig.canvas.get_renderer())
        return ax, legendBbox

    def makeVdLaanPlot(self, usefullness, satisfying, title = ''):
        useBoxPlotPercentiles = True  #  Use box plots percentiles for the whiskers, instead of the mean+variance

        # Settings
        capsize = 0.02  # of the whiskers
        plotlim = 2.1
        markersize = 4
        linewidth = 1.5

        dashspacing = 8

        #axLimits = [-.1,2,-.1,2]
        axLimits = [-2,2,-2,2]
        ticks = [-2,-1,0,1,2]
        #ticks = [-0.5,0,0.5,1,1.5,2]

        # Set class variables
        self.title = title

        ### Plot Data
        # Throw some boxplots in a new plot to get the values of the whiskers
        tempfig = self.newfig()
        # Start a new figure
        vdlfig = self.newfig(None)
        #colors = [self.colors[0],self.colors[1],self.colors[0],self.colors[1]]
        colors = ['royalblue','darkred','deepskyblue','indianred']
        alphas = [1,1,1,1]
        markers = ["^","^","s","s"]

        plt.figure(vdlfig.number)
        # x,y=0 middle lines
        plt.plot([-plotlim,plotlim],[0,0],"k--", linewidth=linewidth/2,dashes=(dashspacing, dashspacing))
        plt.plot([0,0],[-plotlim,plotlim],"k--", linewidth=linewidth/2,dashes=(dashspacing, dashspacing))

        if self.doSwarm: #  Make the scatter plots (these are plotted in a separate loop, so that the other information is plotted on top and isnt obtsructed)
            for i in range(4):
                sati = [satisfying[x][i] for x in range(len(satisfying))]
                usei = [usefullness[x][i] for x in range(len(usefullness))]
                swarmplt, = plt.plot(sati,usei,linewidth=0,marker=markers[i],mfc = "w", color = colors[i], alpha = alphas[i], markersize = markersize)

        for i in range(4):
            sati = [satisfying[x][i] for x in range(len(satisfying))]
            usei = [usefullness[x][i] for x in range(len(usefullness))]
            plt.figure(vdlfig.number)

            if useBoxPlotPercentiles:
                midLabel = 'Median'
                # Get the boxplot data
                plotmetric = 'quartile'  #  'quartile' or 'whisker'
                labels = ["sat","use"]

                plt.figure(tempfig.number)
                bp = plt.boxplot([sati,usei],labels=labels)
                bpdata = self.get_box_plot_data(bp,labels)
                plt.figure(vdlfig.number)

                # Plot Medians
                plt.plot(bpdata['median'][0],bpdata['median'][1],marker=markers[i], color = colors[i], alpha = alphas[i])
                # Plot middle line
                mid, = plt.plot([bpdata['lower_{}'.format(plotmetric)][0],bpdata['upper_{}'.format(plotmetric)][0]],
                         [bpdata['median'][1],bpdata['median'][1]], color = colors[i], linewidth=linewidth)
                plt.plot([bpdata['median'][0],bpdata['median'][0]],
                         [bpdata['lower_{}'.format(plotmetric)][1],bpdata['upper_{}'.format(plotmetric)][1]], color = colors[i],linewidth=linewidth, alpha = alphas[i])

                # Plot Whiskers
                plt.plot([bpdata['median'][0]-capsize,bpdata['median'][0]+capsize],
                         [bpdata['lower_{}'.format(plotmetric)][1],bpdata['lower_{}'.format(plotmetric)][1]], color = colors[i],linewidth=linewidth, alpha = alphas[i])
                plt.plot([bpdata['median'][0]-capsize,bpdata['median'][0]+capsize],
                         [bpdata['upper_{}'.format(plotmetric)][1],bpdata['upper_{}'.format(plotmetric)][1]], color = colors[i],linewidth=linewidth, alpha = alphas[i])

                plt.plot([bpdata['lower_{}'.format(plotmetric)][0],bpdata['lower_{}'.format(plotmetric)][0]],
                         [bpdata['median'][1]-capsize,bpdata['median'][1]+capsize], color = colors[i],linewidth=linewidth, alpha = alphas[i])
                plt.plot([bpdata['upper_{}'.format(plotmetric)][0],bpdata['upper_{}'.format(plotmetric)][0]],
                         [bpdata['median'][1]-capsize,bpdata['median'][1]+capsize], color = colors[i],linewidth=linewidth, alpha = alphas[i])
            else:
                midLabel = 'Mean'  # used for the legend
                # Plot Means
                means = [np.mean(sati),np.mean(usei)]
                variances = [np.std(sati),np.std(usei)]
                plt.plot(means[0],means[1], marker=markers[i],markersize = markersize * 2, color = colors[i])

                # Plot Variances
                mid, = plt.plot([means[0]+variances[0],means[0]-variances[0]],
                         [means[1],means[1]], color = colors[i], linewidth=linewidth, alpha = alphas[i])
                plt.plot([means[0],means[0]],
                         [means[1]+variances[1],means[1]-variances[1]], color = colors[i],linewidth=linewidth, alpha = alphas[i])

                # Plot Whiskers
                plt.plot([means[0]-capsize,means[0]+capsize],
                         [means[1]-variances[1],means[1]-variances[1]], color = colors[i],linewidth=linewidth, alpha = alphas[i])
                plt.plot([means[0]-capsize,means[0]+capsize],
                         [means[1]+variances[1],means[1]+variances[1]], color = colors[i],linewidth=linewidth, alpha = alphas[i])

                plt.plot([means[0]-variances[0],means[0]-variances[0]],
                         [means[1]-capsize,means[1]+capsize], color = colors[i],linewidth=linewidth, alpha = alphas[i])
                plt.plot([means[0]+variances[0],means[0]+variances[0]],
                         [means[1]-capsize,means[1]+capsize], color = colors[i],linewidth=linewidth, alpha = alphas[i])

        # Get the axes object
        ax = vdlfig.axes[0]

        ### Format Plot
        ax.set_xlim(axLimits[0],axLimits[1])
        ax.set_ylim(axLimits[2],axLimits[3])
        plt.xticks(ticks)
        plt.yticks(ticks)
        ax.set_aspect(1)

        ax.set_title(title, fontsize = self.titlefontsize)

        #vdlfig.suptitle(title)
        ax.set_xlabel('Satisfying')
        ax.set_ylabel('Usefulness')

        # Legend
# Old version, kept here for reference
# =============================================================================
#         legend_elements = [plt.Line2D([0], [0], color=colors[0], lw=3, label='Monitor'),
#                            plt.Line2D([0], [0], color=colors[1], lw=3, label='HMD'),
#                            plt.Line2D([0], [0], lw = 0, mfc = "w", marker=markers[0], color='k', label='With Cues', markersize=6),
#                            plt.Line2D([0], [0], lw = 0, mfc = "w", marker=markers[2], color='k', label='Without Cues', markersize=6),
#                            plt.Line2D([0], [0], lw = 0, marker="o", color='k', label=midLabel, markersize=6)]
# =============================================================================
        legendlabels = ["Mon+Cues","HMD+Cues","Mon, No cues","HMD, No cues"]
        legend_elements = [plt.Line2D([0], [0], color=colors[i], marker=markers[i], lw=2, label=legendlabels[i]) for i in range(4)]
        legend_elements = self.swapPositions(legend_elements,1,2)
        leg = plt.legend(handles=legend_elements)
        legendBbox = leg.get_tightbbox(vdlfig.canvas.get_renderer())
        #self.annotateboxplot(bpfig, anova_table_u, legendBbox_u, a)  #  Annotate the significant p values into the plot
        if self.save:
            self.saveFig(title, vdlfig)

        return vdlfig , ax

    def swapPositions(self, list, pos1, pos2):
        list[pos1], list[pos2] = list[pos2], list[pos1]
        return list

    def interactionPlot(self,df, title, valuecol = 'value'):
        fig = self.newfig()
        ax = plt.gca()
        fig = interaction_plot(df['Cues'],df['Display'], df[valuecol], ax = ax, markers=['^','s'], ms=10)
        ax.set_title(title)

        fig.savefig(self.outputPath + "Figures/Interactions/" + title + "_interactionplot" + ".png", bbox_inches='tight',dpi=300)

        return ax

    def qqPlot(self, idata, title, label, pvalue):
        fig = self.newfig()
        ax = plt.gca()
        pg.qqplot(idata, dist='norm', ax = ax)
        ititle = title + ", qq-plot: " + label
        ax.set_title(ititle)
        plt.text(0.8, 0.2,"p = {}".format(np.round(pvalue,5)), ha='center', va='center', transform = ax.transAxes)

        fig.savefig(self.outputPath + "Figures/" + str(title) + ".png", bbox_inches='tight',dpi=300)

        return ax

    ########### ANNOTATIONS ###########
    def annotateboxplot(self, fig, anova_table, legendBbox = [], a = 0.05):  #  Annotate p-values into the boxplots
        # Get the correct figure
        plt.figure(fig.number)
        ax = plt.gca()

        pkey = "Pr > F"
        p1 = anova_table[pkey]["Cues"]
        p2 = anova_table[pkey]["Display"]
        p3 = anova_table[pkey]["Cues:Display"]

        if p1 < a:
            self.topannotation(0, 1, p1, ax)
        if p2 < a:
            self.legendannotation(legendBbox, p2, ax)#-.25, .25, maxval, p2)
        if p3 < a:  # Turns out this wasnt necessary, but left it in the code for later reference
            log.debug("Cues*Display effect found! Perform interction analysis")
            '''t1_statistic, t1_pvalue = stats.ttest_rel([dataList[x][0] for x in range(len(dataList))],[dataList[x][2] for x in range(len(dataList))])
            log.info("Paired sample t-test, Cues:\n  t-statistic: {}, p: {} \n".format(t1_statistic, t1_pvalue))
            t2_statistic, t2_pvalue = stats.ttest_rel([dataList[x][1] for x in range(len(dataList))],[dataList[x][3] for x in range(len(dataList))])
            log.info("Paired sample t-test, Display:\n  t-statistic: {}, p: {} \n".format(t2_statistic, t2_pvalue))
            # Write t-test to text file
            with open(self.filename, "a") as text_file:
                print("Paired sample t-test, Cues:\n  t-statistic: {}, p: {} \n".format(t1_statistic, t1_pvalue), file=text_file)
                print("Paired sample t-test, Display:\n  t-statistic: {}, p: {} \n".format(t2_statistic, t2_pvalue), file=text_file)'''


    def sigdots(self,p):
        #● BLACK CIRCLE        25CF
        #⚫ MEDIUM BLACK CIRCLE 26AB
        #⬤ BLACK LARGE CIRCLE  2B24
        dotstr = ""
        if p<0.05:
            dotstr += '●'
        if p<0.01:
            dotstr += '●'
        if p<0.001:
            dotstr += '●'

        return dotstr

    def topannotation(self,x1,x2, p, ax):
        maxval = self.maxvalue
        ylim = ax.get_ylim()
        ax.set_ylim(ylim[0],ylim[1]*1.3) # Taller to fit the p-values above to it

        h = (maxval/10)
        y = maxval + 2*h
        col = 'k'
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*.5, y+1.2*h, self.sigdots(p), ha='center', va='center', color=col)  #  "p = {}".format(np.round(p,5))

    def legendannotation(self, legendBbox, p, ax):
        arrowfraction = 0.3

        #width   = legendBbox.x1-legendBbox.x0
        #height  = legendBbox.y1-legendBbox.y0
        arrowx  = 0.78#legendBbox.x0-2.5*width
        arrowy0 = 0.51#legendBbox.y0+height/8
        arrowy1 = 0.59#legendBbox.y1-height/12
        arrowheight = arrowy1-arrowy0
        arrowwidth  = arrowheight * arrowfraction

        ann     = ax.annotate('', xy=(arrowx, arrowy0), xycoords='figure fraction',
                  xytext=(arrowx, arrowy1), textcoords='figure fraction',
                  annotation_clip=False,
                  arrowprops=dict(arrowstyle="-",
                                  connectionstyle="bar, fraction={}".format(arrowfraction),
                                  ec="k"))
        anndot  = ax.annotate(self.sigdots(p), xy=(arrowx-arrowwidth, arrowy0+arrowheight*.55), xycoords='figure fraction',
                  textcoords='figure fraction',
                  ha="center",
                  va="center",
                  rotation=90,
                  annotation_clip=False)

    ########### Helper Functions ###########
    def findStartTime(self,timeList):  # This function is the same as in the PlotSelection script, but placed here for convenience
        dt = np.diff(timeList)
        negJumps = np.where(dt < 0)[0]

        if negJumps.shape[0] == 0: return 0
        else: return int(negJumps[0])

    def newfig(self, figsize = (4,5) ): # Start a new figure
        self.fig = plt.figure(self.figCount, figsize=figsize)
        self.figCount += 1
        return self.fig

    def saveFig(self, title, fig = None):
        if fig == None: fig = self.fig
        else:
            title = title.replace(" ","").replace(":","_").replace(",","")
        fig.savefig(self.outputPath + "Figures/" + str(title) + ".png", bbox_inches='tight',dpi=300)

    def get_box_plot_data(self, bp, labels):
        rows_list = []

        for i in range(len(labels)):
            dict1 = {}
            dict1['label'] = labels[i]
            dict1['lower_whisker'] = bp['whiskers'][i*2].get_ydata()[1]
            dict1['lower_quartile'] = bp['boxes'][i].get_ydata()[1]
            dict1['median'] = bp['medians'][i].get_ydata()[1]
            dict1['upper_quartile'] = bp['boxes'][i].get_ydata()[2]
            dict1['upper_whisker'] = bp['whiskers'][(i*2)+1].get_ydata()[1]
            rows_list.append(dict1)

        return pd.DataFrame(rows_list)

    def get_aspect(self, ax=None):
        if ax is None:
            ax = plt.gca()
        fig = ax.figure

        ll, ur = ax.get_position() * fig.get_size_inches()
        width, height = ur - ll
        axes_ratio = height / width
        aspect = axes_ratio / ax.get_data_ratio()

        return aspect