import matplotlib.pyplot as pylab
import qMS
import vizLib
import numpy
import glob

def plotDataSets(files, names, num, den, subunits, title, yLabel, colors, saveName=None, 
                 yMax=1.25, figSize=(22,7), median=False, legendLoc='upper left',
                legendCols=3, normProtein=None, markerSize=None):
                    
    ax = vizLib.makePlotWithFileList(files, num, den, AllProteins=subunits, yMax=yMax, 
                                     names=names, colors=colors, figSize=figSize, 
                                     median=median, normProtein=normProtein, markerSize=markerSize)
                                     
    pylab.legend(loc=legendLoc, ncol=legendCols)
    pylab.xticks(numpy.arange(1,len(subunits)+1,1), [item[4:] for item in subunits], rotation=45)
    ax.set_title(title, multialignment='center')
    ax.set_ylabel(yLabel)
    pylab.tight_layout()
    if not (saveName is None):
        pylab.savefig(saveName)
    return ax

#####################
####Program##########
#####################
if __name__ == '__main__':
    pylab.close('all')
    vizLib.setRcs(scale=12, legendScale=15)
    path = '/home/jhdavis/data/2013_07_27-MSUABCD/handfilt/'
    files = qMS.sort_nicely([i.split('/')[-1] for i in glob.glob(path+'*.csv')])
    
    AllProteins = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 
                   'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 
                   'S18', 'S19', 'S20L26', 'S21', 
                   'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7L12','L9', 'L10', 
                   'L11', 'L13', 'L14', 'L15', 'L16', 'L17', 'L18', 'L19', 
                   'L20', 'L21', 'L22', 'L23', 'L24', 'L25', 'L27', 'L28', 
                   'L29', 'L30', 'L31', 'L32', 'L33', 'L34', 'L35', 'L36']

    LargeSubunit = ['BSubL01', 'BSubL02', 'BSubL03', 'BSubL04', 'BSubL05', 'BSubL06', 'BSubL10', 'BSubL11',
                    'BSubL12', 'BSubL13', 'BSubL14', 'BSubL15', 'BSubL16', 'BSubL17', 'BSubL18', 'BSubL19', 'BSubL20', 'BSubL21',
                    'BSubL22', 'BSubL23', 'BSubL24', 'BSubL27', 'BSubL28', 'BSubL29', 'BSubL30', 'BSubL31a', 'BSubL32',
                    'BSubL33a', 'BSubL35', 'BSubL36']
    
    proteinToNormalizeTo = 'BSubL24'
    reds = ['#fee5d9', '#fc9272', '#de2d26', '#93cae1']
    names = ['msuA', 'msuB', 'mscU', 'msu301_unwashed']
    num = ['AMP_U']
    den = ['AMP_U', 'AMP_S']
    median=False
    yMax = 1.25
    figSize=(16,8)
    print files
    filtPlots_70S_10 = plotDataSets([path+i for i in files], names, num, den, LargeSubunit, '"L16 depletion effect"', '14N/[14N+15N]', reds,
                 saveName=None, yMax=yMax, figSize=figSize, median=median, legendLoc='upper left', legendCols=2, normProtein='BSubL24', markerSize=10)
    pylab.show('all')
    


