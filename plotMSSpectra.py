import matplotlib.pyplot as pylab
import glob
import qMS
import vizLib

path = '/home/jhdavis/muspulsePeaks/'

fileList = glob.glob(path+'*.txt')

fileList = qMS.sort_nicely(fileList)
fileList.reverse()
print fileList

pylab.close('all')

vizLib.plotMSSpectra3D(fileList, listOfNames=None, listOfColors=None, gridLines=False, yMin=0.5, yMax=5.5, legend=False, normalizeToN15=False, subtractRef=None, lw=3)
#vizLib.plotMSSpectra3D(fileList, listOfNames=None, listOfColors=None, gridLines=True, yMin=0.5, yMax=8.5, legend=False, normalizeToN15=False, subtractRef=7)

vizLib.plotMSSpectra3D(fileList, listOfNames=None, listOfColors=None, gridLines=False, yMin=0.5, yMax=5.5, legend=False, normalizeToN15=True, lw=3)
pylab.tight_layout()
pylab.show()
