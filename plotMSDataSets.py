import csv
import numpy
import matplotlib.pyplot as pylab

def readCsv(datapath, pulse):
    data = list( csv.reader( open(datapath, 'rU') ) )
    header = data[0]

# 1 - Find the indices for the quantities of interest using list comprehensions

    protein_index = [index for index, item in enumerate(header) if item == "protein"][0]
    ampu_index = [index for index, item in enumerate(header) if item == "AMP_U"][0]
    ampl_index = [index for index, item in enumerate(header) if item == "AMP_L"][0]
    UID_index = [index for index, item in enumerate(header) if item == "isofile"][0]
    if pulse:
        amps_index = [index for index, item in enumerate(header) if item == "AMP_S"][0]

# 2 - Declare the peptide list and the protein set

    peptide_list = []
    protein_set  = set()
    peptide_set = set()

    data_dictionary = {}
    csv_hold_dict = {}
    UID_list = []

# 3 - Loop over data set, collect amplitudes, charge state, peptide sequence, protein id into protein_set

    for line in data[1:]:
        protein = line[protein_index]
        ampu = float(line[ampu_index])
        ampl = float(line[ampl_index])
        if pulse:
            amps = float(line[amps_index])
        else:
            amps = float(0)
        UID = line[UID_index]
        csv_hold_dict[UID] = line
        UID_list.append(UID)
        identifier = {"protein": protein, "ampu": ampu, "ampl": ampl, "amps": amps, "uid": UID}
	#######################
        data_dictionary[UID] = {"ampu": ampu, "ampl": ampl, "amps": amps, "protein": protein}
	#######################
        protein_set.add(protein)
        peptide_list.append(identifier)
    return [data_dictionary, protein_set, peptide_list, peptide_set]

def outputStatsFile(protein_set, dataByProtein, filePrefix, numerator, denominator):
    fileSuffix = "-"
    for i in numerator:
        fileSuffix = fileSuffix+str(i)+"_p_"
    fileSuffix = fileSuffix[:-3]+"_o_"

    for i in denominator:
        fileSuffix = fileSuffix+str(i)+"_p_"
    fileName = filePrefix+fileSuffix[:-3]+".stats"

    outfile = open(fileName, 'w')
    outfile.write("Protein flab +\\- loc nval vals\n")

    ordered = list(protein_set)
    ordered.sort()
    locIter = 1
    for p in ordered:
        outString = str(p) + " " + str(dataByProtein[p]["flab"])[:6] + " "
        outString = outString + str(dataByProtein[p]["loc"])[:6] + " " + str(locIter) + " "
        outString = outString + str(dataByProtein[p]["nvals"]) + " "
            
        for v in numpy.sort(dataByProtein[p]["vals"]):
            outString = outString + str(v)[:6] + " "
        outString = outString[:-1] + "\n"
        outfile.write(outString)
    outfile.close()

def addStatsToPlot(proteins, dataByProtein, ax, name, offset=0.0, markerSize=40, color='#377db8'):
    proteins.sort()
    xAxis = range(1,len(proteins)+1)
    xs = []
    ys = []
    for x in xAxis:
        p = proteins[x-1]
        if p in dataByProtein.keys():
            for v in dataByProtein[p]["vals"]:
                xs.append(x+offset)
                ys.append(v)
    ax.plot(xs, ys, 'o', color=color, markersize=markerSize, label=name)
    return ax

def plotStatsFile(proteins, dataByProtein, name, offset=0.0, markerSize=40, color='#e31a1c'):
    proteins.sort()
    xAxis = range(1,len(proteins)+1)
    yMax = 1.50
    #fig = pylab.figure(figsize=(22,5))
    fig = pylab.figure(figsize=(22,3))
    ax = fig.add_subplot(111)
    xs = []
    ys = []
    for x in xAxis:
        p = proteins[x-1]
        if p in dataByProtein.keys():
            for v in dataByProtein[p]["vals"]:
                xs.append(x+offset)
                ys.append(v)
    ax.plot(xs, ys, 'o', color=color, markersize=markerSize, label=name)
    pylab.xticks(xAxis, [item[4:] for item in proteins], rotation=45, size=15)
    pylab.xlim(1, len(proteins)+1)
    pylab.yticks([0,yMax/4, 2*yMax/4, 3*yMax/4,yMax], size=15)
    pylab.ylim(0, yMax)
    pylab.grid(b=True, which='major', color='grey', linestyle='--', axis='y', linewidth=1.5, alpha=0.5)
    pylab.grid(b=True, which='major', color='grey', linestyle='-', axis='x', linewidth=1.5, alpha=0.75)
    return ax
        
def getInfoToPlotStats(datapath, pulse, numerator, denominator, proteinToNormalizeTo):
    filePrefix = datapath[:-3]
    [data_dictionary, protein_set, peptide_list, peptide_set] = readCsv(datapath, pulse)
    dataByProtein = calculateStatsFile(data_dictionary, protein_set, numerator, denominator, normalization=1.0)
    normMed = determineNormalization(dataByProtein, proteinToNormalizeTo)
    normalizedDataByProtein = calculateStatsFile(data_dictionary, protein_set, numerator, denominator, normalization=normMed)
    #outputStatsFile(protein_set, dataByProtein, filePrefix, numerator, denominator)
    return normalizedDataByProtein

def composeDataSets(comp, toAdd):
    proteinsToAdd = toAdd.keys()
    for i in proteinsToAdd:
        if i in comp.keys():
            comp[i]["vals"] = comp[i]["vals"] + toAdd[i]["vals"]
            comp[i]["nvals"] = len(comp[i]["vals"])
            comp[i]["loc"] = numpy.std(comp[i]["vals"])
            comp[i]["flab"] = numpy.mean(comp[i]["vals"])
        else:
            comp[i] = toAdd[i]
    return comp

def mergeFiles(fileList, pulse, numerator, denominator, proteinToNormalizeTo, LargeSubunit):
    normalizedDataByProtein = getInfoToPlotStats(fileList[0], pulse, numerator, denominator, proteinToNormalizeTo)
    compositeDataSet = normalizedDataByProtein
    for datapath in fileList[1:]:
        normalizedDataByProtein = getInfoToPlotStats(datapath, pulse, numerator, denominator, proteinToNormalizeTo)
        compsiteDataSet = composeDataSets(compositeDataSet, normalizedDataByProtein)
    return compositeDataSet
                    
def makePlotWithFileList(fileList, pulse, numerator, denominator, proteinToNormalizeTo, LargeSubunit):
    offsets = float(len(fileList)+1)
    #reds = ['#ae2221', '#c32625', '#d72c2b', '#db4140', '#df5655', '#e36c6b', '#e78180']
    #reds = ['#25557d', '#3170a4', '#377db8', '#5696cc']
    reds = ['#ae2221', '#d72c2b', '#df5655', '#e78180']
    normalizedDataByProtein = getInfoToPlotStats(fileList[0], pulse, numerator, denominator, proteinToNormalizeTo)
    ax = plotStatsFile(LargeSubunit, normalizedDataByProtein, fileList[0], offset=(1.0/offsets), markerSize = 10/float(len(fileList))+7, color=reds[0])
    i = 1
    for datapath in fileList[1:]:
        normalizedDataByProtein = getInfoToPlotStats(datapath, pulse, numerator, denominator, proteinToNormalizeTo)
        ax = addStatsToPlot(LargeSubunit, normalizedDataByProtein, ax, datapath, offset=(1.0/offsets)*(i+1), markerSize = 10/float(len(fileList))+7, color=reds[i])
        i = i +1
    #pylab.legend(loc='lower left')
    return ax

def makePlotWithDataSets(listOfDataSets, LargeSubunit, name):
    set1 = listOfDataSets[0]
    offsets = float(len(listOfDataSets)+1)
    #ax2 = plotStatsFile(LargeSubunit, set1, name[0], color = '#e31a1c', offset=(.75/offsets), markerSize = 10/float(len(listOfDataSets))+4)
    #ax2 = plotStatsFile(LargeSubunit, set1, name[0], offset=0.3, markerSize = 12)
    ax2 = plotStatsFile(LargeSubunit, set1, name[0], offset=0.5, markerSize = 15)
    i = 1
    for dataSet in listOfDataSets[1:]:
        #ax2 = addStatsToPlot(LargeSubunit, dataSet, ax2, name[i], offset=(1.0/offsets)*(i+1.25), markerSize = 10/float(len(listOfDataSets))+4)
        ax2 = addStatsToPlot(LargeSubunit, dataSet, ax2, name[i], offset=.7, markerSize = 12)
        i = i +1
    return ax2


#####################
####Program##########
#####################
if __name__ == '__main__':
    pylab.close('all')
    pulse = False
    numerator = ["ampu"]
    denominator = ["ampu", "ampl"]
    path = '/home/jhdavis/scripts/python/figures/plotMSDataSets_mcMasterMSU/McMasterMSUCsvs/'
    proteinToNormalizeTo = "BSubL24"
    '''LargeSubunit = ['BSubL01', 'BSubL02', 'BSubL03', 'BSubL04', 'BSubL05', 'BSubL06', 'BSubL07', 'BSubL09', 'BSubL10', 'BSubL11',
                    'BSubL12', 'BSubL13', 'BSubL14', 'BSubL15', 'BSubL16', 'BSubL17', 'BSubL18', 'BSubL19', 'BSubL20', 'BSubL21',
                    'BSubL22', 'BSubL23', 'BSubL24', 'BSubL27', 'BSubL28', 'BSubL29', 'BSubL30', 'BSubL31a', 'BSubL31b', 'BSubL32',
                    'BSubL33a', 'BSubL33b', 'BSubL34', 'BSubL35', 'BSubL36']
    '''
    LargeSubunit = ['BSubL01', 'BSubL02', 'BSubL03', 'BSubL04', 'BSubL05', 'BSubL06', 'BSubL10', 'BSubL11',
                    'BSubL12', 'BSubL13', 'BSubL14', 'BSubL15', 'BSubL16', 'BSubL17', 'BSubL18', 'BSubL19', 'BSubL20', 'BSubL21',
                    'BSubL22', 'BSubL23', 'BSubL24', 'BSubL27', 'BSubL28', 'BSubL29', 'BSubL30', 'BSubL32',
                    'BSubL33a', 'BSubL34', 'BSubL35', 'BSubL36']
    


    McMaster45S = [path+i for i in ["McMaster45S_esi-run1_filt.csv", "McMaster45S_esi-run2.1_filt.csv", "McMaster45S_esi-run2_filt.csv", "McMaster45S_qtof_filt_filtppm.csv"]]
    McMaster50S = [path+i for i in ["McMaster50S_esi-run1_filt.csv", "McMaster50S_esi-run2.1_filt.csv", "McMaster50S_esi-run2_filt.csv", "McMaster50S_qtof_filt_filtppm.csv"]]
    #MSU_lowSalt = ["MSU_45S_lowSalt_filt.csv", "MSU_45S_lowSalt2_filt.csv"]
    #MSU_lowSalt = ["MSU_45S_lowSalt_filt.csv"]
    #MSU_highSalt = ["MSU_45S_highSalt_filt.csv", "MSU301-45S-esi-run2_filt.csv"]
    #MSU_highSalt = ["MSU_45S_highSalt_filt.csv"]
    #MSU_pellet = ["MSU_45S_highSaltPellet_filt.csv"]
    #MSU_50S = ["MSU_50S_filt.csv", "MSU418-50S-esi-run2_filt.csv"]
    
    
    
    
    #fileLists = [MSU_lowSalt, MSU_highSalt, MSU_pellet, McMaster45S, McMaster50S]
    fileLists = [McMaster45S]
    merged = []
    for listOfFiles in fileLists:
        makePlotWithFileList(listOfFiles, pulse, numerator, denominator, proteinToNormalizeTo, LargeSubunit)
        merged.append(mergeFiles(listOfFiles, pulse, numerator, denominator, proteinToNormalizeTo, LargeSubunit))
    #JHD70S = ["JHD70SPulse_esi_run1_filt.csv", "JHD70S-qtof_filt.csv"]
    #makePlotWithFileList(JHD70S, True, ["ampu", "ampl"], ["ampu", "ampl", "amps"], proteinToNormalizeTo, LargeSubunit)
    #JHD70S_merged = mergeFiles(JHD70S, True, ["ampu", "ampl"], ["ampu", "ampl", "amps"], proteinToNormalizeTo, LargeSubunit)
    #merged.append(JHD70S_merged)
    makePlotWithDataSets(merged, LargeSubunit, ["McMaster45S_merged", "McMaster50S_merged"])
    pylab.show()

