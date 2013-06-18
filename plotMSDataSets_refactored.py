import matplotlib.pyplot as pylab
import qMS
import vizLib

#####################
####Program##########
#####################
if __name__ == '__main__':
    pylab.close('all')
    pulse = False
    numerator = ["ampu"]
    denominator = ["ampu", "ampl"]
    path = '/home/jhdavis/scripts/python/figures/plotMSDataSets_mcMasterMSU/data/depletion'
    
    AllProteins = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 
                   'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 
                   'S18', 'S19', 'S20L26', 'S21', 
                   'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7L12','L9', 'L10', 
                   'L11', 'L13', 'L14', 'L15', 'L16', 'L17', 'L18', 'L19', 
                   'L20', 'L21', 'L22', 'L23', 'L24', 'L25', 'L27', 'L28', 
                   'L29', 'L30', 'L31', 'L32', 'L33', 'L34', 'L35', 'L36']
    
    proteinToNormalizeTo = 'S4'
    '''L14 = [path+i for i in ["/dep1_iso.csv", "/dep2_iso.csv", "/dep3_iso.csv"]]
    stats = []
    ax = vizLib.makePlotWithFileList(L14, AllProteins, pulse, numerator, denominator, proteinToNormalizeTo='S4', yMax=1.5)
    ax.set_title('L14')
    pylab.tight_layout()
    pylab.savefig('L14.png')
    
    
    proteinToNormalizeTo = 'L22'
    S17 = [path+i for i in ["/dep4_iso.csv", "/dep5_iso.csv"]]
    stats = []
    ax = vizLib.makePlotWithFileList(S17, AllProteins, pulse, numerator, denominator, proteinToNormalizeTo='S4', yMax=1.5)
    ax.set_title('S17')
    pylab.tight_layout()
    pylab.savefig('S17.png')
    
    proteinToNormalizeTo = 'S4'
    L17 = [path+i for i in ["/dep6_iso.csv", "/dep7_iso.csv", "/dep8_iso.csv"]]
    stats = []
    ax = vizLib.makePlotWithFileList(L17, AllProteins, pulse, numerator, denominator, proteinToNormalizeTo='S4', yMax=2.5)
    ax.set_title('L17')
    pylab.tight_layout()
    pylab.savefig('L17.png')
    
    proteinToNormalizeTo = 'S4'
    L10 = [path+i for i in ["/dep9_iso.csv", "/dep10_iso.csv", "/dep11_iso.csv"]]
    stats = []
    ax = vizLib.makePlotWithFileList(L10, AllProteins, pulse, numerator, denominator, proteinToNormalizeTo='S4', yMax=1.5)
    ax.set_title('L10')
    pylab.tight_layout()
    pylab.savefig('L10.png')
    
    proteinToNormalizeTo = 'S4'
    L14 = [path+i for i in ["/dep12_iso.csv", "/dep13_iso.csv", "/dep14_iso.csv"]]
    stats = []
    ax = vizLib.makePlotWithFileList(L14, AllProteins, pulse, numerator, denominator, proteinToNormalizeTo='S4', yMax=1.5)
    ax.set_title('L14-2')
    pylab.tight_layout()
    pylab.savefig('L14-2.png')
    
    proteinToNormalizeTo = 'S4'
    L28 = [path+i for i in ["/dep15_iso.csv", "/dep16_iso.csv", "/dep17_iso.csv"]]
    stats = []
    ax = vizLib.makePlotWithFileList(L28, AllProteins, pulse, numerator, denominator, proteinToNormalizeTo='S4', yMax=1.5)
    ax.set_title('L28')
    pylab.tight_layout()
    pylab.savefig('L28.png')
    
    proteinToNormalizeTo = 'S4'
    L22 = [path+i for i in ["/dep18_iso.csv", "/dep19_iso.csv"]]
    stats = []
    ax = vizLib.makePlotWithFileList(L22, AllProteins, pulse, numerator, denominator, proteinToNormalizeTo='S4', yMax=1.5)
    ax.set_title('L22')
    pylab.tight_layout()
    pylab.savefig('L22.png')
    '''
    L17_5 = qMS.getInfoToPlotStats(path+"/dep6_iso.csv", pulse, numerator, denominator, proteinToNormalizeTo = 'S4')
    L17_p4 = qMS.getInfoToPlotStats(path+"/dep7_iso.csv", pulse, numerator, denominator, proteinToNormalizeTo = 'S4')
    L17_0 = qMS.getInfoToPlotStats(path+"/dep8_iso.csv", pulse, numerator, denominator, proteinToNormalizeTo = 'S4')    
    
    L14_5 = qMS.getInfoToPlotStats(path+"/dep1_iso.csv", pulse, numerator, denominator, proteinToNormalizeTo = 'S4')
    L14_p4 = qMS.getInfoToPlotStats(path+"/dep2_iso.csv", pulse, numerator, denominator, proteinToNormalizeTo = 'S4')
    L14_0 = qMS.getInfoToPlotStats(path+"/dep3_iso.csv", pulse, numerator, denominator, proteinToNormalizeTo = 'S4')    
    
    L28_5 = qMS.getInfoToPlotStats(path+"/dep15_iso.csv", pulse, numerator, denominator, proteinToNormalizeTo = 'S4')
    L28_p4 = qMS.getInfoToPlotStats(path+"/dep16_iso.csv", pulse, numerator, denominator, proteinToNormalizeTo = 'S4')
    L28_0 = qMS.getInfoToPlotStats(path+"/dep17_iso.csv", pulse, numerator, denominator, proteinToNormalizeTo = 'S4')    
    
    L17list = [L17_5, L17_p4, L17_0]
    L14list = [L14_5, L14_p4, L14_0]
    L28list = [L28_5, L28_p4, L28_0]
    
    L17 = [qMS.getMedian(i, 'L17') for i in L17list]
    L14 = [qMS.getMedian(i, 'L14') for i in L14list]
    L28 = [qMS.getMedian(i, 'L28') for i in L28list]
    
    xs = [5.0, 0.4, 0.0]
    f = pylab.figure()
    ax = f.add_subplot(111)
    ax.plot(xs, L17, marker='o', c='r', markersize=12, label='L17')
    ax.plot(xs, L14, marker='o', c='b', markersize=12, label='L14')
    ax.plot(xs, L28, marker='o', c='g', markersize=12, label='L28')
    ax.set_xlim((-1,6))
    ax.set_ylim((0.5,2.5))
    ax.legend(loc='upper left')
    ax.set_title('Hsl dependent strains')
    ax.set_xlabel('Hsl [nM]')
    ax.set_ylabel('median protein abundance (normalized to S4)')
    pylab.tight_layout()
    pylab.savefig('scatter.png')
    pylab.show('all')
    
    


