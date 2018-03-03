'''
# =============================================================================
#      FileName: data.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-04-24 15:30:17
#       History: YI lI
# =============================================================================
'''

class Segment:

    def __init__(self):
        self.name = ""
        self.chromIdx = -1
        self.chromName = ""
        self.start = -1
        self.end = -1
        self.nReadNum = -1
        self.tReadNum = -1
        self.gc = -1

        self.LOHFrac = -1
        self.LOHStatus = 'NONE'
        self.APMFrac = -1
        self.APMStatus = 'NONE'

        self.pairedCounts = None
        self.BAFCounts = None

        self.baselineLabel = 'FALSE'
        self.tag = '0'
        self.stripeIdx = -1
        self.alleleType = 'NONE'

        self.copyNumber = -1
